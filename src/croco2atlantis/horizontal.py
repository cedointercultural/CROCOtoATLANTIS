"""Horizontal mapping from model cells to ATLANTIS boxes.

Every water cell is assigned to the box that contains its centre. Water cells whose
centre falls outside every box or inside a box without water (an island of the BGM
that the model grid does not resolve) are moved to the nearest box with water, so
that no model water is lost. The value of a box in a layer is the mean over its
cells weighted by cell area and by the thickness of the cell's water inside the
layer, i.e. a volume-weighted mean. Box-layers that end up without data (a box with
no cell centre, or a layer deeper than any of its columns) are filled with the
profile of the nearest box that has data in all of them.
"""

from dataclasses import dataclass

import numpy as np
import shapely

from . import vertical


@dataclass
class CellBoxMap:
    """Box of every model cell. `box` is (eta, xi), -1 for land and cells outside all boxes.

    `moved` marks the water cells that were moved to the nearest box with water (see
    `move_to_water_boxes`).
    """

    box: np.ndarray
    nbox: int
    moved: np.ndarray | None = None

    def __post_init__(self):
        flat = self.box.ravel()
        order = np.argsort(flat, kind="stable")
        sorted_box = flat[order]
        edges = np.searchsorted(sorted_box, np.arange(-1, self.nbox + 1))
        # columns[b] = flat indices of the cells of box b
        self.columns = [order[edges[b + 1] : edges[b + 2]] for b in range(self.nbox)]

    @property
    def counts(self):
        """Number of cells in every box."""
        return np.array([len(c) for c in self.columns])


def assign_cells(x, y, mask, polygons):
    """Assign each water cell to the polygon that contains its centre.

    `x`, `y` are the cell centres in the projection of the polygons, shaped like
    `mask`. A centre on the edge between two polygons, or inside overlapping
    polygons, goes to the first one in the list.
    """
    box = np.full(np.shape(mask), -1, dtype=int)
    free = np.asarray(mask, dtype=bool).copy()
    for b, poly in enumerate(polygons):
        inside = np.zeros_like(free)
        inside[free] = shapely.intersects_xy(poly, np.asarray(x)[free], np.asarray(y)[free])
        box[inside] = b
        free &= ~inside
    return box


def move_to_water_boxes(box, x, y, mask, polygons, water, max_distance=None):
    """Move water cells outside every box or inside a box without water to the nearest box with water.

    `box` is the output of `assign_cells` and `water` a boolean array over the
    polygons, True for the boxes that have water (layers). The distance is from the
    cell centre to the polygon. With `max_distance` (metres) the cells farther than
    that from every box with water are left with box -1. Returns (box, moved), the
    second marking the cells that were moved.
    """
    box = np.array(box)
    mask = np.asarray(mask, dtype=bool)
    water = np.asarray(water, dtype=bool)
    lost = mask & ((box < 0) | ~water[np.maximum(box, 0)])
    moved = np.zeros(box.shape, dtype=bool)
    if not lost.any():
        return box, moved
    target = np.flatnonzero(water)
    tree = shapely.STRtree([polygons[b] for b in target])
    points = shapely.points(np.asarray(x)[lost], np.asarray(y)[lost])
    (point, nearest), dist = tree.query_nearest(points, return_distance=True, all_matches=False)
    new = np.full(len(points), -1)
    near = np.ones(len(point), bool) if max_distance is None else dist <= max_distance
    new[point[near]] = target[nearest[near]]
    box[lost] = new
    moved[lost] = new >= 0
    return box, moved


def assign_cells_from_grid(grid, bgm, max_distance=None):
    """Assign the water cells of a `CrocoGrid` to the boxes of a `Bgm`.

    Cells outside every box or inside an island are moved to the nearest box with
    water (see `move_to_water_boxes`, which also explains `max_distance`).
    """
    x, y = bgm.project(grid.lon, grid.lat)
    polygons = bgm.box_polygons()
    box = assign_cells(x, y, grid.mask, polygons)
    box, moved = move_to_water_boxes(box, x, y, grid.mask, polygons, bgm.depth > 0, max_distance)
    return CellBoxMap(box, bgm.nbox, moved)


class BoxLayerWeights:
    """Volume weights of every model cell in every layer of its box, for one time step.

    Build it with `box_layer_weights`; it depends on `zeta`, so it is rebuilt at each
    time step and reused for every variable of that step.
    """

    def __init__(self, entries, nbox, nlayer):
        self.entries = entries  # per box: None or (cell flat indices, weights (L, n_s, ncells))
        self.nbox = nbox
        self.nlayer = nlayer

    @property
    def volume(self):
        """Volume of model water in every box and layer, m3 (0 where there is none)."""
        out = np.zeros((self.nbox, self.nlayer))
        for b, e in enumerate(self.entries):
            if e is not None:
                out[b] = e[1].sum(axis=(1, 2))
        return out

    def apply(self, values):
        """Volume-weighted mean of `values` (n_s, eta, xi) in every box and layer.

        Returns (n_box, n_layer), NaN where the box has no water in the layer. A
        missing value (NaN) in a water cell is an error, not a gap to skip.
        """
        values = np.asarray(values, dtype=float)
        flat = values.reshape(values.shape[0], -1)
        out = np.full((self.nbox, self.nlayer), np.nan)
        for b, e in enumerate(self.entries):
            if e is None:
                continue
            cols, w = e
            v = flat[:, cols]
            if not np.isfinite(v).all():
                raise ValueError(f"missing values in the water cells of box {b}")
            den = w.sum(axis=(1, 2))
            num = (w * v[None]).sum(axis=(1, 2))
            ok = den > 0
            out[b, ok] = num[ok] / den[ok]
        return out


def box_layer_weights(z_w, cellmap, cell_area, depth, bottoms, deepest_to_bottom=True):
    """Weights for averaging model cells into box layers at one time step.

    `z_w` are the level interfaces (n_s + 1, eta, xi), negative downwards, index 0 at
    the bottom. `depth` is the mean depth of every box (positive). See `vertical` for
    the layer rule.
    """
    top, bot = vertical.box_layer_bounds(depth, bottoms, deepest_to_bottom)
    count = vertical.n_layers(depth, bottoms)
    zw = np.asarray(z_w).reshape(z_w.shape[0], -1)
    area = np.asarray(cell_area).ravel()
    entries = []
    for b in range(cellmap.nbox):
        cols = cellmap.columns[b]
        if len(cols) == 0 or count[b] == 0:
            entries.append(None)
            continue
        ov = vertical.overlap(zw[:, cols], top[b], bot[b])
        entries.append((cols, ov * area[cols]))
    return BoxLayerWeights(entries, cellmap.nbox, len(bottoms))


def fill_missing(means, depth, bottoms, centroids, max_distance=None):
    """Fill box-layers without data with the profile of the nearest box that has them.

    `means` is (n_box, n_layer) with NaN where there is no data; only box-layers that
    the box really has (see `vertical.n_layers`) are filled. All the missing layers of
    a box come from one source, the nearest box (by `centroids`, (n_box, 2) in
    metres) with data in every one of them, so that a filled column is a real
    profile and not layers stitched from different boxes. Layers with no data in any
    box cannot be filled and stay NaN. With `max_distance` (metres) a box whose
    source is farther than that stays NaN, so boxes outside the model domain are not
    given values copied from far away. Returns (filled, was_filled), the second a
    boolean array.
    """
    filled = means.copy()
    was_filled = np.zeros(means.shape, dtype=bool)
    have = np.isfinite(means)
    has_layer = np.arange(means.shape[1])[None, :] < vertical.n_layers(depth, bottoms)[:, None]
    need = has_layer & ~have & have.any(axis=0)
    for b in np.flatnonzero(need.any(axis=1)):
        src = np.flatnonzero(have[:, need[b]].all(axis=1))
        if len(src) == 0:
            continue
        dist = np.linalg.norm(centroids[src] - centroids[b], axis=1)
        i = dist.argmin()
        if max_distance is not None and dist[i] > max_distance:
            continue
        filled[b, need[b]] = means[src[i], need[b]]
        was_filled[b, need[b]] = True
    return filled, was_filled
