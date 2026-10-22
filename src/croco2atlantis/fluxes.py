"""Water transport between boxes, from the C-grid velocities of the model.

The transport is taken on the faces of the model cells, where `u` and `v` live: every
cell face whose two cells belong to different boxes (see `horizontal`) carries its
velocity times its length and its water thickness from one box to the other. These
are the transports the model itself uses, so the net inflow of a box is the
divergence of the model over the cells of that box, consistent with the volumes of
the box means. Sampling the velocity along the faces of the BGM instead does not
close the budget of the boxes.

Box -1 stands for outside the boxes: the open boundaries of the model (the water
cells of the ghost ring) and any water cell left without a box. Its exchanges are the
water that enters or leaves the boxes through the boundary of the model.

Layers here are the standard ones for every face, independent of the boxes: the first
starts at the free surface and the last extends to the bottom. Mapping the layers of a
face to the layers of its two boxes (a box with fewer layers receives the deeper ones
in its deepest layer) is done when the exchanges are built. Layer index 0 is the
shallowest, as in `vertical`.

With average files the transport is the mean velocity times the thickness from the
mean free surface, which leaves out the correlation between the tidal current and
the tidal elevation (the tidal Stokes transport).
"""

from dataclasses import dataclass

import numpy as np

from . import vertical


@dataclass
class CellFaces:
    """Model cell faces that separate two boxes, and the box pairs they connect.

    `pairs` is (n_pair, 2), each pair once with pairs[:, 0] < pairs[:, 1]; a positive
    transport goes from pairs[:, 0] to pairs[:, 1]. The other arrays have one entry
    per cell face, the u faces first: `index_u` and `index_v` are the flat indices of
    the faces in the u (eta_rho, xi_rho - 1) and v (eta_rho - 1, xi_rho) arrays,
    `cells` (2, n) the flat rho indices of the cells before and after the face (along
    xi for u, eta for v), `length` the length of the face in metres, `pair` its pair
    and `sign` +1 when a positive velocity goes from pairs[pair, 0] to pairs[pair, 1].
    """

    pairs: np.ndarray
    nbox: int
    index_u: np.ndarray
    index_v: np.ndarray
    cells: np.ndarray
    length: np.ndarray
    pair: np.ndarray
    sign: np.ndarray


def cell_faces(cellmap, grid):
    """Find the cell faces between different boxes for a `CellBoxMap` and its `CrocoGrid`.

    A face counts when both of its cells are water (interior or open boundary) and at
    least one is interior. Open boundary cells belong to box -1.
    """
    inner = np.asarray(grid.mask, dtype=bool)
    wet = inner if grid.open_boundary is None else inner | grid.open_boundary
    box = np.where(inner, cellmap.box, -1)
    flat = np.arange(box.size).reshape(box.shape)

    index, cells, length, sides = [], [], [], []
    for axis, metric in ((1, grid.pn), (0, grid.pm)):  # u faces (along xi), then v faces
        before = tuple(slice(None, -1) if a == axis else slice(None) for a in range(2))
        after = tuple(slice(1, None) if a == axis else slice(None) for a in range(2))
        sel = wet[before] & wet[after] & (inner[before] | inner[after]) & (box[before] != box[after])
        index.append(np.flatnonzero(sel))
        cells.append(np.stack([flat[before][sel], flat[after][sel]]))
        length.append((2.0 / (metric[before] + metric[after]))[sel])
        sides.append(np.stack([box[before][sel], box[after][sel]]))

    sides = np.concatenate(sides, axis=1)
    pairs, pair = np.unique(np.sort(sides, axis=0).T, axis=0, return_inverse=True)
    return CellFaces(
        pairs=pairs.reshape(-1, 2),
        nbox=cellmap.nbox,
        index_u=index[0],
        index_v=index[1],
        cells=np.concatenate(cells, axis=1),
        length=np.concatenate(length),
        pair=pair.reshape(-1),
        sign=np.where(sides[0] < sides[1], 1.0, -1.0),
    )


def standard_layer_bounds(bottoms):
    """Top and bottom depth of the layers for a face: surface to bottom, no truncation."""
    bottoms = np.asarray(bottoms, dtype=float)
    top = vertical.layer_tops(bottoms)
    bot = bottoms.copy()
    top[0] = -np.inf
    bot[-1] = np.inf
    return top, bot


def face_transport(faces, u, v, z_w, bottoms):
    """Transport between every pair of boxes in every layer for one time step, m3/s.

    `u` is (n_s, eta_rho, xi_rho - 1), `v` (n_s, eta_rho - 1, xi_rho) and `z_w` the
    level interfaces (n_s + 1, eta, xi), negative downwards, index 0 at the bottom.
    The levels on a face are the mean of those of its two cells, as in the model.
    Returns (n_pair, n_layer), positive from pairs[:, 0] to pairs[:, 1].
    """
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    zw = np.asarray(z_w, dtype=float)
    zw = zw.reshape(zw.shape[0], -1)
    velocity = np.concatenate(
        [u.reshape(u.shape[0], -1)[:, faces.index_u], v.reshape(v.shape[0], -1)[:, faces.index_v]], axis=1
    )  # (n_s, n_faces)
    z_face = 0.5 * (zw[:, faces.cells[0]] + zw[:, faces.cells[1]])

    top, bot = standard_layer_bounds(bottoms)
    ov = vertical.overlap(z_face, top, bot)  # (n_layer, n_s, n_faces)
    per_face = (ov * velocity[None]).sum(axis=1) * (faces.length * faces.sign)  # (n_layer, n_faces)

    out = np.zeros((len(faces.pairs), len(bottoms)))
    for k in range(len(bottoms)):
        out[:, k] = np.bincount(faces.pair, weights=per_face[k], minlength=len(faces.pairs))
    return out


def net_inflow(faces, transport):
    """Net horizontal inflow of every box in every layer, (n_box, n_layer), m3/s.

    `transport` is the output of `face_transport`. Exchanges with box -1 (the model
    boundary) count as inflow or outflow of the box on the other side.
    """
    out = np.zeros((faces.nbox + 1, transport.shape[1]))  # the extra last row takes box -1
    np.add.at(out, faces.pairs[:, 1], transport)
    np.add.at(out, faces.pairs[:, 0], -transport)
    return out[:-1]


def pairs_without_face(faces, bgm):
    """Pairs of boxes (both >= 0) that exchange water but share no face in the BGM.

    Cells are assigned to boxes by their centre, so two boxes that only touch at a
    corner, or are separated by a thin box, can end up with neighbouring cells.
    """
    shared = {tuple(p) for p in np.sort(bgm.face_lr(), axis=1).tolist()}
    missing = [p for p in faces.pairs.tolist() if p[0] >= 0 and tuple(p) not in shared]
    return np.array(missing, dtype=int).reshape(-1, 2)
