"""Vertical remapping from model levels to the fixed ATLANTIS layers.

Depths here are positive downwards. Layers are given by their bottom depths, for
example (18, 54, 270, 684, 918, 1782, 3100). Each box keeps the layers whose top is
above its mean depth, and the deepest kept layer is truncated at that depth. For the
integration, the deepest layer of each box can instead take the whole water column
below its top (`deepest_to_bottom`), because the mean depth of a box is an average
and some model columns are deeper.

Layer index 0 is the shallowest here; the ATLANTIS order (0 = bottom of each box) is
applied later, when the files are written.
"""

import numpy as np


def layer_tops(bottoms):
    """Top depth of every layer for the given layer bottoms."""
    bottoms = np.asarray(bottoms, dtype=float)
    return np.concatenate([[0.0], bottoms[:-1]])


def n_layers(depth, bottoms):
    """Number of layers of each box: the layers whose top is above the mean depth.

    `depth` is the mean depth of each box in metres, positive. Islands (depth 0) get
    no layers.
    """
    tops = layer_tops(bottoms)
    depth = np.asarray(depth, dtype=float)
    return (tops[None, :] < depth[:, None]).sum(axis=1)


def nominal_thickness(depth, bottoms):
    """Thickness of every layer of every box, truncated at the mean depth.

    Returns an array (n_box, n_layer) with 0 where the box has no such layer.
    """
    bottoms = np.asarray(bottoms, dtype=float)
    tops = layer_tops(bottoms)
    depth = np.asarray(depth, dtype=float)
    bot = np.minimum(bottoms[None, :], depth[:, None])
    return np.clip(bot - tops[None, :], 0.0, None)


def box_layer_bounds(depth, bottoms, deepest_to_bottom=True):
    """Top and bottom depth of the integration interval of every layer of every box.

    Returns (top, bot), both (n_box, n_layer). Layers a box does not have are
    empty (top = bot = inf). The first layer starts at -inf so that the elevation of
    the free surface is counted. With `deepest_to_bottom` the deepest layer of each
    box extends to +inf; otherwise it ends at its nominal (truncated) bottom.
    """
    bottoms = np.asarray(bottoms, dtype=float)
    tops = layer_tops(bottoms)
    depth = np.asarray(depth, dtype=float)
    n_box, n_layer = len(depth), len(bottoms)

    top = np.broadcast_to(tops, (n_box, n_layer)).copy()
    bot = np.minimum(bottoms[None, :], depth[:, None])
    top[:, 0] = -np.inf

    count = n_layers(depth, bottoms)
    k = np.arange(n_layer)[None, :]
    last = k == (count[:, None] - 1)
    if deepest_to_bottom:
        bot = np.where(last, np.inf, bot)
    empty = k >= count[:, None]
    top = np.where(empty, np.inf, top)
    bot = np.where(empty, np.inf, bot)
    return top, bot


def overlap(z_w, top, bot):
    """Length of the overlap between each model level and each layer.

    `z_w` are the level interfaces in metres, negative below the surface, shaped
    (n_s + 1,) + columns with index 0 at the bottom (the CROCO order). `top` and `bot`
    are the layer bounds in depth (positive down), either (n_layer,) for all columns
    or (n_layer,) + columns for per-column bounds. Returns (n_layer, n_s) + columns.
    """
    depth = -np.asarray(z_w, dtype=float)
    lev_bot = depth[:-1]  # deeper face of each level
    lev_top = depth[1:]  # shallower face of each level
    top, bot = np.asarray(top, dtype=float), np.asarray(bot, dtype=float)
    if top.ndim == 1:
        shape = (-1,) + (1,) * depth.ndim
        top, bot = top.reshape(shape), bot.reshape(shape)
    else:
        top, bot = top[:, None], bot[:, None]
    lo = np.maximum(lev_top[None], top)
    hi = np.minimum(lev_bot[None], bot)
    return np.clip(hi - lo, 0.0, None)


def remap_scalar(values, ov, valid=None):
    """Thickness-weighted mean of `values` (n_s,) + columns in every layer.

    `ov` is the output of `overlap`. `valid` is an optional boolean array over the
    columns (for example the water mask). Returns (n_layer,) + columns, NaN where a
    layer has no water.
    """
    w = ov if valid is None else ov * np.asarray(valid)[None, None]
    total = w.sum(axis=1)
    num = (w * np.asarray(values)[None]).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(total > 0, num / total, np.nan)
