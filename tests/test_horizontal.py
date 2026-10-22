import numpy as np
import pytest
from shapely.geometry import box

from croco2atlantis import horizontal as hz

BOTTOMS = (20, 70, 200)


def grid4():
    """4 x 4 cells of 1 km, centres at 500, 1500, 2500, 3500 m."""
    c = np.arange(4) * 1000.0 + 500.0
    y, x = np.meshgrid(c, c, indexing="ij")
    return x, y, np.ones((4, 4), bool)


def test_assign_cells_and_counts():
    x, y, mask = grid4()
    mask[0, 0] = False  # land
    polys = [box(0, 0, 2000, 4000), box(2000, 0, 4000, 4000), box(10, 10, 20, 20)]
    cm = hz.CellBoxMap(hz.assign_cells(x, y, mask, polys), 3)
    assert cm.box[0, 0] == -1 and cm.box[0, 1] == 0 and cm.box[0, 3] == 1
    np.testing.assert_array_equal(cm.counts, [7, 8, 0])
    assert set(cm.columns[0]) == {i for i in range(16) if (i % 4) < 2 and i != 0}


def test_overlapping_polygons_first_wins():
    x, y, mask = grid4()
    polys = [box(0, 0, 2000, 4000), box(0, 0, 4000, 4000)]
    box_ = hz.assign_cells(x, y, mask, polys)
    assert (box_[:, :2] == 0).all() and (box_[:, 2:] == 1).all()


def _weights(depth):
    x, y, mask = grid4()
    cm = hz.CellBoxMap(hz.assign_cells(x, y, mask, [box(0, 0, 2000, 4000), box(2000, 0, 4000, 4000)]), 2)
    # 4 levels, interfaces at 0, -10, -30, -60, -100 m everywhere, but the right half is 50 m deep
    z = np.array([0, -10, -30, -60, -100], float)[::-1][:, None, None] * np.ones((5, 4, 4))
    z[:, :, 2:] = np.maximum(z[:, :, 2:], -50.0)
    area = np.full((4, 4), 1e6)
    return hz.box_layer_weights(z, cm, area, np.asarray(depth, float), BOTTOMS), z


def test_volume_weighted_mean_and_volume():
    w, z = _weights([100.0, 50.0])
    values = np.zeros((4, 4, 4))
    values[:] = np.array([1.0, 2.0, 3.0, 4.0])[:, None, None]  # level 0 is the bottom
    out = w.apply(values)
    # box 0, layer 0 (to 20 m): levels 0-10 (value 4) and 10-30 split at 20 -> (4*10 + 3*10)/20
    assert out[0, 0] == pytest.approx(3.5)
    # volume of box 0 layer 0 = 8 cells * 1 km2 * 20 m
    assert w.volume[0, 0] == pytest.approx(8 * 1e6 * 20)
    # box 1 is 50 m deep: layers 0 and 1 only (to 50 m)
    assert np.isfinite(out[1, :2]).all() and np.isnan(out[1, 2])


def test_deepest_layer_takes_the_whole_column():
    w, _ = _weights([100.0, 50.0])
    # box 0 mean depth 100 -> 3 layers; the last one (70 m down) takes everything below
    assert w.volume[0, 2] == pytest.approx(8 * 1e6 * 30)  # columns end at 100 m
    w2, _ = _weights([30.0, 50.0])
    assert w2.volume[0, 1] == pytest.approx(8 * 1e6 * 80)  # layer 2 (20 m down) to the bottom at 100 m


def test_fill_missing_uses_nearest_box_with_data():
    means = np.array([[1.0, 2.0, np.nan], [np.nan, np.nan, np.nan], [5.0, np.nan, np.nan], [np.nan, np.nan, np.nan]])
    depth = np.array([100.0, 100.0, 30.0, 0.0])  # box 3 is an island: no layers, never filled
    centroids = np.array([[0.0, 0], [10, 0], [100, 0], [1000, 0]])
    filled, was = hz.fill_missing(means, depth, BOTTOMS, centroids)
    # box 1 is nearest to box 0 for every layer that box 0 has data in
    np.testing.assert_allclose(filled[1, :2], [1.0, 2.0])
    # layer 2 has no data anywhere, so nothing to fill from
    assert np.isnan(filled[1, 2]) and np.isnan(filled[0, 2])
    # box 2 (30 m deep) has 2 layers; its layer 1 is filled from box 0, the only box with data there
    assert filled[2, 1] == 2.0
    np.testing.assert_array_equal(was[1], [True, True, False])
    np.testing.assert_array_equal(was[2], [False, True, False])
    assert not was[0].any() and not was[3].any()
    assert np.isnan(filled[3]).all()


def test_fill_missing_max_distance_leaves_far_boxes_empty():
    means = np.array([[1.0], [np.nan], [np.nan]])
    centroids = np.array([[0.0, 0], [10, 0], [5000, 0]])
    filled, was = hz.fill_missing(means, np.full(3, 10.0), BOTTOMS, centroids, max_distance=100.0)
    assert filled[1, 0] == 1.0 and was[1, 0]
    assert np.isnan(filled[2, 0]) and not was[2, 0]


def test_fill_missing_takes_the_whole_profile_from_one_box():
    # box 2 needs layers 0 and 1; box 0 is nearer but only has layer 0, box 1 has both
    means = np.array([[1.0, np.nan, np.nan], [3.0, 4.0, np.nan], [np.nan, np.nan, np.nan]])
    depth = np.array([10.0, 100.0, 50.0])
    centroids = np.array([[0.0, 0], [100, 0], [10, 0]])
    filled, was = hz.fill_missing(means, depth, BOTTOMS, centroids)
    np.testing.assert_allclose(filled[2, :2], [3.0, 4.0])
    np.testing.assert_array_equal(was[2], [True, True, False])


def test_centre_on_a_shared_edge_goes_to_the_first_polygon():
    x, y = np.array([[2000.0]]), np.array([[500.0]])
    box_ = hz.assign_cells(x, y, np.ones((1, 1), bool), [box(0, 0, 2000, 4000), box(2000, 0, 4000, 4000)])
    assert box_[0, 0] == 0


def test_cells_outside_or_on_islands_move_to_the_nearest_water_box():
    x = np.array([[500.0, 1500.0, 2500.0, 4200.0, 9000.0]])
    y = np.full_like(x, 500.0)
    mask = np.ones(x.shape, bool)
    mask[0, 0] = False  # land stays -1
    polys = [box(0, 0, 2000, 1000), box(2000, 0, 3000, 1000), box(3000, 0, 4000, 1000)]
    water = np.array([True, False, True])  # box 1 is an island
    b = hz.assign_cells(x, y, mask, polys)
    np.testing.assert_array_equal(b, [[-1, 0, 1, -1, -1]])
    moved_box, moved = hz.move_to_water_boxes(b, x, y, mask, polys, water, max_distance=1000.0)
    # island cell at 2500 is 500 m from box 0 and box 2: either is right; 4200 -> box 2; 9000 is too far
    assert moved_box[0, 2] in (0, 2)
    np.testing.assert_array_equal(moved_box[0, [0, 1, 3, 4]], [-1, 0, 2, -1])
    np.testing.assert_array_equal(moved, [[False, False, True, True, False]])


def test_apply_rejects_missing_values_in_water():
    w, _ = _weights([100.0, 50.0])
    values = np.ones((4, 4, 4))
    values[0, 1, 1] = np.nan
    with pytest.raises(ValueError, match="box 0"):
        w.apply(values)
