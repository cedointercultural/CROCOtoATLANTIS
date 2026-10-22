import numpy as np
import pytest

from croco2atlantis import vertical

BOTTOMS = (18, 54, 270, 684, 918, 1782, 3100)


def test_layer_tops():
    np.testing.assert_allclose(vertical.layer_tops(BOTTOMS), [0, 18, 54, 270, 684, 918, 1782])


def test_n_layers_and_thickness_follow_the_box_rule():
    depth = np.array([0.0, 10.0, 712.45, 648.09, 3101.9])
    np.testing.assert_array_equal(vertical.n_layers(depth, BOTTOMS), [0, 1, 5, 4, 7])
    th = vertical.nominal_thickness(depth, BOTTOMS)
    np.testing.assert_allclose(th[0], 0)
    np.testing.assert_allclose(th[1], [10, 0, 0, 0, 0, 0, 0])
    # box 67 of the Excel: 18, 36, 216, 414 and a last layer of 28.45 m
    np.testing.assert_allclose(th[2], [18, 36, 216, 414, 28.45, 0, 0], atol=1e-9)
    # box 86: 18, 36, 216 and a last layer of 378.09 m
    np.testing.assert_allclose(th[3], [18, 36, 216, 378.09, 0, 0, 0], atol=1e-9)
    # a box deeper than the last bottom still has 7 layers
    np.testing.assert_allclose(th[4], [18, 36, 216, 414, 234, 864, 1318])


def test_bounds_deepest_layer_takes_the_whole_column():
    depth = np.array([712.45])
    top, bot = vertical.box_layer_bounds(depth, BOTTOMS, deepest_to_bottom=True)
    assert top[0, 0] == -np.inf
    np.testing.assert_allclose(top[0, 1:5], [18, 54, 270, 684])
    np.testing.assert_allclose(bot[0, :4], [18, 54, 270, 684])
    assert bot[0, 4] == np.inf
    assert np.isinf(top[0, 5:]).all() and np.isinf(bot[0, 5:]).all()
    _, bot_nominal = vertical.box_layer_bounds(depth, BOTTOMS, deepest_to_bottom=False)
    assert bot_nominal[0, 4] == pytest.approx(712.45)


def _columns():
    # two columns, 4 levels each, interfaces at 0, -10, -30, -60, -100 and 0, -5, -15, -30, -50
    return np.array([[0, 0], [-10, -5], [-30, -15], [-60, -30], [-100, -50]], float)[::-1]


def test_overlap_adds_up_to_level_thickness():
    z_w = _columns()
    bottoms = (20, 70, 200)
    top, bot = vertical.box_layer_bounds(np.array([500.0]), bottoms)
    ov = vertical.overlap(z_w, top[0], bot[0])
    assert ov.shape == (3, 4, 2)
    thickness = np.diff(z_w, axis=0)
    np.testing.assert_allclose(ov.sum(axis=0), thickness)


def test_overlap_per_layer_values():
    z_w = _columns()
    ov = vertical.overlap(z_w, np.array([-np.inf, 20.0]), np.array([20.0, np.inf]))
    # column 0, levels from the bottom: 60-100, 30-60, 10-30, 0-10; the 20 m cut splits 10-30
    np.testing.assert_allclose(ov[0, :, 0], [0, 0, 10, 10])
    np.testing.assert_allclose(ov[1, :, 0], [40, 30, 10, 0])


def test_remap_scalar_is_thickness_weighted_and_conservative():
    z_w = _columns()
    values = np.array([[4.0, 4.0], [3.0, 3.0], [2.0, 2.0], [1.0, 1.0]])[::-1]  # surface warmest
    ov = vertical.overlap(z_w, np.array([-np.inf, 20.0]), np.array([20.0, np.inf]))
    out = vertical.remap_scalar(values, ov)
    # column 0, layer 0: (4*10 + 3*10) / 20; layer 1: (3*10 + 2*30 + 1*40) / 80
    assert out[0, 0] == pytest.approx(3.5)
    assert out[1, 0] == pytest.approx((30 + 60 + 40) / 80)


def test_remap_constant_field_stays_constant_and_empty_layers_are_nan():
    z_w = _columns()
    values = np.full((4, 2), 7.5)
    top, bot = vertical.box_layer_bounds(np.array([500.0]), (20, 70, 200, 400))
    ov = vertical.overlap(z_w, top[0], bot[0])
    out = vertical.remap_scalar(values, ov)
    np.testing.assert_allclose(out[:2], 7.5)
    assert out[2, 0] == 7.5 and np.isnan(out[2, 1])  # the second column ends at 50 m
    assert np.isnan(out[3]).all()  # nothing reaches 200-400 m


def test_remap_respects_the_valid_mask():
    z_w = _columns()
    values = np.array([[1.0, 5.0]] * 4)
    ov = vertical.overlap(z_w, np.array([-np.inf]), np.array([np.inf]))
    out = vertical.remap_scalar(values, ov, valid=np.array([True, False]))
    assert out[0, 0] == pytest.approx(1.0)
    assert np.isnan(out[0, 1])
