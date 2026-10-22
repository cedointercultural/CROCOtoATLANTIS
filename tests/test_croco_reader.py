import numpy as np
import pytest
import xarray as xr

from croco2atlantis.readers import croco


def params(vtransform, hc=15.0):
    n = 8
    s_w = np.linspace(-1, 0, n + 1)
    s_rho = (s_w[:-1] + s_w[1:]) / 2
    # any monotonic stretching curve with Cs(-1) = -1 and Cs(0) = 0; Vtransform 1 needs hc <= min(h)
    return croco.VerticalParams(vtransform, hc, s_rho, s_w, s_rho * np.abs(s_rho), s_w * np.abs(s_w))


@pytest.mark.parametrize("vtransform", [1, 2])
def test_s_depths_reach_the_surface_and_the_bottom(vtransform):
    p = params(vtransform)
    h = np.array([[20.0, 500.0], [3000.0, 18.0]])
    zeta = np.array([[0.3, -0.2], [0.0, 0.1]])
    z = croco.s_depths(h, zeta, p, "w")
    assert z.shape == (9, 2, 2)
    np.testing.assert_allclose(z[0], -h, atol=1e-9)
    np.testing.assert_allclose(z[-1], zeta, atol=1e-9)
    assert (np.diff(z, axis=0) > 0).all()
    zr = croco.s_depths(h, zeta, p, "rho")
    assert zr.shape == (8, 2, 2)
    assert ((zr > z[:-1]) & (zr < z[1:])).all()
    with pytest.raises(ValueError):
        croco.s_depths(h, zeta, p, "u")


def test_read_grid_masks_land_and_ghost_ring(tmp_path):
    ny, nx = 6, 5
    mask = np.ones((ny, nx))
    mask[3, 3] = 0
    ds = xr.Dataset(
        {
            "lon_rho": (("eta_rho", "xi_rho"), np.zeros((ny, nx))),
            "lat_rho": (("eta_rho", "xi_rho"), np.zeros((ny, nx))),
            "h": (("eta_rho", "xi_rho"), np.full((ny, nx), 100.0)),
            "mask_rho": (("eta_rho", "xi_rho"), mask),
            "angle": (("eta_rho", "xi_rho"), np.full((ny, nx), 0.5)),
            "pm": (("eta_rho", "xi_rho"), np.full((ny, nx), 1 / 2000.0)),
            "pn": (("eta_rho", "xi_rho"), np.full((ny, nx), 1 / 1000.0)),
        }
    )
    grid = croco.read_grid(ds)
    assert grid.mask.sum() == (ny - 2) * (nx - 2) - 1
    assert not grid.mask[0].any() and not grid.mask[:, -1].any() and not grid.mask[3, 3]
    # the water cells of the ghost ring are the open boundary
    assert grid.open_boundary.sum() == 2 * ny + 2 * nx - 4
    assert not (grid.open_boundary & grid.mask).any()
    np.testing.assert_allclose(grid.cell_area, 2e6)


def test_read_vertical_params_and_times():
    p = params(2, hc=50.0)
    ds = xr.Dataset(
        {
            "Vtransform": ((), np.float32(2)),
            "hc": ((), np.float32(50)),
            "s_rho": (("s_rho",), p.s_rho),
            "s_w": (("s_w",), p.s_w),
            "Cs_rho": (("s_rho",), p.cs_rho),
            "Cs_w": (("s_w",), p.cs_w),
            "time": (("time",), np.array([0.0, 43200.0])),
        }
    )
    got = croco.read_vertical_params(ds)
    assert got.vtransform == 2 and got.hc == 50.0
    np.testing.assert_allclose(got.cs_w, p.cs_w)
    t = croco.read_times(ds, "2000-01-01")
    assert str(t[1]) == "2000-01-01 12:00:00"
