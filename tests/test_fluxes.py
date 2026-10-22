from types import SimpleNamespace

import numpy as np

from croco2atlantis import fluxes
from croco2atlantis.horizontal import CellBoxMap
from croco2atlantis.readers.croco import CrocoGrid

BOTTOMS = [18, 54, 270]


def make_grid(ny, nx, d=1000.0, open_boundary=None):
    """Grid of `ny` x `nx` square cells of side `d` m, all water, closed ghost ring."""
    mask = np.zeros((ny, nx), bool)
    mask[1:-1, 1:-1] = True
    return CrocoGrid(
        lon=np.zeros((ny, nx)),
        lat=np.zeros((ny, nx)),
        h=np.full((ny, nx), 100.0),
        mask=mask,
        angle=np.zeros((ny, nx)),
        pm=np.full((ny, nx), 1 / d),
        pn=np.full((ny, nx), 1 / d),
        open_boundary=open_boundary,
    )


def flat_levels(h, n_s=4, zeta=0.0):
    """Interfaces (n_s + 1, ...) of evenly spaced levels from -h to zeta."""
    frac = np.linspace(-1, 0, n_s + 1).reshape((-1,) + (1,) * np.ndim(h))
    return zeta + (h + zeta) * frac


def halves(ny, nx):
    """Interior cells: box 0 on the west half, box 1 on the east half."""
    box = np.full((ny, nx), -1)
    box[1:-1, 1 : nx // 2] = 0
    box[1:-1, nx // 2 : -1] = 1
    return CellBoxMap(box, 2)


def test_uniform_current_gives_transport_per_layer():
    ny, nx = 5, 6  # 3 interior rows of 1 km
    grid = make_grid(ny, nx)
    faces = fluxes.cell_faces(halves(ny, nx), grid)
    np.testing.assert_array_equal(faces.pairs, [[0, 1]])
    u = np.ones((4, ny, nx - 1))  # 1 m/s towards +xi, from box 0 to box 1
    v = np.zeros((4, ny - 1, nx))
    q = fluxes.face_transport(faces, u, v, flat_levels(grid.h), BOTTOMS)
    # 3 faces of 1 km; layers of 18, 36 and the remaining 46 m of the 100 m column
    np.testing.assert_allclose(q[0], 3000.0 * np.array([18.0, 36.0, 46.0]))
    inflow = fluxes.net_inflow(faces, q)
    np.testing.assert_allclose(inflow[1], q[0])
    np.testing.assert_allclose(inflow[0], -q[0])


def test_sign_follows_the_pair_order():
    ny, nx = 5, 6
    grid = make_grid(ny, nx)
    cm = halves(ny, nx)
    swapped = CellBoxMap(np.where(cm.box >= 0, 1 - cm.box, -1), 2)  # box 1 west, box 0 east
    faces = fluxes.cell_faces(swapped, grid)
    u = np.ones((4, ny, nx - 1))
    q = fluxes.face_transport(faces, u, np.zeros((4, ny - 1, nx)), flat_levels(grid.h), BOTTOMS)
    assert (q < 0).all()  # the flow goes from box 1 to box 0


def test_face_depth_is_the_mean_of_its_two_cells_and_counts_the_free_surface():
    ny, nx = 3, 4  # one interior row
    grid = make_grid(ny, nx)
    h = np.full((ny, nx), 100.0)
    h[:, nx // 2 :] = 50.0
    box = np.full((ny, nx), -1)
    box[1, 1], box[1, 2] = 0, 1
    faces = fluxes.cell_faces(CellBoxMap(box, 2), grid)
    u = np.ones((4, ny, nx - 1))
    q = fluxes.face_transport(faces, u, np.zeros((4, ny - 1, nx)), flat_levels(h, zeta=2.0), BOTTOMS)
    # the face is 75 m deep plus 2 m of free surface: 20, 36 and 21 m
    np.testing.assert_allclose(q[0], 1000.0 * np.array([20.0, 36.0, 21.0]))


def test_nondivergent_flow_closes_every_box():
    rng = np.random.default_rng(0)
    ny, nx, d, h = 9, 11, 1000.0, 100.0
    grid = make_grid(ny, nx, d)
    # stream function on the cell corners, zero on the edge of the interior water
    psi = np.zeros((ny + 1, nx + 1))
    psi[2:-2, 2:-2] = rng.normal(size=(ny - 3, nx - 3)) * 1e5
    U = psi[1:-1, 1:-1][1:] - psi[1:-1, 1:-1][:-1]  # transport through the u faces of rows 1..ny-2
    V = -(psi[1:-1, 1:-1][:, 1:] - psi[1:-1, 1:-1][:, :-1])
    u = np.zeros((4, ny, nx - 1))
    v = np.zeros((4, ny - 1, nx))
    u[:, 1:-1, :] = (U / (d * h))[None]
    v[:, :, 1:-1] = (V / (d * h))[None]
    box = np.where(grid.mask, rng.integers(0, 4, size=(ny, nx)), -1)
    faces = fluxes.cell_faces(CellBoxMap(box, 4), grid)
    q = fluxes.face_transport(faces, u, v, flat_levels(np.full((ny, nx), h)), BOTTOMS)
    assert np.abs(q).max() > 1e3
    np.testing.assert_allclose(fluxes.net_inflow(faces, q), 0.0, atol=1e-6)


def test_open_boundary_exchanges_with_box_minus_one():
    ny, nx = 5, 6
    ob = np.zeros((ny, nx), bool)
    ob[1:-1, 0] = True  # western open boundary
    grid = make_grid(ny, nx, open_boundary=ob)
    box = np.where(grid.mask, 0, -1)
    faces = fluxes.cell_faces(CellBoxMap(box, 1), grid)
    np.testing.assert_array_equal(faces.pairs, [[-1, 0]])
    u = np.ones((4, ny, nx - 1))  # eastwards: into the domain through the west
    q = fluxes.face_transport(faces, u, np.zeros((4, ny - 1, nx)), flat_levels(grid.h), BOTTOMS)
    np.testing.assert_allclose(q.sum(), 3000.0 * 100.0)
    np.testing.assert_allclose(fluxes.net_inflow(faces, q)[0], q[0])


def test_land_faces_carry_nothing():
    ny, nx = 5, 6
    grid = make_grid(ny, nx)
    grid.mask[2, 3] = False  # land on box 1 side of the middle row
    faces = fluxes.cell_faces(halves(ny, nx), grid)
    assert len(faces.length) == 2


def test_pairs_without_face():
    box = np.array([[0, 1, 2]])
    grid = make_grid(3, 5)
    full = np.full((3, 5), -1)
    full[1, 1:4] = box
    faces = fluxes.cell_faces(CellBoxMap(full, 3), grid)
    bgm = SimpleNamespace(face_lr=lambda: np.array([[1, 0]]))
    np.testing.assert_array_equal(fluxes.pairs_without_face(faces, bgm), [[1, 2]])
