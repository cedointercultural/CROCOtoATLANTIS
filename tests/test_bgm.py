import numpy as np

from croco2atlantis.geometry.bgm import read_bgm

# Two 1 km squares side by side sharing one face at x = 1000.
BGM = """# test
projection +proj=laea +lat_0=26 +lon_0=-110 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs
nbox 2
nface 1
maxwcbotz -50.0
box0.label\tBox0
box0.inside\t500 500
box0.nconn\t1
box0.iface\t0
box0.ibox\t1
box0.botz\t-20.0
box0.area\t1000000
box0.vertmix\tNaN
box0.horizmix\tNaN
box0.vert 0 0
box0.vert 1000 0
box0.vert 1000 1000
box0.vert 0 1000
box0.vert 0 0
box1.label\tBox1
box1.inside\t1500 500
box1.nconn\t1
box1.iface\t0
box1.ibox\t0
box1.botz\t0.0
box1.area\t1000000
box1.vertmix\tNaN
box1.horizmix\tNaN
box1.vert 1000 0
box1.vert 2000 0
box1.vert 2000 1000
box1.vert 1000 1000
box1.vert 1000 0
face0.p1 1000 0
face0.p2 1000 1000
face0.length 1000
face0.cs 0 1
face0.lr 0 1
"""


def test_read_bgm(tmp_path):
    path = tmp_path / "test.bgm"
    path.write_text(BGM)
    bgm = read_bgm(path)
    assert (bgm.nbox, bgm.nface) == (2, 1)
    assert bgm.maxwcbotz == -50.0
    np.testing.assert_allclose(bgm.depth, [20.0, 0.0])
    np.testing.assert_allclose(bgm.area, [1e6, 1e6])
    assert np.isnan(bgm.vertmix).all() and np.isnan(bgm.horizmix).all()
    assert bgm.boxes[0].vert.shape == (5, 2)
    assert [p.area for p in bgm.box_polygons()] == [1e6, 1e6]


def test_face_normal_points_to_left_box(tmp_path):
    path = tmp_path / "test.bgm"
    path.write_text(BGM)
    bgm = read_bgm(path)
    # p1 -> p2 goes north, so the left side is west, where box 0 (lr[0]) is.
    np.testing.assert_allclose(bgm.face_normals(), [[-1.0, 0.0]])
    np.testing.assert_array_equal(bgm.face_lr(), [[0, 1]])
    np.testing.assert_allclose(bgm.face_midpoints(), [[1000.0, 500.0]])


def test_mixing_coefficients(tmp_path):
    path = tmp_path / "test.bgm"
    text = BGM.replace("box0.vertmix\tNaN", "box0.vertmix\t1.0E-6").replace("box0.horizmix\tNaN", "box0.horizmix\t1.0")
    path.write_text(text)
    bgm = read_bgm(path)
    np.testing.assert_allclose(bgm.vertmix[:1], [1e-6])
    np.testing.assert_allclose(bgm.horizmix[:1], [1.0])
    assert np.isnan(bgm.vertmix[1])


def test_project_unproject_roundtrip(tmp_path):
    path = tmp_path / "test.bgm"
    path.write_text(BGM)
    bgm = read_bgm(path)
    x, y = bgm.project(-110.5, 26.5)
    lon, lat = bgm.unproject(x, y)
    np.testing.assert_allclose([lon, lat], [-110.5, 26.5], atol=1e-8)
