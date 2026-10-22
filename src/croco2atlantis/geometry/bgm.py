"""Read an ATLANTIS BGM file.

A BGM describes the boxes (polygons) and the faces shared between them. All
coordinates are in the projected system declared in the file (metres).

Face convention: `lr` holds the boxes on the left and right of the face looking from
`p1` to `p2`. `Bgm.face_normals` returns the unit normal pointing to the left, so a
positive flux along it goes from the right box to the left box.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class Box:
    id: int
    label: str
    inside: np.ndarray  # (2,) interior point
    botz: float  # mean depth, negative below the surface
    area: float  # m2
    vertmix: float  # vertical mixing coefficient used by ATLANTIS, NaN if not set
    horizmix: float  # horizontal mixing coefficient used by ATLANTIS, NaN if not set
    iface: np.ndarray  # ids of the faces of the box
    ibox: np.ndarray  # box on the other side of each face
    vert: np.ndarray  # (n, 2) closed ring of vertices


@dataclass
class Face:
    id: int
    p1: np.ndarray
    p2: np.ndarray
    length: float
    cs: np.ndarray  # components of the face angle, as stored in the file
    lr: np.ndarray  # (left box, right box)


@dataclass
class Bgm:
    projection: str
    boxes: list[Box]
    faces: list[Face]
    maxwcbotz: float | None = None
    bnd_vert: np.ndarray = field(default_factory=lambda: np.empty((0, 2)))

    @property
    def nbox(self):
        return len(self.boxes)

    @property
    def nface(self):
        return len(self.faces)

    @property
    def depth(self):
        """Mean depth of each box in metres, positive downwards (0 for islands)."""
        return -np.array([b.botz for b in self.boxes])

    @property
    def area(self):
        return np.array([b.area for b in self.boxes])

    @property
    def vertmix(self):
        return np.array([b.vertmix for b in self.boxes])

    @property
    def horizmix(self):
        return np.array([b.horizmix for b in self.boxes])

    def box_polygons(self):
        """Shapely polygon of every box, in the projected system of the file."""
        from shapely.geometry import Polygon

        return [Polygon(b.vert) for b in self.boxes]

    def face_midpoints(self):
        return np.array([(f.p1 + f.p2) / 2 for f in self.faces])

    def face_lengths(self):
        return np.array([f.length for f in self.faces])

    def face_lr(self):
        return np.array([f.lr for f in self.faces], dtype=int)

    def face_normals(self):
        """Unit normal of every face pointing to the left of the p1 -> p2 direction."""
        d = np.array([f.p2 - f.p1 for f in self.faces])
        n = np.stack([-d[:, 1], d[:, 0]], axis=1)
        return n / np.linalg.norm(n, axis=1, keepdims=True)

    def project(self, lon, lat):
        """Project longitude and latitude (degrees) to the coordinates of the BGM."""
        from pyproj import Transformer

        tr = Transformer.from_crs("EPSG:4326", self.projection, always_xy=True)
        return tr.transform(lon, lat)

    def unproject(self, x, y):
        """Inverse of `project`: coordinates of the BGM to longitude and latitude (degrees)."""
        from pyproj import Transformer

        tr = Transformer.from_crs(self.projection, "EPSG:4326", always_xy=True)
        return tr.transform(x, y)


def _floats(tokens):
    return np.array([float(t) for t in tokens])


def read_bgm(path):
    """Read a BGM file and return a `Bgm`."""
    header = {}
    box_attrs = {}
    face_attrs = {}
    box_verts = {}
    bnd_vert = []

    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)  # keys and values are separated by spaces or tabs
        key = parts[0]
        rest = parts[1].strip() if len(parts) > 1 else ""
        if key == "bnd_vert":
            bnd_vert.append(_floats(rest.split()))
        elif key.startswith("box") and "." in key:
            num, attr = key[3:].split(".", 1)
            i = int(num)
            if attr == "vert":
                box_verts.setdefault(i, []).append(_floats(rest.split()))
            else:
                box_attrs.setdefault(i, {})[attr] = rest
        elif key.startswith("face") and "." in key:
            num, attr = key[4:].split(".", 1)
            face_attrs.setdefault(int(num), {})[attr] = rest
        else:
            header[key] = rest

    nbox = int(header["nbox"])
    nface = int(header["nface"])
    if sorted(box_attrs) != list(range(nbox)):
        raise ValueError(f"expected boxes 0..{nbox - 1}, found {len(box_attrs)} boxes")
    if sorted(face_attrs) != list(range(nface)):
        raise ValueError(f"expected faces 0..{nface - 1}, found {len(face_attrs)} faces")

    boxes = []
    for i in range(nbox):
        a = box_attrs[i]
        boxes.append(
            Box(
                id=i,
                label=a.get("label", f"Box{i}"),
                inside=_floats(a["inside"].split()),
                botz=float(a["botz"]),
                area=float(a["area"]),
                vertmix=float(a.get("vertmix", "nan")),
                horizmix=float(a.get("horizmix", "nan")),
                iface=np.array(a["iface"].split(), dtype=int),
                ibox=np.array(a["ibox"].split(), dtype=int),
                vert=np.array(box_verts[i]),
            )
        )

    faces = []
    for i in range(nface):
        a = face_attrs[i]
        faces.append(
            Face(
                id=i,
                p1=_floats(a["p1"].split()),
                p2=_floats(a["p2"].split()),
                length=float(a["length"]),
                cs=_floats(a["cs"].split()),
                lr=np.array(a["lr"].split(), dtype=int),
            )
        )

    return Bgm(
        projection=header["projection"],
        boxes=boxes,
        faces=faces,
        maxwcbotz=float(header["maxwcbotz"]) if "maxwcbotz" in header else None,
        bnd_vert=np.array(bnd_vert) if bnd_vert else np.empty((0, 2)),
    )
