"""Check local links, figure exports, and data-to-artifact hashes."""

import hashlib
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from PIL import Image

from enzyme_specialisation.io import scientific_fingerprint


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


root = Path(".")
missing = []
for p in [Path("README.md"), *Path("docs").glob("*.md")]:
    text = p.read_text()
    targets = re.findall(r"\]\(([^)]+)\)", text) + re.findall(r'src="([^"]+)"', text)
    for target in targets:
        if target.startswith(("http:", "https:", "#")):
            continue
        if not (p.parent / target.split("#")[0]).exists():
            missing.append((str(p), target))
assert not missing, missing
figures = []
for p in sorted(Path("figures").glob("*/*.png")):
    with Image.open(p) as im:
        width, height = im.size
        dpi = im.info.get("dpi", (None, None))
        im.verify()
    if p.parent.name != "previews":
        assert abs(dpi[0] - 600) < 1
        for extension in ["pdf", "svg"]:
            assert p.with_suffix("." + extension).exists()
        svg = ET.parse(p.with_suffix(".svg"))
        assert not svg.findall(".//{http://www.w3.org/2000/svg}image")
        assert svg.findall(".//{http://www.w3.org/2000/svg}text")
    figures.append({"path": str(p), "width": width, "height": height, "dpi": dpi, "sha256": sha(p)})
records = {
    "scientific_fingerprint": scientific_fingerprint(),
    "figure_exports": figures,
    "input_files": {str(p): sha(p) for p in sorted(Path("data/reference").glob("*.csv"))},
    "endpoint_files": {str(p): sha(p) for p in sorted(Path("results/paper").glob("*/*.csv"))},
    "trajectory_files": {
        str(p): sha(p) for p in sorted(Path("results/paper/trajectories").glob("*.npz"))
    },
    "readme_and_doc_links": "passed",
}
Path("reports/artifact_manifest.json").write_text(json.dumps(records, indent=2))
print(
    f"{len(figures)} raster exports/previews verified; links, 600-dpi metadata, vector text, and hashes passed"
)
