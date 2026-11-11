# Figure specifications

Main figures use a 7.2-inch double-column layout with 8-point base typography, a restrained colour palette, and consistent panel lettering. Colours distinguish incompatibility values in endpoint plots; line styles and markers distinguish mutation kernels. Exposure curves use colours for redundancy and facets for kernel and incompatibility.

Replicate points, means, and sample-standard-deviation error bars are shown in the endpoint grids. Exposure curves show replicate means and sample-standard-deviation bands after stepwise alignment. Corresponding panels in the two endpoint grids share y-axis limits.

Every assembled presentation has PDF, SVG, and 600-dpi PNG exports. SVG text remains editable, and plot elements are vector objects. Lower-resolution previews support browsing in the README. The exposure grid is designed for double-column placement; captions are provided separately in `docs/figure_captions.md`.

`scripts/check_release_artifacts.py` checks file availability, image dimensions, resolution metadata, vector exports, relative links, and checksums. `artifact_manifest.json` associates exported artwork with endpoint tables and trajectory files.
