# Segmentation methods panel (Panel A) — authoring spec

The anatomical L3–L5 segmentation panel is the one methods figure that cannot be
generated from tabular data. It is authored from real imaging, not synthesized, so
it is documented here rather than faked in `methods_figures.py`.

## Content
A representative axial (and optional sagittal) lumbar MRI at the L4 level with the
three analyzed muscle groups overlaid as color masks, matching the JAMA palette:

| Muscle group | Color (pal_jama) | Role |
|---|---|---|
| Iliopsoas (bilateral) | blue `#00A1D5` | reference muscle; strongest signal predictor |
| Deep paraspinal / erector+multifidus (bilateral) | green `#79AF97` | opposite-direction volume finding |
| Gluteus medius (bilateral) | orange `#DF8F44` | multi-muscle block |

## Provenance to state in the caption
- Segmentation: **3D Slicer** with **TotalSegmentator v5.8.1**, at **L3–L5**.
- Metrics extracted per region: labeled-mask volume (LM) and sub-volume after
  fat-threshold exclusion (sv); mean/median T2 signal intensity; the sv/LM ratio is
  a lean-fraction (muscle-quality) proxy. Volumes normalized to L4 vertebral body
  volume.
- Shared with the companion "fatproxy" paper (same segmentation methodology; Paper A
  uses iliopsoas + deep paraspinal, Paper B adds gluteus medius).

## How to produce
Export representative masked slices from 3D Slicer (de-identified, no burned-in
PHI), assemble at JAMA single/1.5-column width, 300 dpi. Place the file at
`figures/MethodsA_segmentation.(png|pdf|svg)`. A single de-identified illustrative
case is sufficient; do not include identifiers or accession numbers.
