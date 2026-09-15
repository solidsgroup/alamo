# Propellant regression animation

[View the GIF](propellant_regression_q500_ap40.gif).

This is the selected cp=2130 J/(kg K) simulation at 40% AP by volume and
q=500 cal/(cm² s), case `q500_t0.400_e0.25_n8_p1e+08_ap_frozen`.
Its 41 saved snapshots span 0–4.763 ms and show 58.6 µm of recession.
The GIF loops over 7.7 seconds, including pauses at the beginning and end.

Color shows the saved propellant temperature on a fixed 300–900 K scale;
opacity follows the solid fraction. The moving line is the interpolated
solid-fraction=0.5 surface. The dashed line marks its initial position.
The transverse direction is stretched for visibility, with physical
coordinate labels retained. Gas chemistry is frozen. No flames, particle
motion or intermediate simulation snapshots are synthesized.

The neighboring CSV contains each frame's time, recession, surface
temperature and plotfile path; the JSON records input/parameter hashes,
frame timing and rendering details. PNG files provide a final still and
a beginning/middle/end preview.

Reproduce from the repository root with the existing Python environment:

```bash
/home/esandall/Software/anaconda3/bin/python validation_homogeneous_chen2002/analysis/render_regression_gif.py \
  --case validation_homogeneous_chen2002/runs/q500_t0.400_e0.25_n8_p1e+08_ap_frozen \
  --output validation_homogeneous_chen2002/animations/propellant_regression_q500_ap40.gif
```
