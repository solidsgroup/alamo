# Native-face isolation

Classification: **native feature present, but native-smooth/nodal-smooth**.
In the plan's categorical decision table this follows the “native feature
remains” branch, not the “native-smooth/nodal-rough” branch.

On the synchronized state, both symmetry-axis native normal tractions are
strictly monotone through `0.1 <= phi <= 0.9`; monotone-envelope overshoot is
zero. `eta` is exactly `1.0` throughout both sampled bands and its neighboring
face delta is zero.

| Axis | Native stress span | Native Nyquist | Nyquist/span |
|---|---:|---:|---:|
| bottom | `1.4378 MPa` | `49.19 Pa` | `3.42e-5` |
| right | `1.7847 MPa` | `90.24 Pa` | `5.06e-5` |

The saved polar stress and exact face-average nodal reconstruction agree to
`2.35e-6 Pa` on the bottom axis, `1.49e-7 Pa` on the right axis, and
`1.88e-5 Pa` over the full finest-level `phi` band. Therefore neither the
VisIt polar projection nor the saved nodal tensor introduces the visible
line.

A single native face and its low/high face average can differ by as much as
`60.2 kPa` (bottom) or `144.9 kPa` (right), as expected when comparing a
face-centered value with a node-centered average across a steep smooth
gradient. That centering difference is smooth and does not create alternating
overshoot.

Conclusion: solution 2 would change the centering and name of the displayed
quantity, but would not remove the visible `phi`-following transition. Adding
a production “native radial” scalar would also be misleading off-axis because
the coordinate faces are not radial there. No production diagnostic is
retained.

Independent artifacts:

- [`native-face/profiles.csv`](native-face/profiles.csv)
- [`native-face/profiles.png`](native-face/profiles.png)
- [`native-face/metrics.json`](native-face/metrics.json)

