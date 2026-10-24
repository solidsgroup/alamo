# Article-only reference extraction

Scientific source: the user-supplied *Homogenization issues and the combustion
of heterogeneous solid propellants*, Chen et al. (2002), pp. 2923–2929.
No additional publication is used for these artifacts.

`digitize_figure4.py` extracts the PDF's actual vector coordinates, using
`mutool draw -F trace` on PDF page 5 (printed page 2927). It reads the solid
Eq. 15 curves, dashed single-temperature curves, 12 Eq. 19 square markers,
24 two-dimensional DNS circles, and 24 three-dimensional DNS asterisks.
`figure4_provenance.json` records the source hash, coordinate calibration,
selected paths, and uncertainty convention. `figure4_extracted.pdf` is the
visually checked reconstruction. Endpoint DNS markers overlap and are not
independent heterogeneous validation cases.

The CSV decimals preserve vector geometry, not scientific accuracy. One full
printed stroke corresponds to approximately 0.005 cm/s or 0.0011 in volume
fraction. This is a conservative graphical resolution, not a confidence
interval; the source supplies no uncertainty estimates for the DNS points.

## What the article identifies

The captions of Figs. 5 and 6 specify E/R = 7500 K and A = 1036 cm/s for
binder; E/R = 11000 K and A = 94800 cm/s for AP. Eq. 12 gives
Ts = (E/R)/ln(A/r). Applying this to the three pure-component endpoints in
Fig. 4 and rearranging Eq. 14 gives

    q/r = alpha Ts - alpha beta,
    alpha = rho c,  beta = T0 + Q/c.

`infer_article_combinations.py` fits these two combinations only. It yields
approximately alpha_B = 1.155 MJ/(m³ K), beta_B = 143 K, alpha_AP =
2.435 MJ/(m³ K), beta_AP = −38 K. These are calibrated to Fig. 4 and cannot
provide independent endpoint validation. The article's approximate thermal
lengths (33 µm AP and 47 µm binder at r = 0.5 cm/s, p. 2927) can constrain
conductivity through k = alpha r delta, but remain approximate.

The separate rho, c, Q, and T0 are not uniquely identified. A chosen density
and T0 allow a numerical realization c = alpha/rho, Q = c(beta−T0), but the
choice must be disclosed. Graphical perturbation sensitivity is substantial,
especially for beta_AP; see `article_inferred_combinations.json`.

## Blend heat-capacity consistency

The article text on p. 2926 says that c, rho, and Q are mass-weighted.
With density ratio d = rho_AP/rho_B and volume fraction t, that rule gives

    alpha_blend = [(1−t)+d²t]/[(1−t)+dt]² * [(1−t)alpha_B+t alpha_AP].

The corresponding offset is

    beta_blend = [(1−t)alpha_B beta_B+t alpha_AP beta_AP]
                 / [(1−t)alpha_B+t alpha_AP].

Using the published Eq. 19 squares to infer their surface temperatures gives
a required multiplier of about 1.0004–1.0012 on the simple volume-weighted
alpha for the interior points. See `figure4_blend_capacity_audit.csv`.
Thus, conditional on these digitized endpoint fits, the squares are
consistent to graphical resolution with the simple volume-weighted
volumetric heat capacity. A substantially unequal assumed density ratio
under arithmetic mass-weighted density gives a larger multiplier and changes
predicted mixed-composition speeds. This is a sensitivity of the selected
numerical realization, not proof that the authors used an incorrect rule:
separate material parameters are unidentified and endpoint fits have
graphical uncertainty. The article alone does not resolve the ambiguity.
It must be kept separate from numerical discretization error and from the
measured agreement of a chosen homogeneous closure.

The solid curves are harmonic means (Eq. 15) of pure-component speeds at the
same q but distinct surface temperatures. The squares use the geometric
Arrhenius closure (Eq. 19) evaluated at the blend surface temperature from
Eq. 14. The two calculations have different definitions even though they
nearly coincide in the printed figure. The circles and asterisks are resolved
heterogeneous DNS. A planar homogeneous solver study compares a closure with
those DNS results; it does not reproduce their particle morphology.

## Simulation postprocessing

`extract_runs.py` measures each x-column's eta = 0.1, 0.5, and 0.9 crossings
directly from raw AMReX fields. It uses the eta = 0.5 position for a late-time
linear fit, records surface-temperature variation and roughness, and obtains
an independent rate from integrated solid volume. It reports first/second
halves of the late interval to expose transient drift. Regression standard
errors describe fit scatter only, not total simulation uncertainty.

`compare_study.py` reads successful case metadata and raw outputs, rejects
cases that end before their specified duration, caches extracted timeseries,
and writes comparison tables and plots. Analytic targets in metadata are
never used to define measured front positions. The default fit uses the last
40% of available physical time. A final simulation report must discuss
actual convergence and the heating surrogate after those cases are complete.
