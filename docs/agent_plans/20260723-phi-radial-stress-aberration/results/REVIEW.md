# Adversarial review

No production implementation was selected, so the implementation-specific
review items (event ordering, duplicate solves, device captures, and feedback
into physics) are not applicable. The analysis was reviewed for ways it could
falsely dismiss the visible feature.

1. **Could cubic detrending hide a real oscillation?** The independent
   monotone-envelope metric is zero on all native and saved profiles. The
   manufactured alternating-mode test recovers `1234.5 Pa` and `4321 Pa` to
   floating-point accuracy, while its clean mode is `7.2e-11 Pa`.
2. **Could polar projection be wrong?** On the symmetry axes, radial stress is
   exactly the coordinate-normal component, so no off-axis co-location is
   required for the primary result. The full saved polar field agrees with
   the exact face-average reconstruction to `1.88e-5 Pa` after synchronization.
3. **Could `eta` pressure balancing explain the band?** No. Both axis bands
   have `eta_min = eta_max = 1.0` and zero neighbor delta. Adding `p eta I`
   is only a constant offset there.
4. **Could the controlled solve be a no-op?** No. Every displacement and
   stress hash changes, while every frozen non-mechanics hash and the hierarchy
   remain identical.
5. **Could the residual comparison be overstated?** Yes if the offline
   band `Linf` were compared directly with the Newton composite norm. The
   result reports them separately and uses the paired radial profiles—not
   residual equality—to reject the visual fix.
6. **Could the conclusion be generalized beyond the sampled axes?** The
   native-face statement is limited to the two symmetry axes. The conclusion
   about output reconstruction over the arc is supported separately by the
   full finest-level saved-versus-face-average comparison. A genuinely native
   off-axis radial traction is not available on a Cartesian face without an
   explicitly labeled reconstruction.

No finding changes the decision-table outcome.

