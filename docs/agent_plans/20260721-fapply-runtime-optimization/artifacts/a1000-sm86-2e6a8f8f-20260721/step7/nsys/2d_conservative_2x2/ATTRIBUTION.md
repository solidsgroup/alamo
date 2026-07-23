# Step 7 trace attribution

The trace SQLite `ProcessStreams` table embeds the exact target executable path:

`docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/bin/alamo_gpu-2d-fast-sm86-baseline`

That file is read-only, has build timestamp `2026-07-21 15:29:24 -0500`, and its
SHA-256 is `ee2d5669c904620e9bc0c1cee596375dd8459843b20ec3991a5bc2033877c4da`,
matching `baseline/BINARY_SHA256SUMS`. The trace session began at
`2026-07-22T21:37:33Z`; embedded target GPU identity is NVIDIA RTX A1000 UUID
`ff00e057-b36d-c833-9da1-8f71516e7d73`.

The input SHA-256 is
`5bc70d56a2fa7443096eef4659d9685bb41f5ca91da7824426488f016730d022`.
The trace embeds the process command and arguments, including
`tests/ElasticSoftVoid/input`, `max_step=2`, the 64x64 two-level conservative
SoftVoid overrides, `elastic.solver.pre_smooth=2`, and
`elastic.solver.post_smooth=2`. The plot metadata independently records the
same effective two-step/2x2 configuration.

Artifact hashes:

- `trace.nsys-rep`: `cc4192acc2febef0ceda80137f751ebae05f0072bffc8b053d282ff6dee576c3`
- `trace.sqlite`: `a425f79913f9caebb8a6fd8f93a7caa94790596b6730a3e3df10942c6da172f8`
- `stats_cuda_api_sum.csv`: `9d1f158bbf35976d0563260676264c5702b7660140378eb55f5344197d2aa8dc`
- `stats_cuda_api_trace.csv`: `06d42f0213a055ddcb6d67a2515665c2f115abb4fea985c01364bee3c7a16c1f`
- `stats_cuda_gpu_kern_sum.csv`: `28dbeb65305c786e50b333f8612d2deb085678d17b44b22756ed1e09a69d169e`
- `stats_cuda_gpu_trace.csv`: `6a4a34eb295d68057201e22dd8ee28d422a8d3e5b53652da909244f72bb86cb0`
