# NOVA/A100 retained-2/2 confirmation

Date: 2026-07-22

Device: NVIDIA A100 80GB PCIe, sm_80,
`GPU-c6a01dea-8df2-84e4-2e50-7bde877e70cc`

Build job `11751553` completed in 8:24 with exit `0:0`. The exact binary is
`bin/alamo_gpu-3d-profile-cuda80-g++`, SHA-256
`036c7ec520c6ae301b31d1e75a107b5dff2e343f02729c880f7be3ceeda98c84`, at
HEAD `2e6a8f8f5430a58dd6a85b07b08d472e2bd1d5cd`.

Warmups are excluded below. Values are median ± median absolute deviation over
five matched measurements per arm.

| metric | 4/4 | 2/2 | reduction |
|---|---:|---:|---:|
| external wall | 15.58 ± 0.02 s | 13.00 ± 0.01 s | 16.560% |
| inclusive MLMG solve | 11.42 ± 0.01 s | 8.846 ± 0.002 s | 22.539% |
| inclusive FApply | 8.388 ± 0.008 s | 6.434 ± 0.001 s | 23.295% |
| FApply calls | 15,477 | 12,660 | 18.201% |

Timing job `11751561` completed all warmups and five alternating measurements,
then failed only when the first version of the campaign-local exact validation
wrapper called a newer helper API absent from the isolated checkout. Validation
tail job `11751579` reused the preserved timing directory, ran both physics
bundles, and completed in 1:53 with exit `0:0` and marker
`NOVA_A100_VALIDATION_TAIL_PASS`.

The 4/4-versus-2/2 physics comparison has overall and gate verdicts `PASS` for
the frozen `max_step=2` case. The final residual is `4.459071293e-09` for 4/4
and `9.721808305e-09` for 2/2, both within the `1e-08` absolute gate. This does
not establish longer-horizon stability.

Nsight Compute job `11751588` captured one retained-2/2 FApply launch with
application replay and completed with exit `0:0`. The 256-thread launch used
254 registers/thread, was limited to one block/SM by registers, and measured
12.50% theoretical / 12.02% achieved occupancy. Default kernel-replay attempts
`11751583` (32 GiB) and `11751584` (64 GiB) hit Slurm host-memory cgroups while
snapshotting the managed arena; neither is accepted profile evidence.

Primary raw evidence is under `nova-results-11751561/`. The compact local copy
includes all timing logs, commands, manifests, field norms, metrics, validation
reports, profiler report/CSV, and job logs. It intentionally omits the 3 GB of
raw `000*node`/`000*cell` validation plot payloads retained on NOVA. Of the
remote artifact manifest entries, all 103 copied files re-hash correctly; the
72 omitted plotfile entries cannot be revalidated from the compact copy.
