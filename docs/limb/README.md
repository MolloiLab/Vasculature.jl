# The leg, carried to the full limb

This folder takes the leg of `docs/notebooks/02_limb_pipeline.jl` past its first stages:

- both trees grown in full to 7 µm terminals;
- the flow at rest and at maximal dilation;
- the contrast transport at rest;
- a narrowed superficial femoral artery;
- the iodine in the whole body.

Nothing outside this folder is changed. `src/`, the heart, lung and leg notebooks in `docs/notebooks/` and the rendered pages are as they were. The docs site does not render this folder.

- `limb_pipeline.jl` is the notebook. It is `02_limb_pipeline.jl` with the changes listed below. It runs the stages that take seconds or minutes and shows what the scripts save.
- The scripts run the long steps, patch by patch, and cache everything in `LIMB_CACHE_DIR`. Each script includes the notebook, so it reuses the notebook's functions and cached stages.

## Licensing

XCAT is licensed. The XCAT leg file, and everything made from it, stay on lab storage. That covers trees, caches, results, figures and pages. None of it belongs in git, an issue or a pull request.

This folder holds code and text only. The scripts read XCAT from the path in `XCAT_LEG_NRB` and write only to `LIMB_CACHE_DIR`.

## Before running

| Variable | What it is |
|---|---|
| `XCAT_LEG_NRB` | XCAT's both-legs NURBS file (licensed, on lab storage) |
| `LIMB_CACHE_DIR` | where every stage, result and figure is written (lab storage or local scratch, never the repository) |
| `LIMB_NOTEBOOK` | optional: another copy of the notebook for the scripts to include |
| `STAGE2_LIMIT_S` | optional: `grow_stage2.jl`'s time limit per patch, in seconds |

The scripts stop if `XCAT_LEG_NRB` or `LIMB_CACHE_DIR` is not set. Set them even where the notebook's defaults would do.

The notebook and the scripts use the docs environment (`docs/Project.toml`). Run `Pkg.instantiate()` there once.

The scripts were run on a Windows workstation with 28 threads and no GPU. The times below are from that machine.
- On Windows, `full_speed()` keeps a process off the efficiency cores. Without it, growth is about 6× slower there. On other systems it does nothing.
- `grow_stage2.jl` stops a stuck worker with `taskkill`, or `kill -9` elsewhere.

## Run order

Run from this folder: `cd docs/limb`.

| Step | Command | Time | Writes to `LIMB_CACHE_DIR` |
|---|---|---|---|
| 1. The notebook: geometry, pass 1, patches, a sample of pass 2 | open `limb_pipeline.jl` in Pluto | minutes the first time | the cached stages (`*.jls`) |
| 2. Pass 2, arteries | `julia grow_stage2.jl arteries` | about 1 h 40 min | `patch_<k>.jls`, `stage2_progress.log` |
| 3. Pass 2, veins | `julia grow_stage2.jl veins` | about 1 h 45 min | `vein_patch_<k>.jls` |
| 4. A look at both trees, with checks | `julia -t auto view_stage2.jl`, then `julia -t auto view_stage2.jl veins`, then `julia view_stage2.jl both` | a few minutes each | `stage2_view*/` |
| 5. Flow, rest and maximal dilation | `julia -t auto flow_leg.jl` | about 20 min | `flow/` |
| 6. Contrast transport at rest | `julia -t auto transport_leg.jl` | about 11 min | `transport_rest/` |
| 7. The transport's page and slide figures | `julia -t auto transport_view.jl`, then `julia transport_figs.jl` | about 1 min each | `transport_rest/` |
| 8. A narrowed superficial femoral artery, 0–95% | `julia -t auto stenosis_leg.jl` | about 2 h | `stenosis/` |
| 9. Its figures | `julia stenosis_figs.jl` | about 1 min | `stenosis/figures/` |
| 10. The iodine in the whole body | `julia -t auto whole_body.jl` | 1–2 min | `whole_body/` |
| 11. Its figures | `julia whole_body_figs.jl` | about 1 min | `whole_body/figures/` |

- Each step skips what is already cached, so a stopped run carries on where it left off.
- `flow_leg.jl test`, `transport_leg.jl test` and `stenosis_leg.jl test` check a few patches and only write to the step's `progress.log`.
- `inputcurve_check.jl` reproduces a precision issue in the package (see *Known issues in the package*). It uses a made-up input, so it does not need XCAT.

## What `limb_pipeline.jl` changes from `02_limb_pipeline.jl`

- **Terminals.** 7 µm terminal arterioles instead of 10 µm; section 3 gives the reasons.
- **Bones.** Five foot bones that XCAT names without a side are added to the bones that block growth.
- **Imaged lumens.** The lumens of the imaged artery and vein take no terminal, and the patches' tissue leaves them out.
- **Veins.** The veins grow in two passes too:
  - Their pass 1 has its own seed (`SEEDS`).
  - Their root is the measured femoral vein's radius (`ROOT_RADII`), because XCAT's femoral vein is narrow.
  - Pass 2 grows each venous patch with `grow_patch(k; system = :veins)`.
- **Patches blocked by bone.** `growth_finds_clear` runs growth's own test. A patch whose feeding segment growth would find blocked by bone keeps that bone in its inlet's lumen, so it can still grow.
- **The foot's maximal flow.** `FOOT_MAX_SHARE`, from the share of muscle below the ankle.
- **`LIMB_FLOW`.** `venous_scale = :grown`.
- **`DEVICE`.** A GPU backend is used only when its hardware is present. Metal also loads on machines that cannot run it.
- **Section 7.** It reads the flow, the transport, the narrowing and the whole body from `CACHE_DIR`. Sections 1–9 are updated to match.

## Copies of package code

The package is not changed. Where the leg needs something it does not do, the script carries a copy:

- `transport_leg.jl` runs the loop of the package's `transport`, with the leg in the body's right-leg place. `ORGAN_SLOTS` takes only `:coronary` and `:lung`.
- `whole_body.jl` copies `Body.simulate` step for step so it can record every segment's iodine. Each run checks the copy against `Body.simulate`.
- `stenosis_leg.jl` takes the narrowing's loss from the package circuit's `stenosis` (Young–Tsai) and applies it on the pass-1 network.

A copy does not follow later changes to the package. Compare these copies whenever those functions change.

## Not done yet


- Compartments, a phantom slab, idealized first-pass perfusion against ground truth, and the CT scan (BasisSimulator.jl).

## Known issues in the package

Found while building the leg. They are not fixed here.

1. **`Flow.InputCurve` loses digits at small k.** Its tables are computed in Float32, which loses digits when k·dt is below about 10⁻³: the leg's slow back-flux from the interstitium. `inputcurve_check.jl` shows it. `transport_leg.jl` uses its own step responses in Float64.
2. **Growth has no early exit when nothing can attach.** An example is a patch whose feeding segment ends in a bone voxel. Growth then runs all its steps, about an hour, and returns the feeding segment alone with only a warning. `grow_stage2.jl` guards against this with a time limit and a check of the grown count.
3. **Growth's voxel lookup truncates coordinates (`itrunc3`).** The voxel it picks can differ from the nearest one, so an end that is in tissue can be in bone for growth. That is why `growth_finds_clear` asks growth itself.
4. **`grow_patch` reports the target terminal count**, not the grown one.
5. **The flow solver can stop at its 60-iteration limit.** At maximal dilation, with a patch's inlet pressure below about 35 mmHg, it stops there with the flow still cycling: per patch by a few percent, per part of the leg by at most 1%. This happened in the stenosis runs.

## Results so far

In words only. The numbers and figures are made from XCAT, so they stay with the results (see *Licensing*).

### What they show

- **Both trees.** Every arterial and venous patch grew its full number of terminals; none failed or ran out of time. The terminal radii are close to their targets, and virtually no terminal ends inside an imaged artery or vein.
- **Flow.** Every check in section 7 passes.
  - The flow at maximal dilation meets its target.
  - At rest, the flow agrees with the leg's oxygen balance and with measured values, and the oxygen left in the femoral vein is on target.
  - No bed's flow reverses, and at rest no vessel is at full tone.
  - At maximal dilation, the flow per gram and the pressure reaching the tissue rise from the foot to the thigh.
- **Contrast at rest.** Contrast reaches the thigh first and the foot last. The femoral vein's outflow peaks well after the arterial input, and much lower. At the end of the run, 120 s after the injection, every part of the leg holds about the same iodine per mL of tissue, still rising slowly, most of it in the capillary beds. What the leg holds matches what the body sends it, within the body model's own bookkeeping.
- **A narrowed superficial femoral artery.**
  - At rest, the flow below the narrowing holds until the narrowing is severe, then falls steeply; the foot loses the most. Once the flow falls, the tissue below runs short of oxygen and its arterioles open fully.
  - At maximal dilation, the flow below falls from the start and is far lower already at moderate narrowing.
  - The thigh above the narrowing keeps its flow, and at maximal dilation takes slightly more.
- **The iodine in the whole body.** The iodine passes the large veins, the heart chambers and the lungs first, then the aorta and the organs. The grown right leg fills slowly, peaks late and keeps nearly all of its iodine to the end. It holds more than the body's lumped left leg, because its exchange with the tissue, its flow and its blood volume are all higher. The liver and gut and the rest of the body are still filling at the end.

### Where the numbers and figures are

In `LIMB_CACHE_DIR`:

| Step | Numbers | Figures |
|---|---|---|
| Both trees | `stage2_view*/summary.txt` | `stage2_view/overview.html` |
| Flow | `flow/summary.txt` | tables in section 7 of the notebook |
| Contrast at rest | `transport_rest/summary.txt` | `transport_rest/figures/figures.html`, `transport_rest/transport_view.html` |
| Narrowed artery | `stenosis/summary.txt` | `stenosis/figures/figures.html` |
| Whole body | `whole_body/summary.txt` | `whole_body/figures/figures.html` |

Section 7 of the notebook also shows the flow, the transport, the narrowing and the whole body, when `LIMB_CACHE_DIR` holds them.

They are on the leg workstation, in its `LIMB_CACHE_DIR` (about 31 GB). They are to be moved to lab storage, at: **(location to be filled in)**.

### How they were made and checked

- **Computed with:** the package at `68c0b75`, as of 2026-10-28.
- **Checked on:** the current `main` (`c1ea484`) on 2026-11-24. The long steps (pass 2, flow, transport, stenosis) were not run again.
  - The notebook loads, with the new parts of section 7.
  - `stenosis_leg.jl test` gives the same numbers as before.
  - `whole_body.jl` writes a byte-identical result, its copy still matches `Body.simulate` exactly, and `whole_body_figs.jl` redraws the same figures.
  - `transport_leg.jl test` gives the same body and patches as before; its differences from the package's own `frame_masses` and `sinus_outflow` moved only in their second digit.
  - The `InputCurve` issue is still there: its tables are still Float32.

