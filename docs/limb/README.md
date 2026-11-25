# The leg, carried to the full limb

This folder takes the leg of `docs/notebooks/02_limb_pipeline.jl` past its first stages:

- both trees grown in full to 7 µm terminals;
- the flow at rest and at maximal dilation;
- the contrast transport at rest;
- a narrowed superficial femoral artery;
- the iodine in the whole body.

Nothing outside this folder is changed. `src/`, the heart, lung and leg notebooks in `docs/notebooks/` and the rendered pages are as they were. The docs site does not render this folder.

- `limb_pipeline.jl` is the notebook. It is `02_limb_pipeline.jl` with the changes listed below. It runs the stages that take seconds or minutes and shows what the scripts save.
- The scripts run the long steps, patch by patch, and cache everything in `LIMB_CACHE_DIR`. The scripts for the long steps include the notebook, so they reuse its functions and cached stages.

## Quick start: use the saved results

The results are made from XCAT: keep them on lab storage or a local disk, never in git (see *Licensing*).

1. Unzip `limb_results.zip` from share drive (see *Where they are*) to a local disk, with Windows Explorer or 7-Zip. Unzipped, it takes 31 GB.
2. Set `LIMB_CACHE_DIR` to the unzipped `limb_pipeline_cache` folder, and `XCAT_LEG_NRB` to the XCAT leg file (see *The variables*).
3. To look at the results, open these pages in a browser. No Julia is needed.
   - `stage2_view/overview.html`, `stage2_view_veins/overview.html` and `stage2_view_both/overview.html`: the arteries, the veins, and both trees together.
   - `transport_rest/transport_view.html` and `transport_rest/figures/figures.html`: the contrast at rest.
   - `stenosis/figures/figures.html`: the narrowed artery.
   - `whole_body/figures/figures.html`: the whole body.
   - The flow has no page: its numbers are in `flow/summary.txt`, and section 7 of the notebook shows them. Every folder's `summary.txt` holds its step's numbers.
4. To work with them, open `limb_pipeline.jl` (see *Opening the notebook*). Every stage loads from the results instead of being computed again, and section 7 shows the flow, the transport, the narrowing and the whole body.
5. To run a script on them, see *Run order*. Each step skips what is already there.

## Licensing

XCAT is licensed. The XCAT leg file, and everything made from it, stay on lab storage. That covers trees, caches, result files, figures and pages. None of it belongs in git, an issue or a pull request.

This folder holds code and text only. The scripts read XCAT from the path in `XCAT_LEG_NRB` and write only to `LIMB_CACHE_DIR`.

## Before running

### What you need

- **Julia 1.12.** The docs environment (`docs/Project.toml`) was resolved with 1.12.7. Instantiate it once, from the repository root: `julia --project=docs -e "import Pkg; Pkg.instantiate()"`. The notebook and the scripts activate it themselves.
- **A machine like the leg workstation:** Windows, 28 threads, 32 GB of memory, no GPU. The times in *Run order* are from it.
- **About 35 GB of disk** for `LIMB_CACHE_DIR`.
- **On Windows,** `full_speed()` keeps a process off the efficiency cores. Without it, growth is about 6× slower there. `grow_stage2.jl` and `flow_leg.jl`, with the scripts built on it, call it themselves. On other systems it does nothing.
- `grow_stage2.jl` stops a stuck worker with `taskkill`, or `kill -9` elsewhere.

### The variables

| Variable | What it is |
|---|---|
| `XCAT_LEG_NRB` | XCAT's both-legs NURBS file, `both_legs_1.nrb` (licensed). On lab storage: `Yishu\both_legs_1.nrb` on the lab's `Molloilab` share (`Z:\Yishu\both_legs_1.nrb` on the leg workstation). It is read only when the leg is built from scratch: the saved results already hold what the notebook takes from it (`leg_surfaces.jls`, `muscle.jls`) |
| `LIMB_CACHE_DIR` | where every stage, result and figure is written (lab storage or local scratch, never the repository). To start from the saved results: the unzipped `limb_pipeline_cache` (see *Where they are*) |
| `LIMB_NOTEBOOK` | optional: another copy of the notebook for the scripts to include |
| `STAGE2_LIMIT_S` | optional: `grow_stage2.jl`'s time limit per patch, in seconds |

The scripts stop if `XCAT_LEG_NRB` or `LIMB_CACHE_DIR` is not set. Set them even where the notebook's defaults would do. A program sees the variables only if it starts after they are set.

- **Windows,** for every new terminal: `setx LIMB_CACHE_DIR "D:\limb_results\limb_pipeline_cache"` (an example path), and the same for `XCAT_LEG_NRB`. For the current PowerShell only: `$env:LIMB_CACHE_DIR = "D:\limb_results\limb_pipeline_cache"`.
- **macOS and Linux:** `export LIMB_CACHE_DIR=/path/to/limb_pipeline_cache`, and the same for `XCAT_LEG_NRB`, in the shell that starts Julia or in `~/.zshrc`.

### Opening the notebook

- **With Pluto,** from the repository root, in Julia:

  ```julia
  import Pluto   # once, if Pluto is missing: import Pkg; Pkg.add("Pluto")
  Pluto.run(notebook = "docs/limb/limb_pipeline.jl")
  ```

- **With SpaceStation** (the package's *Getting Started*), open the file: `spacestation docs/limb/limb_pipeline.jl`. Opening the folder drops the environment variables, and the notebook then falls back to its defaults: `~/Documents/dev/MolloiLab/xcat-data/legs/both_legs_1.nrb` and `~/xcat-perfused/limb_pipeline_cache`. Unzipping the results into `~/xcat-perfused/` puts them in that default folder, which works either way.

## Run order

Run from this folder: `cd docs/limb`.

| Step | Command | Time | Writes to `LIMB_CACHE_DIR` |
|---|---|---|---|
| 1. The notebook: geometry, pass 1, patches, a sample of pass 2 | open `limb_pipeline.jl` (see *Opening the notebook*) | minutes the first time | the cached stages (`*.jls`) |
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
- From scratch, all the steps take about 6.5 hours on the leg workstation.

## Files in this folder

| File | What it does | Needs XCAT or its results |
|---|---|---|
| `README.md` | this page | no |
| `limb_pipeline.jl` | the notebook: every stage up to a sample of pass 2, and section 7's results | yes |
| `grow_stage2.jl` | pass 2: grows every arterial or venous patch | yes |
| `view_stage2.jl` | draws and checks both trees | yes |
| `flow_leg.jl` | the flow, at rest and at maximal dilation | yes |
| `transport_leg.jl` | the contrast transport at rest, in the body model | yes |
| `transport_view.jl` | builds the transport's page | yes |
| `transport_view_page.jl` | the page's template: draws whatever data it is given | no |
| `transport_figs.jl` | the transport's slide figures | yes |
| `transport_figs_plot.jl` | draws those figures; the other figures share its helpers and colours | no |
| `stenosis_leg.jl` | the narrowed superficial femoral artery | yes |
| `stenosis_figs.jl` | its figures and tables | yes |
| `stenosis_figs_plot.jl` | draws them | no |
| `whole_body.jl` | the iodine in every organ of the body | yes |
| `whole_body_figs.jl` | its figures and tables | yes |
| `whole_body_figs_plot.jl` | draws them | no |
| `inputcurve_check.jl` | shows the package's `InputCurve` precision issue on a made-up input | no |

The drawing files (`*_plot.jl`, `transport_view_page.jl`) only take numbers and draw. A figure can be changed there without running anything else, and made-up numbers can check its layout without XCAT.

## When you change something

Every stage is kept in `LIMB_CACHE_DIR`. Whenever its file is there, it is loaded instead of computed, and nothing checks whether it is out of date. After a change, delete what depends on it, then run again:

| You change | Delete from `LIMB_CACHE_DIR` | Then run again from |
|---|---|---|
| A notebook stage in sections 1–5 | that stage's file and every later one (`leg_surfaces`, `leg`, `muscle`, `inlets`, `pass1`, then `patches` and `vein_patches`), every `patch_<k>.jls` and `vein_patch_<k>.jls`, and every folder | step 1 |
| Pass 2 of one tree (section 6) | that tree's `patch_<k>.jls` or `vein_patch_<k>.jls`, the `stage2_view*` folders, and `flow`, `transport_rest`, `stenosis` and `whole_body` | step 2 or 3 |
| `LIMB_FLOW`, or `flow_leg.jl` | `flow`, `transport_rest`, `stenosis` and `whole_body` | step 5 |
| The injection, or `transport_leg.jl` | `transport_rest` and `whole_body` | step 6 |
| The narrowing, or `stenosis_leg.jl` | `stenosis` | step 8 |
| `whole_body.jl` | `whole_body` | step 10 |
| A figure only | nothing: the figure scripts draw again each time | that figure's step |

## Common problems

| Problem | Why | What to do |
|---|---|---|
| The zip will not open | `unzip` in Git Bash cannot read a zip this large | Use Windows Explorer or 7-Zip |
| Everything is slow | The results are read over the network | Unzip them to a local disk and set `LIMB_CACHE_DIR` there |
| Growth is about 6× slower than in *Run order* | Windows runs the process on its efficiency cores | See `full_speed()` in *What you need*; the scripts warn when it fails |
| The notebook computes stages that should be saved | It does not see your `LIMB_CACHE_DIR`: it started before the variable was set, or SpaceStation opened the folder | Restart the terminal, Pluto or SpaceStation. In SpaceStation, open the file |
| A script stops at once | `XCAT_LEG_NRB` or `LIMB_CACHE_DIR` is not set | Set both, even where the notebook's defaults would do |

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

The main numbers of each step are below. The figures and the result files stay on lab storage (see *Licensing* and *Where they are*).

### Both trees

Every patch grew its full number of terminals; none failed or ran out of time. Virtually no terminal ends inside an imaged artery or vein.

| | Arteries | Veins |
|---|---|---|
| Patches | 1588 | 1715 |
| Terminals | 171.26 million | 171.26 million |
| Segments | 342.5 million | 342.5 million |
| Terminal radius, median (target) | 6.77 µm (7) | 11.88 µm (12.2) |

### Flow

Every check in section 7 passes.

| | Result | Target |
|---|---|---|
| Maximal flow | 3704 mL/min (0.291 mL/min/g) | 3700 mL/min |
| Resting flow | 325 mL/min (0.0256 mL/min/g) | 0.0255 mL/min/g from the oxygen balance; measured 300 ± 30 and 353 mL/min |
| Flow reserve | 11.4 | |
| Oxygen saturation in the femoral vein at rest | 0.662 | 0.662 |
| Oxygen use at rest | 18.6 mL/min | 18–28 mL/min |
| Beds with reversed flow, at rest and at maximal flow | 0 | 0 |
| Patches with a vessel at full tone, at rest | 0 | 0 |

At maximal dilation, from the foot to the top of the thigh, the flow rises from 0.156 to 0.311 mL/min/g, and the pressure at the patches' inlets from 64 to 80 mmHg.

### Contrast at rest

The leg takes 325 mL/min of the body's 4896 mL/min cardiac output (6.65%).

| | Arrives | Peak | Peak at |
|---|---|---|---|
| Input, from the abdominal aorta | 12.4 s | 15.42 mg/mL | 19.3 s |
| Outflow, into the femoral vein | 54.5 s | 1.82 mg/mL | 91.8 s |

Iodine per mL of tissue, in mg/mL. The contrast reaches the thigh first and the foot last:

| Part | Tissue | At 30 s | At 60 s | At 120 s |
|---|---|---|---|---|
| Foot | 899 mL | 0.001 | 0.086 | 0.097 |
| Calf | 1798 mL | 0.025 | 0.094 | 0.098 |
| Knee | 3589 mL | 0.049 | 0.095 | 0.098 |
| Thigh | 6201 mL | 0.060 | 0.095 | 0.099 |

- At 120 s the leg holds 1482 mg of the 15,924 mg injected (9.3%), 1019 mg of it in the capillary beds. That matches what the body sends it, within the body model's own bookkeeping.
- Plasma in the leg: 247 mL. Mean transit time: 82.8 s.

### A narrowed superficial femoral artery

The narrowing is 20 mm long, in the adductor canal. Its middle is 13 cm above the knee joint.

| Narrowing (diameter) | Flow through it, rest | Flow through it, maximal | Pressure below it, rest | Pressure below it, maximal |
|---|---|---|---|---|
| None | 113 mL/min | 1139 mL/min | 89.1 mmHg | 83.8 mmHg |
| 50% | 113 | 775 | 88.0 | 63.2 |
| 70% | 113 | 326 | 79.4 | 36.6 |
| 80% | 108 | 131 | 41.1 | 23.7 |
| 85% | 59 | 61 | 21.6 | 18.3 |
| 90% | 17 | 17 | 14.0 | 14.1 |
| 95% | 1.3 | 1.4 | 11.2 | 11.7 |

- At rest the flow holds until the narrowing is severe, then falls steeply. At maximal dilation it falls from the start.
- The foot loses the most at rest: at 85% it keeps 44% of its flow, the other parts below the narrowing 54–56%.
- At 90%, at rest, the oxygen saturation of the blood leaving the parts below falls from 0.66 to 0.12–0.25, and their arterioles open fully.
- The thigh above the narrowing keeps its flow: 100% at rest, and up to 102% at maximal dilation.

### The iodine in the whole body

15,924 mg of iodine injected into the right arm at rest, as in the body model, with the grown leg as the body's right leg. Iodine as % of the dose:

| | Blood flow (mL/min) | Peak | Peak at | At 120 s |
|---|---|---|---|---|
| Large veins | | 43% | 7.0 s | 1.8% |
| Heart chambers | 4896 | 20% | 10.2 s | 2.5% |
| Lungs | 4896 | 57% | 12.5 s | 8.4% |
| Aorta | | 15% | 18.0 s | 1.5% |
| Brain | 588 | 5.6% | 22.6 s | 0.9% |
| Heart muscle | 162 | 1.6% | 24.7 s | 1.0% |
| Kidneys | 930 | 10% | 26.5 s | 3.0%, and 6.2% in urine |
| Grown right leg | 325 | 9.3% | 86.2 s | 9.3% |
| Liver and gut | 1224 | still rising | | 29.6% |
| Rest of the body | 1177 | still rising | | 25.1% |
| Left leg (the body model's own) | 294 | still rising | | 6.3% |
| Arms | 196 | still rising | | 4.2% |

The grown right leg holds 1.5 times the left leg's iodine at 120 s (1482 against 998 mg). Its flow is about 10% higher, and its exchange with the tissue and its blood volume are higher too.

### Where they are

On lab storage, in one zip: `Yishu\limb_results.zip` on the lab's `Molloilab` share (`Z:\Yishu\limb_results.zip` on the leg workstation). It is 21 GB, and 31 GB unzipped.

To use the results, unzip it and set `LIMB_CACHE_DIR` to the `limb_pipeline_cache` folder inside. The notebook and the scripts then skip every step already done, and section 7 of the notebook shows the flow, the transport, the narrowing and the whole body.

That folder holds every cached stage, both trees among them (`patch_<k>.jls` and `vein_patch_<k>.jls`), and each step's numbers and figures:

| Step | Numbers | Figures |
|---|---|---|
| Both trees | `stage2_view*/summary.txt` | `stage2_view*/overview.html` |
| Flow | `flow/summary.txt` | tables in section 7 of the notebook |
| Contrast at rest | `transport_rest/summary.txt` | `transport_rest/figures/figures.html`, `transport_rest/transport_view.html` |
| Narrowed artery | `stenosis/summary.txt` | `stenosis/figures/figures.html` |
| Whole body | `whole_body/summary.txt` | `whole_body/figures/figures.html` |

### How they were made and checked

- **Computed with:** the package at `68c0b75`, as of 2026-10-28.
- **Checked on:** the current `main` (`c1ea484`) on 2026-11-24. The long steps (pass 2, flow, transport, stenosis) were not run again.
  - The notebook loads, with the new parts of section 7.
  - `stenosis_leg.jl test` gives the same numbers as before.
  - `whole_body.jl` writes a byte-identical result, its copy still matches `Body.simulate` exactly, and `whole_body_figs.jl` redraws the same figures.
  - `transport_leg.jl test` gives the same body and patches as before; its differences from the package's own `frame_masses` and `sinus_outflow` moved only in their second digit.
  - The `InputCurve` issue is still there: its tables are still Float32.

