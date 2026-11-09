# Changelog

All notable changes to Vasculature.jl. Each release section links its diff; see RELEASING.md for how a release is cut.

## [Unreleased]

### Changed
- Formatting and lint by [Fatou](https://fatou.dev) with its defaults, replacing Runic: `fatou format .` for the code,
  `docs/fatou_notebooks.sh` for the Pluto notebooks; the `Fatou` CI workflow checks both.
- The docs' notebook export fingerprint no longer covers `src/` or `Project.toml`, so a code change, compat bump or
  release no longer stales the rendered notebook pages; re-render the notebooks a change affects by hand.

## [0.1.0](https://github.com/MolloiLab/Vasculature.jl/releases/tag/v0.1.0) (2026-11-09)

First registered release.

### Added
- Tree growth (`grow`): a port of HCT-CCO onto AcceleratedKernels.jl/KernelAbstractions.jl, running on multithreaded CPUs
  and CUDA, ROCm, oneAPI and Metal GPUs, bit-identical to the upstream C++ on the pinned sample case.
- Blood flow from the pressures (`flow`), contrast transport through a whole-body circulation (`transport`) and voxel
  CT phantoms with their ground-truth perfusion (`phantom`).
- Organ specs for the heart, the lungs and the leg, each with a Pluto notebook pipeline in `docs/notebooks/`.
