# Changelog

All notable changes to Vasculature.jl. Each release section links its diff; see RELEASING.md for how a release is cut.

## [Unreleased]

## [0.1.0](https://github.com/MolloiLab/Vasculature.jl/releases/tag/v0.1.0) (2026-11-09)

First registered release.

### Added
- Tree growth (`grow`): a port of HCT-CCO onto AcceleratedKernels.jl/KernelAbstractions.jl, running on multithreaded CPUs
  and CUDA, ROCm, oneAPI and Metal GPUs, bit-identical to the upstream C++ on the pinned sample case.
- Blood flow from the pressures (`flow`), contrast transport through a whole-body circulation (`transport`) and voxel
  CT phantoms with their ground-truth perfusion (`phantom`).
- Organ specs for the heart, the lungs and the leg, each with a Pluto notebook pipeline in `docs/notebooks/`.
