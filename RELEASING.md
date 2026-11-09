# Releasing Vasculature.jl

A release is one pull request and one registry comment. The docs deploy from `main` on their own.

## Rules that keep `main` releasable

- **`version` in `Project.toml` stays at the last registered version** until a release pull request changes it. The
  `version` job in `.github/workflows/CI.yml` fails any commit whose `version` is neither that nor a valid next one.
- **Every user-visible change adds a line under `## [Unreleased]`** in `CHANGELOG.md`, in the same pull request.

## Cutting a release

1. **Choose the version.** Below 1.0, a breaking change bumps the minor version (0.1 → 0.2) and anything else the patch
   (0.1.0 → 0.1.1). Breaking means removing or renaming an export, changing what an argument means, or changing the
   numbers a caller gets back unasked.
2. **Branch** `release/<version>` from an up-to-date `main`.
3. **Bump** `version` in `Project.toml` and `version:` in `CITATION.cff`.
4. **Date the changelog:** rename `## [Unreleased]` to
   `## [<version>](https://github.com/MolloiLab/Vasculature.jl/compare/v<previous>...v<version>) (<YYYY-MM-DD>)` and add
   an empty `## [Unreleased]` above it. Group entries as `Breaking` / `Added` / `Changed` / `Fixed` / `Performance`.
5. **Test:** `julia --project=. -t 8 -e 'using Pkg; Pkg.test()'`, the GPU paths on a GPU (CI has none), and loading on
   the newest released Julia, because the registry's AutoMerge loads the package there whatever the `julia` compat says.
6. **Merge the pull request, then comment on its merge commit** on GitHub:

   ```
   @JuliaRegistrator register
   ```

   A breaking release must add release notes, or AutoMerge blocks it:

   ```
   @JuliaRegistrator register

   Release notes:

   ## Breaking changes
   - <the Breaking entries of the changelog section, one line each>
   ```

   When AutoMerge merges (about 15 minutes; a package's **first** registration waits 3 days), TagBot creates the
   `v<version>` tag and the GitHub Release.

## If something goes wrong

- **"Not a valid version increment":** `version` skipped a number. Set it to the next valid one and comment again.
- **AutoMerge fails to load the package:** a dependency does not resolve on the newest Julia. Fix the compat bound (or
  release the dependency) and comment again.
- **TagBot did not tag:** run it from the Actions tab (`TagBot`, "Run workflow").
