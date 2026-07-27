# Repository guidance

## Project shape

This repository supports the study described in `README.md`, testing whether ageing or population decline affects national socio-economic performance. It is a data-driven R analysis, not an R package or application:

- `scripts/ageingSEperformanceGH.R` is the sole, sequential analysis script. It defines helper functions, imports and cleans country-level data, derives population and dependency-ratio measures, then runs all analyses and exports figures/CSV results.
- `data/` contains the versioned source datasets. Most tables are country cross-sections or panels; `popXage.csv` is the population-by-age, country-year panel from which population trajectories and age structure are derived.
- `www/` contains only README assets.

The script first combines population, wealth, GDP, inequality, labour productivity, healthy-life-expectancy, wellbeing, HDI/planetary-pressure-adjusted HDI, and other indicator data. It analyses both long-term (1950--2021) and recent (2012--2021) population trajectories. The primary modelling work uses tuned, cross-validated boosted regression trees (BRTs) with repeated regional resampling; later sections add country time-series, regional relationships, and lag scans. Results are written as CSV and JPG files.

## Running and validation

There is no build system, test suite, linter, or single-test command. The smallest non-executing validation is:

```sh
Rscript -e 'parse("scripts/ageingSEperformanceGH.R")'
```

The full script command is:

```sh
Rscript scripts/ageingSEperformanceGH.R
```

Do not run the full script from a fresh clone without first making its paths and inputs available. It uses many hard-coded `setwd()` locations under `~/Documents/Papers/Health/pop trend & wealth/`, writes to that external `out/` directory, and expects additional unversioned inputs such as `areakm2.csv`, `TESpc.csv`, and several intermediate CSVs. Keep the analysis sequential when modifying it because later sections consume objects created earlier.

The script loads `boot`, `dismo`, `gbm`, `ggplot2`, `ggpubr`, `ggrepel`, `ggtext`, `usdm`, `sandwich`, `lmtest`, `car`, and `nlme`. The package list in `README.md` is useful but does not include all packages loaded later in the script.

## Analysis conventions

- Use the three-letter `cntry.code` as the country join key. Retain `cntry`/`country` as labels only, and preserve `year` when handling panel data. `continent.country2.csv` maps these codes to `cont` and the analysis-specific `region`.
- Preserve the script's source-specific cleaning before joins: select the documented/latest observation or time window, remove duplicate country-year records where required, and state whether an `all.x` merge intentionally retains missing outcomes. Do not assume indicators share their most recent year.
- Age columns in `popXage.csv` are `a0` through `a100`. Existing derived measures use these annual age groups to compute total population and the older-age dependency ratio; maintain their definitions and transformation choices when extending models.
- Reuse the BRT helper stack at the top of the script rather than calling `dismo::gbm.step()` directly. `gbm.step.patched` dynamically adds `n.minobsinnode`, and the tuning helpers rank fits by cross-validation correlation, convergence jitter, standard error, and tree-count targets.
- Repeated BRT sections create a `parallel` cluster using `detectCores() - 1`, then stop it explicitly. Follow that lifecycle and avoid nesting parallel work. BRT output naming is based on the response/train-column name plus a period suffix such as `Pred1221` or `Pred5021`.
- The regional and lag analyses are near the end of the script. They share `regional.output.dir` and `lag.output.dir`, use grouped regions to maintain sample sizes, and emit consistently prefixed `regional_depratio_*` and `lag_*` artefacts. Preserve those output contracts if changing the relevant specifications or summary functions.
