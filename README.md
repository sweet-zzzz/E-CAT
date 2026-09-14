# E-CAT

Standalone R implementation of E-CAT for instrument selection and estimation
of a constant treatment effect in the presence of candidate invalid and weak
instruments. This release contains the joint distance-rank and residual-moment
implementation used by the original E-CAT comparison experiments.

## Repository structure

```text
example_data/
    simulation_01.csv
    simulation_02.csv
utils/
    ecat_core.R
E_CAT.R
README.md
example.R
```

## Requirements

Tested with R 4.4.1. The method requires `energy`. The example checks for it
and installs it from CRAN if missing (internet access is needed only for
installation). RStudio editor-path detection uses `rstudioapi`, which the
example also installs if missing in an interactive RStudio session.
No manual `library()` calls are needed. For manual installation:

```r
install.packages("energy")
```

## Quick start

### RStudio: open and click Run

Keep `example.R`, `E_CAT.R`, and both subfolders together in the release folder.
Open the saved `example.R` in RStudio. Place the cursor inside its outer
`{ ... }` block, clear any text selection, and click **Run** once. The entire
example is one expression, so both designs run together. **Source** and
selecting all code followed by **Run** are also supported.

The script detects the saved editor file's folder using the
[RStudio editor API](https://rstudio.github.io/rstudioapi/reference/rstudio-editors.html).
Changing the working directory is unnecessary in RStudio. If editor detection
is unavailable, set the working directory to the release folder and run again.
Opening a file alone does not load its functions; the example loads `E_CAT.R`
and its numerical routines automatically.

The console prints only the seven-column table described below.
`ecat_summary` contains exactly those seven columns,
and `ecat_results` contains both complete fits, including candidate
scores. For example:

```r
ecat_results[["simulation_02.csv"]]$joint_info$candidates
```

### Command line

From the repository directory:

```sh
Rscript example.R
Rscript E_CAT.R
Rscript E_CAT.R example_data/simulation_02.csv 3
```

`example.R` runs the first dataset with K = 2 and the second with K = 3,
and prints the seven-column evaluation table.
Candidate scores remain available in `ecat_results`. `E_CAT.R` accepts an optional CSV path
and subset size K. Without arguments it uses the first dataset and K = 2.
Scripts also work when invoked by their full path from another directory.
Relative input CSV paths are resolved against the current working directory.

In R or RStudio:

```r
source("E_CAT.R")
data <- read.csv("example_data/simulation_01.csv")
result <- ecat(data, K = 2)
result$selected_instruments
result$beta_hat
result$confidence_interval
```

## Main function

```r
ecat(data, instruments = grep("^IV", names(data), value = TRUE),
     treatment = "Treatment", outcome = "Outcome",
     covariates = character(0), K = 2,
     lambda_moment = 1, tuning.1st = 1, alpha = 0.05,
     max_combs = 2000, verbose = FALSE)
```

Input is a data frame with one observation per row. Selected columns must be
numeric, finite, and distinct. Resolve missing values before calling `ecat`.
The instrument matrix must have full column rank after covariate adjustment.

| Argument | Meaning |
| --- | --- |
| `instruments` | Names of candidate instrument columns, in the desired order. |
| `treatment`, `outcome` | Names of the treatment and outcome columns. |
| `covariates` | Optional numeric adjustment columns; default: none. |
| `K` | Target subset size, an integer from 2 to the number of candidates. |
| `lambda_moment` | Nonnegative weight on the residual-moment rank. |
| `tuning.1st` | Positive first-stage threshold multiplier; experiment default: 1. |
| `alpha` | Nominal interval significance level; default: 0.05. |
| `max_combs` | Maximum number of subsets retained for joint scoring. |
| `verbose` | Print intermediate selection details. |

Custom column names and covariates can be supplied explicitly:

```r
# my_data contains numeric columns dose, response, z_a, z_b, z_c, and age.
result <- ecat(my_data, treatment = "dose", outcome = "response",
               instruments = c("z_a", "z_b", "z_c"),
               covariates = "age", K = 2)
```

Main return fields:

| Field | Meaning |
| --- | --- |
| `beta_hat` | Estimated constant treatment effect. |
| `selected_instruments` | Selected instrument column names. |
| `screened_instruments` | Column names retained by the first-stage screen. |
| `confidence_interval` | Nominal interval from the original variance estimator. |
| `Valid_IV_subset`, `Strong_IV_set` | One-based indices into the input instrument list. |
| `joint_info$candidates` | Candidate subsets, names, raw scores, ranks, and joint score. |
| `estimate_detail` | Original effect, variance, and interval output. |
| `n` | Number of observations used. |
| `runtime_sec`, `cpu_time_sec` | Elapsed and CPU seconds for the E-CAT call. |

## Evaluation metrics

The printed table and `ecat_summary` contain only these columns, in this order:

| Field | Definition |
| --- | --- |
| `filename` | Input dataset filename. |
| `n` | Sample size. |
| `K` | Target instrument subset size. |
| `beta_hat` | Estimated treatment effect. |
| `true_hat` | Known true treatment effect, equal to 1 for both examples. |
| `mse` | Mean squared estimation error within the design. |
| `mean_runtime_sec` | Mean elapsed E-CAT runtime in seconds. |

Each design contains one fixed sample, so `mse = (beta_hat - true_hat)^2`
and `mean_runtime_sec` is the runtime of that single fit. These are not
Monte Carlo averages over independently generated samples. `true_hat` is
the requested column name for the known truth, not an estimated quantity.
Timing excludes package installation, file loading, and console printing.
Runtime varies with the computer and session.

Candidate `subset` indices refer to the screened instrument pool;
`instrument_names` maps them back to input names. `joint_info` is `NULL` when
the screen retains at most K instruments and no subset ranking is needed.

## Method

1. Optionally regress the treatment, outcome, and instruments on covariates
   with an intercept and use their residuals.
2. Screen instruments using the original first-stage coefficient threshold.
3. Compute pairwise auxiliary-residual distance correlations and average them
   within each candidate subset.
4. Compute standardized residual-moment scores using the original selected-IV
   effect estimator. Rank each score across candidate subsets and minimize
   `distance_rank01 + lambda_moment * moment_rank01`.
5. Estimate the effect with the selected subset and the original covariance
   weighting routine. The TSHT-named internal routines support E-CAT screening
   and estimation; a standalone TSHT comparison method is not included.

The numerical core is extracted from `Main_New_functions.R`, specifically
`ECAT_JointRank_TSHT` and its dependencies. The public interface adds input
validation, portable loading, and original-column-name mapping.

If all instruments fail the screen, the original algorithm warns and falls
back to all candidates. If at most K survive, it uses all survivors, possibly
fewer than K. Tied joint scores select the first subset in enumeration order.
If more than `max_combs` subsets exist, only the lowest distance-score subsets
receive joint scoring. All combinations are still enumerated first, so this
limit does not prevent combinatorial memory growth. Without covariates the
original estimator does not add an intercept or automatically center inputs.

## Example data

Both files are synthetic demonstration data generated from the existing
`make_configurable_iv_data` experiment generator with its default structural
parameters and Gaussian noise. Only the following two designs are included,
with one fixed sample per design. The true constant treatment effect is 1.

| Dataset | Observations | Strong valid IVs | Weak valid IVs | Invalid IVs | K | Seed |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `simulation_01.csv` | 200 | 2 | 2 | 2 | 2 | 1 |
| `simulation_02.csv` | 500 | 4 | 3 | 7 | 3 | 2 |

The generator uses nonlinear instrument contributions, a shared unobserved
confounder, and direct outcome contributions from invalid instruments.
The files contain observed variables only; no personal data are included.

Both datasets contain `Treatment` and `Outcome`. Instrument columns are:

| Dataset | Strong valid | Weak valid | Invalid |
| --- | --- | --- | --- |
| `simulation_01.csv` | IV1-IV2 | IV3-IV4 | IV5-IV6 |
| `simulation_02.csv` | IV1-IV4 | IV5-IV7 | IV8-IV14 |

Expected results using each design's K, rounded to six decimals:

| Dataset | Selected instruments | Estimated effect |
| --- | --- | ---: |
| `simulation_01.csv` | IV5, IV6 | 1.889287 |
| `simulation_02.csv` | IV1, IV3, IV4 | 0.623103 |

These are individual execution examples rather than an aggregate benchmark.
The first sample selects invalid instruments; the second selects three strong
valid instruments. Selection does not certify instrument validity, and the
nominal interval does not account for instrument selection uncertainty.
