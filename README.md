# CausalCompetingRisks

> **Status**: pre-implementation, under active development.
> 
R package for causal inference on **competing-event** survival outcomes, using discrete-time pooled logistic regression. Implements the **separable effects** framework (Stensrud et al. 2020) to decompose a treatment's total effect into a separable **direct** effect (path A → Y) and a separable **indirect** effect (path A → D → Y, through the competing event). Provides parametric g-formula and inverse probability weighting estimators for the cumulative incidence function (CIF) of each event under static, baseline-only treatment regimes, with bootstrap confidence intervals and identifying-assumption accessors.

The treatment `A` is conceptually split into a component `A_Y` acting on the event of interest and a component `A_D` acting on the competing event. The four counterfactual arms `(a_Y, a_D) ∈ {0,1}²` (`arm_11`, `arm_00`, `arm_10`, `arm_01`) yield two algebraically-equivalent decompositions of the total effect (Decomposition A and B).

Second package of a two-package ecosystem: [CausalSurvival](https://github.com/MourerAlex/RCausalSurvival) handles single-event survival and supplies the shared discrete-time machinery (person-time expansion, hazard fitting, IPW/IPCW weights, bootstrap engine) that this package imports and extends.

## Installation

Not yet on CRAN. Development version:

```r
# install.packages("remotes")
remotes::install_github("MourerAlex/RCausalSurvival")        # required dependency
remotes::install_github("MourerAlex/RCausalCompetingRisks")
```

## Usage

```r
library(CausalCompetingRisks)

# --- Real data: Byar & Green (1980) prostate trial (ships with package) -----
data(prostate_data)
# A: 1 = DES, 0 = placebo
# event_type: 0 = censored, 1 = prostate-cancer death (Y), 2 = other-cause death (D)

pt <- to_person_time_competing(prostate_data,
                               id = "id", time = "event_time", event = "event_type",
                               treatment = "A",
                               covariates = c("normal_act", "age_cat", "cv_hist", "hemo_bin"),
                               event_y = 1, event_d = 2, event_c = 0,
                               cut_points = 12)

# --- g-formula --------------------------------------------------------------
fit_g <- causal_competing_risks(pt, method = "gformula")
print(fit_g)
summary(fit_g)
print(assumptions(fit_g))

# --- IPW (separable swap weights, Rep 1 + Rep 2) ----------------------------
fit_i <- causal_competing_risks(pt, method = "ipw", truncate = c(0.01, 0.99))
print(fit_i)

# --- Accessors --------------------------------------------------------------
print(risk(fit_g))             # long-format CIF, one row per (arm, k)
print(contrast(fit_g, method = "gformula"))  # emits the loud ci = NULL warning
print(diagnostic(fit_g))

# --- Bootstrap + contrast with CI -------------------------------------------
# contrast() reports the total effect plus the separable direct and indirect
# effects under Decomposition A and B (RD and RR).
boot_g <- bootstrap(fit_g, n_boot = 500, alpha = 0.05, seed = 1)
print(boot_g)
print(contrast(fit_g, method = "gformula", ci = boot_g))
summary(fit_g, ci = boot_g)

# --- Risk-table accessor ----------------------------------------------------
print(risk_table(fit_g, count = "at_risk"))

# --- Plot -------------------------------------------------------------------
# plot() takes `method` (a fit can hold several). Each fit needs its OWN
# bootstrap — the replicates are method-specific, so reusing one fit's
# bootstrap on another's plot would silently mislabel the bands.
boot_i <- bootstrap(fit_i, n_boot = 500, alpha = 0.05, seed = 1)

plot(risk(fit_g, ci = boot_g), method = "gformula")   # CIF curves + CI ribbons
plot(risk(fit_i, ci = boot_i), method = "ipw_rep1")   # IPW methods: ipw_rep1 / ipw_rep2

# one or more risk-table panels stacked below the curves (cumulative
# events/censoring; at_risk is a snapshot):
plot(risk(fit_g, ci = boot_g), method = "gformula",
     risk_table = c("events_y", "events_d", "censored"))

# plot(contrast(...)) and plot(diagnostic(...)) are planned for a later release.
```

See `vignette("getting-started")` once available.

## Method

- **Estimand**: separable direct and indirect effects on the risk of the event of interest, in the presence of a competing event.
- **Identification**: dismissible component conditions plus the usual exchangeability, positivity, and consistency assumptions (Stensrud et al. 2020; Stensrud & Young 2021). Inspect them with `assumptions(fit)`.
- **Estimation**: discrete-time pooled logistic hazards for the event of interest, the competing event, and censoring; the CIF follows the sequential `C → D → Y` within-interval convention (`F_Y = Σ h_Y (1 − h_D) S(k−1)`).
  
## Notes

A full vignette is planned.

**What are `A_Y` and `A_D`: why four arms?**
The separable-effects framework splits the treatment `A` into a component
`A_Y` acting on the event of interest and `A_D` acting on the competing
event. Setting them independently gives four arms `(a_Y, a_D) ∈ {0,1}²`. The
off-diagonal arms (`arm_10`, `arm_01`) are the counterfactuals where
treatment acts on only one pathway — they are what make the direct/indirect
decomposition possible. The `A_y` / `A_d` columns start equal to the observed
`A`; the estimator overwrites them per arm.

**Decomposition A vs B: which do I report?**
`contrast()` returns both. They are algebraically-equivalent rearrangements of
the same total effect and agree on it exactly; they differ in which arm
anchors the separable direct vs indirect split. Report whichever matches your
scientific question.

**Is it pooled logistic, or a cubic in time?**
Both — same as CausalSurvival. Each hazard (`Y`, `D`, `C`) is a pooled
logistic regression over all person-time rows; time enters as
`k + I(k^2) + I(k^3)` rather than as interval dummies. Default:
`<event> ~ <arm/treatment> + k + I(k^2) + I(k^3) + covariates`. Override any
model with `formulas = list(y = ..., d = ..., c = ..., A = ..., A_num = ...)`.

**Are the IPW weights stabilized? Is truncation on by default?**
Neither by default. The treatment weight is the unstabilized `1 / P(A | L)`
unless you supply a numerator model `formulas$A_num` (e.g. `A ~ 1`), which
switches it to the stabilized ratio. And `truncate` defaults to `NULL` (no
truncation) — pass `truncate = c(0.01, 0.99)` to clip extreme weights
(Cole & Hernán 2008).

**Time-varying treatment / covariates?**
Not in v1. Treatment is **point (baseline)**: `A_0` is carried across all
intervals. Time-varying treatment is planned for a future release.

**Is bootstrap the only way to get confidence intervals?**
In v1, yes (analytic / influence-function CIs are planned). Replicates are
method-specific — a g-formula bootstrap stores g-formula CIFs — so **each fit
needs its own bootstrap**, and reusing one fit's bootstrap on another's plot
would silently mislabel the bands.

## References

- Stensrud MJ, Young JG, Didelez V, Robins JM, Hernán MA (2020). *Separable Effects for Causal Inference in the Presence of Competing Events.* JASA.
- Stensrud MJ, Young JG (2021). *Identification and estimation of separable effects.*

## License

MIT © Alex Mourer
