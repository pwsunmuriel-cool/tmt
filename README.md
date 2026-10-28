# Continuous-response IRT for a smartphone Trail Making Test

## About this project

This project explores **how to run item response theory (IRT) in R when the responses are continuous**: here, the time taken for each segment of the Trail Making Test (TMT).

To have data to work with, I asked AI assistants (Claude and Gemini) to build a smartphone TMT web app and to simulate 2,000 participants. The R analysis in this repository was also written with AI help.

- You can try the app here: https://pwsunmuriel-cool.github.io/tmt/tmt_mobile_en.html
- A summary of the whole project, including the simulation results: https://pwsunmuriel-cool.github.io/tmt/site/tmt_project_en.html

In this app **each participant gets a randomly generated layout**, so a given transition (e.g. circle 7 → 8) has a different distance, direction and set of nearby circles for each person. The "item" is the position in the sequence (and, in Part B, number→letter vs letter→number), not a fixed stimulus.

## Analyses in `irt_tmt.R`

The script treats the 24 log transition times of each part as continuous items:

| Step | Method | R package |
|---|---|---|
| Distance adjustment | One common within-person slope of log RT on distance (with transition fixed effects, error-free transitions only); slope × centred distance is subtracted before IRT. Each person's overall speed is kept | base R |
| Log RT one-factor model (exploratory) | One-factor model on log RT. Item means = time intensity. The script also reports α = 1 / residual SD (from the uniquenesses), which is van der Linden's (2006) discrimination. The standardized loadings are **not** that discrimination; they relate to the item-specific speed slopes of the extended lognormal model. A dedicated fit (e.g. `LNIRT::LNRT()`) is not done yet | `factanal()` |
| Graded response model | Each transition's time cut into 5 speed bands, then GRM | `mirt` |
| Dimensionality | One speed factor vs two correlated factors (A, B) | `mirt` |
| DIF | Device (iOS vs Android) on the binned data | `mirt` |
| Errors | Rasch model on "error on this transition" (Part B only; Part A too sparse) | `mirt` |
| Explanatory IRT | `log RT ~ part × (distance + position) + switch + (1 + part_B \| person)`; Poisson model for errors | `lme4` |
| Measurement invariance | Multi-group CFA on log RT: configural → metric → scalar → strict, using ΔCFI ≥ −.010 and ΔRMSEA ≤ .015 (Chen, 2007) as reference cut-offs, not as standards validated for this model | `lavaan` |

Not yet done: distribution checks, local dependence (Yen's Q3, with data-specific critical values rather than a fixed .20), person fit for response times, censoring at time limits, a joint speed–accuracy model (e.g. `LNIRT::LNIRT()`), and sensitivity to the number of bins.

## Important: the data are simulated

- **No real participants.** 2,000 simulated people aged 55–89, healthy controls (HC) 1,489, MCI 355, dementia (labelled AD) 156.
- Completion times are calibrated to Tombaugh (2004) paper-and-pencil norms. MCI and dementia rates follow Sun et al. (2014) for Taiwan; ages 55–64 are extrapolated. Errors follow Ashendorf et al. (2008).
- Each person's 24 transition times are built from one latent time per part, so **IRT results on these data largely reflect the simulation's own assumptions** (high reliability, clean fit). They show that the code runs, not that the test is good.

## Results on the simulated data

| | Part A | Part B |
|---|---|---|
| Log RT one-factor model RMSEA / omega | .003 / .94 | .009 / .97 |
| GRM marginal reliability | .94 | .96 |
| GRM discrimination range | 1.33–1.59 | 1.80–2.07 |

- Two speed factors fit far better than one (BIC about 5,850 lower); latent correlation A–B = .70.
- Number→letter transitions take 15.3% longer after adjusting for distance and position.
- No device DIF (0 of 48 transitions); the simulation contains no device effect.
- Measurement invariance: configural, metric and scalar hold across device, age and education; strict fails across age.
- HC vs MCI, AUC: raw Part B time .80, IRT Part B score .77, B−A .77, B/A .64.

## Files

| File | What it is |
|---|---|
| `tmt_mobile_en.html` | Smartphone TMT app. Tap or swipe; exports summary and per-transition CSV |
| `simulate_tmt_en.R` | Simulation |
| `tmt_sim_participants.csv` | One row per participant (2,000 rows) |
| `tmt_sim_transitions.csv` | One row per transition (96,000 rows): `rt_ms`, `distance_norm`, `errors` |
| `plot_tmt.R` → `figures/` | Six descriptive figures |
| `irt_tmt.R` → `irt_output/` | All IRT analyses above; tables, figures and `run_log.txt` |
| `site/tmt_project_en.html` | Project summary page |

## How to run

R 4.6 with `ggplot2`, `dplyr`, `mirt`, `lme4`, `lavaan`.

```r
setwd("path/to/this/folder")
source("simulate_tmt_en.R")   # about 3 seconds
source("plot_tmt.R")
source("irt_tmt.R")           # about 3 minutes; set RUN_DIF <- FALSE to skip the slowest step
```

## Known issues in the app (not yet fixed)

- In swipe mode, crossing an unvisited circle on the way to the correct one counts as an error. Paper-test instructions do not clearly cover this case, so the app's own rule still has to be decided (for example, only a selected circle counts).
- No time limits. Stop rules differ across protocols (for example, CENTER-TBI stops Part A at 100 s and Part B at 300 s); 150 s / 300 s would be this project's own choice.
- Tap vs swipe is not recorded per transition. Tapping may change A, B and B−A compared with drawing, but the direction and size are not established for this interface. The simulation does not build in such an effect.
- Part B uses English letters A–L, which may not suit older Taiwanese adults with little schooling.

## Key references

- Tombaugh, T. N. (2004). Trail Making Test A and B: Normative data stratified by age and education. *Archives of Clinical Neuropsychology, 19*(2), 203–214. https://doi.org/10.1016/S0887-6177(03)00039-8
- Ashendorf, L., et al. (2008). Trail Making Test errors in normal aging, mild cognitive impairment, and dementia. *Archives of Clinical Neuropsychology, 23*(2), 129–137. https://doi.org/10.1016/j.acn.2007.11.005
- Sun, Y., et al. (2014). A nationwide survey of mild cognitive impairment and dementia, including very mild dementia, in Taiwan. *PLoS ONE, 9*(6), e100303. https://doi.org/10.1371/journal.pone.0100303
- Woods, D. L., Wyma, J. M., Herron, T. J., & Yund, E. W. (2015). The effects of aging, malingering, and traumatic brain injury on computerized trail-making test performance. *PLoS ONE, 10*(6), e0124345. https://doi.org/10.1371/journal.pone.0124345 (path length vs segment time; mouse-based computerized TMT, not touch)
- van der Linden, W. J. (2006). A lognormal model for response times on test items. *Journal of Educational and Behavioral Statistics, 31*(2), 181–204. https://doi.org/10.3102/10769986031002181
- van der Linden, W. J. (2007). A hierarchical framework for modeling speed and accuracy on test items. *Psychometrika, 72*(3), 287–308. https://doi.org/10.1007/s11336-006-1478-z
- Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of measurement invariance. *Structural Equation Modeling, 14*(3), 464–504. https://doi.org/10.1080/10705510701301834
- Yen, W. M. (1984). Effects of local item dependence on the fit and equating performance of the three-parameter logistic model. *Applied Psychological Measurement, 8*(2), 125–145. https://doi.org/10.1177/014662168400800201
- Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values for Yen's Q3: Identification of local dependence in the Rasch model using residual correlations. *Applied Psychological Measurement*. https://doi.org/10.1177/0146621616677520

Bibliographic details were checked on 2026-10-28 (Crossref / PubMed, and an independent check with Perplexity); see `Verification_PP.md`.

