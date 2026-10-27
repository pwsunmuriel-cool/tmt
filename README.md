# Continuous-response IRT for a smartphone Trail Making Test

## What I want to ask

On Wednesday, October 21, I attended your INS online lecture, and I was too shy to ask my question during the Q&A. Here is my first question: how to do it with continuous variables. So I started a small project with my AI assistants.

More precisely, I want to know **how to run item response theory (IRT) in R when the responses are continuous**: here, the time taken for each segment of the Trail Making Test (TMT).

To have data to work with, I asked AI assistants (Claude and Gemini) to build a smartphone TMT web app and to simulate 2,000 participants. The R analysis in this repository was also written with AI help.

- You can try the app here: https://pwsunmuriel-cool.github.io/tmt/tmt_mobile_en.html
- A summary of the whole project, including the simulation results: https://pwsunmuriel-cool.github.io/tmt/site/tmt_project_en.html

With my AI assistants' help, here are my two questions.

### Question 1. Random layouts: can a transition still be treated as the same item?

On the paper TMT the layout is fixed, so "transition 7" (circle 7 → 8) is the identical jump for everyone. In this app **each participant gets a randomly generated layout**, so transition 7 has a different distance, direction and set of nearby circles for each person. The "item" is really the position in the sequence (and, in Part B, number→letter vs letter→number), not a fixed stimulus.

- Is it acceptable to treat the 24 transitions of each part as items under random layouts?
- Currently the analysis removes the distance between circles (a within-person regression on log time) before fitting IRT. Crowding is not adjusted. Is that enough?
- Or should the app switch to fixed layouts? Options I am considering:
  1. one fixed layout for everyone, like the paper test;
  2. a small set of fixed layouts (e.g. 4 forms with matched path length, no crossing lines), rotated across participants and retests;
  3. keep random layouts and model layout features (distance, crowding) explicitly.

### Question 2. Are the methods the AI used appropriate?

The script `irt_tmt.R` treats the 24 log transition times of each part as continuous items:

| Step | Method | R package |
|---|---|---|
| Distance adjustment | Within-person regression of log RT on distance, subtracted before IRT | base R |
| Lognormal response-time model (van der Linden, 2006) | One-factor model on log RT: loadings = discrimination, means = time intensity | `factanal()` |
| Graded response model | Each transition's time cut into 5 speed bands, then GRM | `mirt` |
| Dimensionality | One speed factor vs two correlated factors (A, B) | `mirt` |
| DIF | Device (iOS vs Android) on the binned data | `mirt` |
| Errors | Rasch model on "error on this transition" (Part B only; Part A too sparse) | `mirt` |
| Explanatory IRT | `log RT ~ part × (distance + position) + switch + (1 + part_B \| person)`; Poisson model for errors | `lme4` |
| Measurement invariance | Multi-group CFA on log RT: configural → metric → scalar → strict, judged by ΔCFI ≥ −.010 and ΔRMSEA ≤ .015 | `lavaan` |

Specific doubts:

- Is a one-factor model on log RT a sound way to fit the lognormal response-time model, or should I use a dedicated package (e.g. `LNIRT`) or a Bayesian model?
- **Two-step vs one-step distance handling:** is subtracting the distance effect first acceptable, or should distance always be a covariate inside the model?
- **Binning for the GRM:** the GRM score separates MCI from healthy controls slightly worse than raw time (AUC .77 vs .80). Should the continuous model be the primary analysis and the GRM only a supplement?
- **Measurement invariance at N = 2,000:** are ΔCFI / ΔRMSEA cut-offs appropriate here? Strict invariance fails across age because older people's times vary more between transitions. I read this as within-person variability, not bias. Is that interpretation right?
- **What is missing?** Not yet done: distribution checks, local dependence (Yen's Q3), person fit for response times, censoring at time limits (150 s / 300 s), a joint speed–accuracy model, and sensitivity to the number of bins.

## Important: the data are simulated

- **No real participants.** 2,000 simulated people aged 55–89, healthy controls (HC) 1,489, MCI 355, dementia (labelled AD) 156.
- Completion times are calibrated to Tombaugh (2004) paper-and-pencil norms. MCI and dementia rates follow Sun et al. (2014) for Taiwan; ages 55–64 are extrapolated. Errors follow Ashendorf et al. (2008).
- Each person's 24 transition times are built from one latent time per part, so **IRT results on these data largely reflect the simulation's own assumptions** (high reliability, clean fit). They show that the code runs, not that the test is good.

## Results on the simulated data

| | Part A | Part B |
|---|---|---|
| Lognormal model RMSEA / omega | .003 / .94 | .009 / .97 |
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

- In swipe mode, crossing an unvisited circle on the way to the correct one counts as an error; standard scoring does not count this.
- No time limits (150 s for Part A, 300 s for Part B).
- Tap vs swipe is not recorded per transition; smartphone tapping is known to shrink B−A compared with drawing.
- Part B uses English letters A–L, which may not suit older Taiwanese adults with little schooling.

## Key references

- Tombaugh, T. N. (2004). Trail Making Test A and B: Normative data stratified by age and education. *Archives of Clinical Neuropsychology, 19*(2), 203–214. https://doi.org/10.1016/S0887-6177(03)00039-8
- Ashendorf, L., et al. (2008). Trail Making Test errors in normal aging, mild cognitive impairment, and dementia. *Archives of Clinical Neuropsychology, 23*(2), 129–137. https://doi.org/10.1016/j.acn.2007.11.005
- Sun, Y., et al. (2014). A nationwide survey of mild cognitive impairment and dementia, including very mild dementia, in Taiwan. *PLoS ONE, 9*(6), e100303. https://doi.org/10.1371/journal.pone.0100303
- Woods, D. L., et al. (2015). *PLoS ONE*. https://doi.org/10.1371/journal.pone.0124345 (distance and segment timing in a computerized TMT; not checked against the source)
- van der Linden, W. J. (2006). A lognormal model for response times on test items. *Journal of Educational and Behavioral Statistics* (not checked against the source)
- van der Linden, W. J. (2007). A hierarchical framework for modeling speed and accuracy on test items. *Psychometrika* (not checked against the source)

