# =============================================================================
# Trail Making Test: item response theory (IRT) analysis of Parts A and B
#
#   Items   = the 24 transitions of each part (A01..A24, B01..B24)
#   Persons = participants
#
#   1. Prepare log RTs and remove the distance effect (layouts are random)
#   2. Lognormal response-time model (van der Linden, 2006) per part
#   3. Graded response model (GRM) on binned speed, per part
#   4. Dimensionality: are A and B one speed trait or two?
#   5. Differential item functioning (DIF) by device
#   6. Accuracy: Rasch model on errors (Part B)
#   7. Validity of the IRT scores against total time and diagnosis
#   8. Explanatory IRT with lme4: distance, position and switch effects
#   9. Measurement invariance on continuous log RT with lavaan
#
#   Input : tmt_sim_participants.csv, tmt_sim_transitions.csv (or app exports)
#   Output: ./irt_output/*.csv and *.png
#   Packages: mirt, ggplot2, dplyr; lme4 for section 8 and lavaan for section 9
#             (each skipped if missing)
#   Real data: aim for N >= 500 before trusting GRM item parameters.
# =============================================================================

library(mirt)
library(ggplot2)
library(dplyr)

N_CAT   <- 5        # speed categories per transition in the GRM
RUN_DIF <- TRUE     # DIF by device takes a few minutes
DIF_VAR <- "device" # any two-level column in the participants file
OUT     <- "irt_output"
EM_CTRL <- list(NCYCLES = 3000)  # mirt's default 500 EM cycles can stop short
# "Observed log-likelihood may still be increasing" with changes of ~0.001 at the
# end of a fit is harmless; the same warning during M2()/itemfit() comes from the
# fit statistics, not from the models whose parameters are reported.
dir.create(OUT, showWarnings = FALSE)
set.seed(20261024)

p  <- read.csv("tmt_sim_participants.csv")
tr <- read.csv("tmt_sim_transitions.csv")
tr <- tr[tr$participant_id %in% p$participant_id, ]

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(), legend.position = "top"))
part_cols <- c(A = "#0072B2", B = "#D55E00")
save_png <- function(plot, name, w = 8, h = 4.5) {
  ggsave(file.path(OUT, name), plot, width = w, height = h, dpi = 200, bg = "white")
}

# ---- 1. Prepare: log RT adjusted for distance ---------------------------------
# Layouts are random, so transition k connects different circles for different
# people. Item identity is "position in the sequence"; distance between the two
# circles is a nuisance that we remove with a within-person (fixed-effects) slope.
adjust_distance <- function(d) {
  d$log_rt <- log(d$rt_ms)
  ok <- d$errors == 0                          # error transitions include correction time
  y  <- d$log_rt        - ave(d$log_rt,        d$participant_id)
  x  <- d$distance_norm - ave(d$distance_norm, d$participant_id)
  slope <- coef(lm(y ~ x + factor(d$transition_index), subset = ok))["x"]
  d$log_rt_adj <- d$log_rt - slope * (d$distance_norm - mean(d$distance_norm))
  attr(d, "slope") <- unname(slope)
  d
}

to_wide <- function(d, value, part) {
  m <- tapply(d[[value]], list(d$participant_id, d$transition_index), function(v) v[1])
  colnames(m) <- sprintf("%s%02d", part, as.integer(colnames(m)))
  m[as.character(p$participant_id), , drop = FALSE]
}

parts <- c("A", "B")
long  <- lapply(parts, function(pt) adjust_distance(tr[tr$part == pt, ]))
names(long) <- parts
Y   <- lapply(parts, function(pt) to_wide(long[[pt]], "log_rt_adj", pt)); names(Y)   <- parts
ERR <- lapply(parts, function(pt) to_wide(long[[pt]], "errors",     pt)); names(ERR) <- parts

cat("Distance slope (log ms per unit of board diagonal):",
    sprintf("A = %.2f, B = %.2f", attr(long$A, "slope"), attr(long$B, "slope")), "\n")

# Letter targets in Part B (odd transitions end on a letter: 1->A, 2->B, ...)
to_letter <- function(k) k %% 2 == 1

# ---- 2. Lognormal response-time model -------------------------------------------
# log T_ik = beta_k - lambda_k * tau_i + e_ik,  e_ik ~ N(0, sigma_k^2)
#   beta_k  : time intensity of transition k (mean log time)
#   lambda_k: loading (how strongly k reflects the person's speed tau)
#   alpha_k : 1 / sigma_k, van der Linden's discrimination
# A one-factor linear model on log RT with free loadings is the extended model
# (item-specific speed slopes lambda_k); van der Linden (2006) fixes lambda_k = 1.
# So std_loading is not alpha; alpha is computed from the uniquenesses below.
lognormal_fit <- function(Ymat, part) {
  fa  <- factanal(Ymat, factors = 1, scores = "regression")
  sdv <- apply(Ymat, 2, sd)
  std <- as.numeric(fa$loadings)
  sgn <- sign(sum(std))                        # factor = slowness
  std <- std * sgn
  n   <- nrow(Ymat); chi <- unname(fa$STATISTIC); df <- unname(fa$dof)
  items <- data.frame(
    part = part, item = colnames(Ymat), transition = seq_len(ncol(Ymat)),
    time_intensity_ms = round(exp(colMeans(Ymat))),
    std_loading = round(std, 3),
    loading_log = round(std * sdv, 3),
    alpha = round(1 / (sqrt(fa$uniquenesses) * sdv), 2)
  )
  fit <- data.frame(
    part = part, chisq = round(chi, 1), df = df,
    rmsea = round(sqrt(max((chi - df) / (df * (n - 1)), 0)), 3),
    omega = round(sum(std)^2 / (sum(std)^2 + sum(fa$uniquenesses)), 3)
  )
  list(items = items, fit = fit, speed = -sgn * as.numeric(fa$scores))
}

LN <- lapply(parts, function(pt) lognormal_fit(Y[[pt]], pt)); names(LN) <- parts
ln_items <- bind_rows(LN$A$items, LN$B$items)
ln_fit   <- bind_rows(LN$A$fit,   LN$B$fit)
write.csv(ln_items, file.path(OUT, "lognormal_items.csv"), row.names = FALSE)
cat("\n== Lognormal RT model: fit ==\n"); print(ln_fit, row.names = FALSE)

# ---- 3. Graded response model on binned speed --------------------------------
# Each transition's adjusted log RT is cut into N_CAT bins with cut points pooled
# across the part's transitions, so slow transitions really are "harder" (higher
# thresholds). Categories are reversed: higher = faster, so theta = speed.
bin_speed <- function(Ymat) {
  cuts <- quantile(Ymat, probs = seq(0, 1, length.out = N_CAT + 1), na.rm = TRUE)
  cuts[c(1, length(cuts))] <- c(-Inf, Inf)
  out <- apply(Ymat, 2, function(x) N_CAT - as.integer(cut(x, cuts, labels = FALSE)))
  dimnames(out) <- dimnames(Ymat)
  out
}
CAT <- lapply(Y, bin_speed)

grm_items <- function(mod, part) {
  cf <- as.data.frame(coef(mod, IRTpars = TRUE, simplify = TRUE)$items)
  cbind(part = part, item = rownames(cf), transition = seq_len(nrow(cf)), round(cf, 3))
}

GRM <- lapply(parts, function(pt) mirt(CAT[[pt]], 1, itemtype = "graded", verbose = FALSE, technical = EM_CTRL))
names(GRM) <- parts

grm_par <- bind_rows(grm_items(GRM$A, "A"), grm_items(GRM$B, "B"))
grm_fit <- bind_rows(lapply(parts, function(pt) {
  m2  <- M2(GRM[[pt]], type = "C2")
  itf <- itemfit(GRM[[pt]])
  data.frame(part = pt, M2 = round(m2$M2, 1), df = m2$df, RMSEA = round(m2$RMSEA, 3),
             CFI = round(m2$CFI, 3), SRMSR = round(m2$SRMSR, 3),
             misfit_items_p01 = sum(itf$p.S_X2 < .01, na.rm = TRUE),
             marginal_rxx = round(marginal_rxx(GRM[[pt]]), 3))
}))
write.csv(grm_par, file.path(OUT, "grm_items.csv"), row.names = FALSE)
cat("\n== Graded response model: fit and reliability ==\n"); print(grm_fit, row.names = FALSE)

# Information curves
theta <- matrix(seq(-4, 4, by = 0.05))
info_test <- bind_rows(lapply(parts, function(pt) {
  ti <- testinfo(GRM[[pt]], theta)
  data.frame(part = pt, theta = theta[, 1], info = ti, se = 1 / sqrt(ti))
}))
info_item <- bind_rows(lapply(parts, function(pt) bind_rows(lapply(seq_len(ncol(CAT[[pt]])), function(k) {
  data.frame(part = pt, transition = k, theta = theta[, 1],
             info = iteminfo(extract.item(GRM[[pt]], k), theta))
}))))
info_item$target <- ifelse(info_item$part == "B" & to_letter(info_item$transition),
                           "to letter", "to number")

max_info <- max(info_test$info)               # SE is drawn on a rescaled second axis
info_test$se_scaled <- info_test$se * max_info / 2
g1 <- ggplot(info_test, aes(theta, info, colour = part)) +
  geom_line(linewidth = 1.1) +
  geom_line(aes(y = se_scaled), linetype = "dashed", linewidth = 0.6) +
  scale_colour_manual(values = part_cols, labels = c(A = "TMT-A", B = "TMT-B")) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 2 / max_info, name = "Standard error (dashed)")) +
  labs(x = "Processing speed θ (higher = faster)", y = "Test information (solid)", colour = NULL,
       title = "Where each part measures speed precisely")
save_png(g1, "irt_test_information.png")

g2 <- ggplot(info_item, aes(theta, info, group = transition, colour = transition)) +
  geom_line(aes(linetype = target), linewidth = 0.5, alpha = 0.85) +
  scale_colour_viridis_c(option = "C", end = 0.9) +
  scale_linetype_manual(values = c("to number" = "solid", "to letter" = "22")) +
  facet_wrap(~ part, labeller = labeller(part = c(A = "TMT-A", B = "TMT-B"))) +
  labs(x = "Processing speed θ", y = "Item information", colour = "Transition", linetype = NULL,
       title = "Item information for each of the 24 transitions")
save_png(g2, "irt_item_information.png", w = 9)

grm_par$target <- ifelse(grm_par$part == "B" & to_letter(grm_par$transition), "to letter", "to number")
g3 <- ggplot(grm_par, aes(transition, a, colour = part, shape = target)) +
  geom_line(aes(group = part), alpha = 0.4) +
  geom_point(size = 2.3) +
  scale_colour_manual(values = part_cols, labels = c(A = "TMT-A", B = "TMT-B")) +
  scale_shape_manual(values = c("to number" = 16, "to letter" = 17)) +
  labs(x = "Transition (1 = first circle to second)", y = "GRM discrimination a", colour = NULL, shape = NULL,
       title = "Which transitions best separate fast from slow participants")
save_png(g3, "irt_discrimination.png")

# ---- 4. Dimensionality: one speed trait or two? -------------------------------
# Compare a unidimensional GRM on all 48 items with a two-factor model where
# A items load on F_A, B items on F_B, and the factors correlate.
CAT_AB <- cbind(CAT$A, CAT$B)
uni <- mirt(CAT_AB, 1, itemtype = "graded", verbose = FALSE, technical = EM_CTRL)
two <- mirt(CAT_AB, mirt.model("FA = 1-24\nFB = 25-48\nCOV = FA*FB"),
            itemtype = "graded", verbose = FALSE, quadpts = 21, technical = EM_CTRL)
dim_cmp <- anova(uni, two)
r_AB <- coef(two, simplify = TRUE)$cov["FA", "FB"]
cat("\n== Dimensionality (1 vs 2 factors) ==\n"); print(round(dim_cmp, 1))
cat(sprintf("Latent correlation between Part A speed and Part B speed: %.3f\n", r_AB))
write.csv(cbind(model = c("one factor", "two factors"), round(dim_cmp, 2)),
          file.path(OUT, "dimensionality.csv"), row.names = FALSE)

# ---- 5. DIF by device ---------------------------------------------------------
# Same speed, different device: does any transition behave differently?
# Group means and variances are free, so true speed differences are not DIF.
dif_res <- NULL
if (RUN_DIF && DIF_VAR %in% names(p) && length(unique(p[[DIF_VAR]])) == 2) {
  grp <- as.character(p[[DIF_VAR]])
  pars <- c("a1", paste0("d", seq_len(N_CAT - 1)))
  dif_res <- bind_rows(lapply(parts, function(pt) {
    mg <- multipleGroup(CAT[[pt]], 1, group = grp, itemtype = "graded", verbose = FALSE,
                        invariance = c("slopes", "intercepts", "free_means", "free_var"))
    d  <- DIF(mg, which.par = pars, scheme = "drop", p.adjust = "BH", verbose = FALSE)
    data.frame(part = pt, item = rownames(d), round(d[, c("X2", "df", "p", "adj_p")], 4))
  }))
  write.csv(dif_res, file.path(OUT, "dif_by_device.csv"), row.names = FALSE)
  cat(sprintf("\n== DIF by %s: %d of %d transitions flagged (BH-adjusted p < .05) ==\n",
              DIF_VAR, sum(dif_res$adj_p < .05), nrow(dif_res)))
}

# ---- 6. Accuracy: Rasch model on errors --------------------------------------
# Item = "an error occurred on this transition". Part A errors are usually too
# rare to estimate; transitions with fewer than 15 errors are dropped.
acc <- lapply(parts, function(pt) {
  E <- (ERR[[pt]] > 0) * 1L
  keep <- colSums(E) >= 15
  if (sum(keep) < 5) return(NULL)
  mod <- mirt(E[, keep], 1, itemtype = "Rasch", verbose = FALSE)
  list(mod = mod, items = cbind(part = pt, item = colnames(E)[keep],
                                as.data.frame(round(coef(mod, IRTpars = TRUE, simplify = TRUE)$items, 3))))
})
names(acc) <- parts
acc_items <- bind_rows(lapply(acc, function(a) a$items))
if (nrow(acc_items)) write.csv(acc_items, file.path(OUT, "rasch_error_items.csv"), row.names = FALSE)
cat("\n== Rasch model on errors: transitions kept ==\n")
print(sapply(acc, function(a) if (is.null(a)) 0L else nrow(a$items)))

# ---- 7. Person scores and validity -------------------------------------------
sc <- data.frame(participant_id = p$participant_id,
                 theta_A = fscores(GRM$A, method = "EAP")[, 1],
                 theta_B = fscores(GRM$B, method = "EAP")[, 1],
                 speed_ln_A = LN$A$speed, speed_ln_B = LN$B$speed)
# B-specific speed: the part of B speed not explained by A speed (cf. B - A)
sc$theta_B_given_A <- resid(lm(theta_B ~ theta_A, data = sc))
if (!is.null(acc$B)) sc$error_proneness_B <- fscores(acc$B$mod, method = "EAP")[, 1]
scores <- merge(p, sc, by = "participant_id")
write.csv(sc, file.path(OUT, "person_scores.csv"), row.names = FALSE)

cat("\n== Validity: correlations ==\n")
print(round(with(scores, c(
  "theta_A vs log A_time" = cor(theta_A, log(A_time_s)),
  "theta_B vs log B_time" = cor(theta_B, log(B_time_s)),
  "theta_B vs lognormal B speed" = cor(theta_B, speed_ln_B),
  "theta_B|A vs B/A" = cor(theta_B_given_A, B_div_A)
)), 3))

auc <- function(score, case) {                 # P(score of a case > score of a control)
  r <- rank(score); n1 <- sum(case); n0 <- sum(!case)
  (sum(r[case]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
if ("group" %in% names(scores)) {
  cat("\n== Group means of IRT speed (higher = faster) ==\n")
  print(aggregate(cbind(theta_A, theta_B, theta_B_given_A) ~ group, data = scores,
                  FUN = function(x) round(mean(x), 2)))
  hc_mci <- scores[scores$group %in% c("HC", "MCI"), ]
  case <- hc_mci$group == "MCI"
  auc_tab <- data.frame(
    score = c("B_time_s", "-theta_B", "B_minus_A_s", "-theta_B_given_A", "B_div_A"),
    AUC_HC_vs_MCI = round(c(auc(hc_mci$B_time_s, case), auc(-hc_mci$theta_B, case),
                            auc(hc_mci$B_minus_A_s, case), auc(-hc_mci$theta_B_given_A, case),
                            auc(hc_mci$B_div_A, case)), 3))
  cat("\n== HC vs MCI discrimination (AUC) ==\n"); print(auc_tab, row.names = FALSE)
  write.csv(auc_tab, file.path(OUT, "auc_hc_vs_mci.csv"), row.names = FALSE)

  g4 <- ggplot(scores, aes(theta_B, B_time_s, colour = group)) +
    geom_point(alpha = 0.35, size = 0.9) +
    scale_y_log10() +
    scale_colour_manual(values = c(HC = "#0072B2", MCI = "#E69F00", Dementia = "#D55E00")) +
    labs(x = "TMT-B speed θ (GRM, higher = faster)", y = "TMT-B completion time (s, log scale)",
         colour = NULL, title = "IRT speed score against raw completion time")
  save_png(g4, "irt_theta_vs_time.png", w = 7)
}

# ---- 8. Explanatory IRT with lme4 --------------------------------------------
# All 96,000 transitions in one model. Item properties replace item parameters
# (an explanatory lognormal RT model, in the spirit of the LLTM):
#   log RT = part + distance + position + number->letter switch
#            + person speed (random intercept) + person B-specific slowness (random slope)
# Errors get the matching Poisson model (Rasch Poisson counts with item covariates).
if (requireNamespace("lme4", quietly = TRUE)) {
  d <- bind_rows(long$A, long$B)
  d$part_B    <- as.integer(d$part == "B")
  d$switch    <- as.integer(d$part == "B" & to_letter(d$transition_index))
  d$position  <- (d$transition_index - 12.5) / 23     # -0.5 (first) to 0.5 (last)
  d$dist_c    <- d$distance_norm - mean(d$distance_norm)
  rt_ok <- d[d$errors == 0, ]                         # error transitions carry correction time

  m_rt  <- lme4::lmer(log_rt ~ part_B * (dist_c + position) + switch + (1 + part_B | participant_id),
                      data = rt_ok, REML = FALSE)
  m_rt0 <- lme4::lmer(log_rt ~ part_B * (dist_c + position) + switch + (1 | participant_id),
                      data = rt_ok, REML = FALSE)
  fx <- summary(m_rt)$coefficients
  rt_fixed <- data.frame(term = rownames(fx), estimate = round(fx[, 1], 4), se = round(fx[, 2], 4),
                         pct_change = round(100 * (exp(fx[, 1]) - 1), 1))
  vc <- as.data.frame(lme4::VarCorr(m_rt))
  cat("\n== Explanatory lognormal model (lme4): fixed effects, % change in RT ==\n")
  print(rt_fixed, row.names = FALSE)
  cat(sprintf("Person SD: general speed %.3f, B-specific %.3f, correlation %.2f (log scale)\n",
              vc$sdcor[1], vc$sdcor[2], vc$sdcor[3]))
  cat("B-specific random slope vs none (LRT; p is conservative at the boundary):\n")
  print(anova(m_rt0, m_rt))
  write.csv(rt_fixed, file.path(OUT, "lmm_rt_fixed_effects.csv"), row.names = FALSE)

  m_err <- lme4::glmer(errors ~ part_B + switch + dist_c + position + (1 | participant_id),
                       data = d, family = poisson)
  fe <- summary(m_err)$coefficients
  err_fixed <- data.frame(term = rownames(fe), estimate = round(fe[, 1], 4), se = round(fe[, 2], 4),
                          rate_ratio = round(exp(fe[, 1]), 3))
  cat("\n== Explanatory Poisson error model (lme4): rate ratios ==\n")
  print(err_fixed, row.names = FALSE)
  cat(sprintf("Person SD of error proneness: %.3f (log rate)\n",
              as.data.frame(lme4::VarCorr(m_err))$sdcor[1]))
  write.csv(err_fixed, file.path(OUT, "glmm_error_fixed_effects.csv"), row.names = FALSE)
} else {
  cat("\nSection 8 skipped: install.packages(\"lme4\") to run the explanatory models.\n")
}

# ---- 9. Measurement invariance on continuous log RT (lavaan) -----------------
# The continuous counterpart of DIF (section 5 works on binned times instead).
# The one-factor lognormal model of section 2 is fitted per group, adding
# equality constraints step by step:
#   configural  same structure in every group
#   metric      equal loadings: transitions reflect speed equally strongly
#   scalar      equal intercepts: same speed -> same expected time per transition
#   strict      equal residual variances: same within-person noise
# With N = 2000 the chi-square difference rejects almost anything, so each step
# is judged by the change in fit: dCFI >= -.010 and dRMSEA <= .015
# (Cheung & Rensvold, 2002; Chen, 2007).
if (requireNamespace("lavaan", quietly = TRUE)) {
  mi_groups <- list()
  if ("device" %in% names(p)) mi_groups$device <- p$device
  if ("age" %in% names(p))
    mi_groups$age_band <- as.character(cut(p$age, c(-Inf, 64, 74, Inf), labels = c("<=64", "65-74", ">=75")))
  if ("education_years" %in% names(p))
    mi_groups$education <- ifelse(p$education_years <= 9, "<=9 years", ">9 years")  # 9 = compulsory schooling

  steps <- list(configural = character(0),
                metric     = "loadings",
                scalar     = c("loadings", "intercepts"),
                strict     = c("loadings", "intercepts", "residuals"))

  # Which equality constraints strain the model most (score tests), labelled by parameter
  flag_constraints <- function(fit, top = 5) {
    st   <- lavaan::lavTestScore(fit)$uni
    ptab <- lavaan::parTable(fit)
    lab  <- setNames(paste0(ptab$lhs, ptab$op, ptab$rhs), ptab$plabel)
    st$parameter <- lab[st$lhs]
    st <- st[order(-st$X2), ][seq_len(min(top, nrow(st))), ]
    data.frame(parameter = st$parameter, X2 = round(st$X2, 1), p = signif(st$p.value, 3))
  }

  mi_rows <- list(); mi_flags <- list()
  for (pt in parts) {
    items <- colnames(Y[[pt]])
    model <- paste("speed =~", paste(items, collapse = " + "))
    for (g in names(mi_groups)) {
      grp <- mi_groups[[g]]
      ok  <- !is.na(grp) & complete.cases(Y[[pt]])
      tab <- table(grp[ok])
      if (length(tab) < 2 || min(tab) < 100) next          # too few people per group
      dat  <- data.frame(Y[[pt]][ok, , drop = FALSE], grp = grp[ok])
      fits <- lapply(steps, function(eq) lavaan::cfa(model, data = dat, group = "grp", group.equal = eq))
      fm   <- t(sapply(fits, lavaan::fitMeasures, fit.measures = c("chisq", "df", "cfi", "rmsea")))
      d_cfi   <- c(NA, diff(fm[, "cfi"]))
      d_rmsea <- c(NA, diff(fm[, "rmsea"]))
      holds   <- c(NA, d_cfi[-1] >= -0.010 & d_rmsea[-1] <= 0.015)
      mi_rows[[paste(pt, g)]] <- data.frame(
        part = pt, grouping = g, groups = paste(names(tab), tab, sep = " n=", collapse = "; "),
        step = names(steps), round(fm, 3), d_cfi = round(d_cfi, 4), d_rmsea = round(d_rmsea, 4), holds = holds)
      if (isFALSE(holds[3]))
        mi_flags[[paste(pt, g)]] <- data.frame(part = pt, grouping = g, flag_constraints(fits$scalar))
    }
  }

  mi <- bind_rows(mi_rows)
  write.csv(mi, file.path(OUT, "measurement_invariance.csv"), row.names = FALSE)
  cat("\n== Measurement invariance (continuous log RT, lavaan) ==\n")
  print(mi[, c("part", "grouping", "step", "cfi", "rmsea", "d_cfi", "d_rmsea", "holds")], row.names = FALSE)
  if (length(mi_flags)) {
    flags <- bind_rows(mi_flags)
    write.csv(flags, file.path(OUT, "invariance_flagged_constraints.csv"), row.names = FALSE)
    cat("\nWhere scalar invariance fails, the constraints that strain the model most:\n")
    print(flags, row.names = FALSE)
  }
} else {
  cat("\nSection 9 skipped: install.packages(\"lavaan\") to test measurement invariance.\n")
}

cat("\nSaved results to ./", OUT, "\n", sep = "")
