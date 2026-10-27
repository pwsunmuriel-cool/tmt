# =============================================================================
# Trail Making Test simulated data generator
#   Output columns match the CSVs exported by tmt_mobile_en.html exactly
#   (plus one extra `group` column), so you can build the analysis pipeline on
#   simulated data first and swap in real data later.
#
#   tmt_sim_participants.csv  one row per participant (summary metrics)
#   tmt_sim_transitions.csv   one row per transition (long format, inter-node RT)
#
# Calibration: age/education trends in completion time roughly follow the
#   Tombaugh (2004) paper-and-pencil norms. Touch tapping is usually faster than
#   drawing; adjust TAP_SPEED_FACTOR to change overall speed.
# Base R only; no packages required.
# =============================================================================

set.seed(20261024)

N                <- 2000   # number of participants
N_NODES          <- 25     # number of circles (same as the app's standard version)
ASPECT_H         <- 4 / 3  # portrait phone board: width = 1, height = 4/3 (same as the app)
TAP_SPEED_FACTOR <- 1.00   # < 1 means faster than paper norms
OUT_DIR          <- "."

# ---- 1. Demographics ---------------------------------------------------------
age <- round(runif(N, 55, 89))
sex <- sample(c("F", "M"), N, replace = TRUE, prob = c(0.53, 0.47))

# Cohort effect: older adults have fewer years of education on average
education_years <- round(rnorm(N, mean = 16 - 0.10 * (age - 20), sd = 2.8))
education_years <- pmin(pmax(education_years, 0), 22)

# Cognitive status: HC / MCI / AD; prevalence rises with age, education is protective
a     <- age - 20
age_mid <- c(57, 62, 67, 72, 77, 82, 87)   # Sun et al. (2014) Table 2; 55-64 extrapolated
p_mci <- approx(age_mid, c(.080, .110, .139, .164, .220, .235, .262), xout = age, rule = 2)$y
p_ad  <- approx(age_mid, c(.010, .020, .034, .035, .072, .130, .219), xout = age, rule = 2)$y
edu_or <- exp(-0.08 * (education_years - (16 - 0.10 * a)))   # education vs same-age peers
u     <- runif(N)
group <- ifelse(u < p_ad * edu_or, "AD", ifelse(u < (p_ad + p_mci) * edu_or, "MCI", "HC"))
mci   <- as.numeric(group == "MCI")
ad    <- as.numeric(group == "AD")

device       <- sample(c("iOS", "Android"), N, replace = TRUE, prob = c(0.55, 0.45))
session_time <- format(as.POSIXct("2026-01-05 09:00:00", tz = "UTC") +
                         runif(N, 0, 290 * 86400), "%Y-%m-%dT%H:%M:%SZ")
participant_id <- sprintf("SIM%04d", seq_len(N))

# ---- 2. Completion time (seconds, before error penalties) --------------------
# Shared processing-speed latent factor z1 correlates A and B (residual r = .30)
z1 <- rnorm(N); z2 <- rnorm(N)
logA <- log(21.8) + 0.004 * a + 0.00015 * a^2 - 0.020 * (education_years - 12) +
        0.25 * mci + 0.60 * ad + 0.25 * z1
logB <- log(49) - 0.0022 * a + 0.00027 * a^2 - 0.030 * (education_years - 12) +
        0.35 * mci + 0.75 * ad + 0.32 * (0.30 * z1 + sqrt(1 - 0.30^2) * z2)
base_s <- list(A = exp(logA) * TAP_SPEED_FACTOR, B = exp(logB) * TAP_SPEED_FACTOR)

# Error counts (Poisson); far more in B than A
lambda <- list(
  A = exp(-3.1 + 0.020 * a + 0.1 * mci + 0.1 * ad),
  B = exp(-1.8 + 0.020 * a + 0.70 * mci + 0.80 * ad - 0.04 * (education_years - 12))
)

# Intra-individual variability (IIV) and attentional lapses: larger with age and impairment
iiv     <- pmax(0.12, 0.22 + 0.003 * a + 0.08 * mci + 0.18 * ad + rnorm(N, 0, 0.05))
p_lapse <- 0.02 + 0.0008 * a + 0.04 * mci + 0.10 * ad

# ---- 3. Layout (same algorithm as the app: random placement + minimum spacing) --
make_layout <- function(n, H = ASPECT_H) {
  rU <- 0.25 * sqrt(H / n)               # app: r = 0.25 * sqrt(w*h/n), in board-width units
  pad <- 1.35 * rU; min_d <- 2.6 * rU
  repeat {
    u <- numeric(0); v <- numeric(0); tries <- 0
    while (length(u) < n && tries < 4000) {
      tries <- tries + 1
      x <- pad + runif(1) * (1 - 2 * pad)
      y <- pad + runif(1) * (H - 2 * pad)
      if (length(u) == 0 || all((u - x)^2 + (v - y)^2 >= min_d^2)) { u <- c(u, x); v <- c(v, y) }
    }
    if (length(u) == n) return(cbind(u, v))
    min_d <- min_d * 0.98
  }
}

make_sequence <- function(part, n) {
  if (part == "A") return(as.character(seq_len(n)))
  i <- seq_len(n) - 1
  ifelse(i %% 2 == 0, as.character(i / 2 + 1), LETTERS[(i - 1) / 2 + 1])
}
seqs <- list(A = make_sequence("A", N_NODES), B = make_sequence("B", N_NODES))

# ---- 4. Transition-level simulation ------------------------------------------
n_tr  <- N_NODES - 1
k     <- seq_len(n_tr)
diag  <- sqrt(1 + ASPECT_H^2)
trans_list <- vector("list", 2 * N)
summ <- data.frame(participant_id)

for (part in c("A", "B")) {
  cols <- matrix(NA_real_, N, 7,
                 dimnames = list(NULL, paste0(part, c("_time_s", "_errors", "_init_ms", "_mean_rt_ms",
                                                      "_rt_sd_ms", "_rt_cv", "_interrupted"))))
  seq_p <- seqs[[part]]
  for (i in seq_len(N)) {
    lay  <- make_layout(N_NODES)
    dist <- sqrt(rowSums((lay[-1, , drop = FALSE] - lay[-N_NODES, , drop = FALSE])^2)) / diag

    init_ms <- rlnorm(1, log(850 + 8 * a[i] + 300 * mci[i] + 700 * ad[i]), 0.30)

    # Weights: distance x remaining distractors (later = easier to find) x within-person noise x lapses
    w <- (0.55 + 0.45 * dist / mean(dist)) *
         (1 + (part == "A") * 0.35 * ((n_tr - k) / n_tr - 0.5)) *
         rlnorm(n_tr, 0, iiv[i])
    if (part == "B") {
      to_letter <- k %% 2 == 1                     # odd transitions are number -> letter
      w <- w * ifelse(to_letter, 1 + 0.012 * k, 1) # later letters are harder to retrieve
    }
    lapse <- runif(n_tr) < p_lapse[i]
    w[lapse] <- w[lapse] * runif(sum(lapse), 1.8, 3.5)

    rt <- w / sum(w) * max(base_s[[part]][i] * 1000 - init_ms, 400 * n_tr)

    # Errors: each adds ~0.6-3 s of correction time to its transition; B errors skew late
    n_err <- min(rpois(1, lambda[[part]][i]), 12)
    err_k <- integer(n_tr)
    if (n_err > 0) {
      where <- sample(k, n_err, replace = TRUE, prob = if (part == "B") k else rep(1, n_tr))
      for (j in where) {
        err_k[j] <- err_k[j] + 1
        rt[j] <- rt[j] + rlnorm(1, log(1400), 0.45)
      }
    }
    rt <- round(rt, 1)

    trans_list[[(part == "B") * N + i]] <- data.frame(
      participant_id = participant_id[i], part = part, n_nodes = N_NODES,
      transition_index = k, from_node = seq_p[-N_NODES], to_node = seq_p[-1],
      distance_norm = round(dist, 4), rt_ms = rt, errors = err_k
    )
    m <- mean(rt); s <- sd(rt)
    cols[i, ] <- c((init_ms + sum(rt)) / 1000, n_err, init_ms, m, s, s / m,
                   as.numeric(runif(1) < 0.01))   # ~1% left the app mid-test
  }
  summ <- cbind(summ, cols)
}

participants <- data.frame(
  participant_id, age, sex, education_years, group, device, session_time, n_nodes = N_NODES,
  summ[, -1]
)
participants$B_minus_A_s <- participants$B_time_s - participants$A_time_s
participants$B_div_A     <- participants$B_time_s / participants$A_time_s
num <- vapply(participants, is.double, logical(1))
participants[num] <- lapply(participants[num], round, 3)

transitions <- do.call(rbind, trans_list)

write.csv(participants, file.path(OUT_DIR, "tmt_sim_participants.csv"), row.names = FALSE)
write.csv(transitions,  file.path(OUT_DIR, "tmt_sim_transitions.csv"),  row.names = FALSE)

# ---- 5. Quick checks ---------------------------------------------------------
cat("Participants:", nrow(participants), " Transition rows:", nrow(transitions), "\n\n")
print(table(group))
cat("Part A > 150 s:", sum(participants$A_time_s > 150), "  Part B > 300 s:", sum(participants$B_time_s > 300), "
")
participants$age_band <- cut(participants$age, c(54, 64, 74, 89))
cat("\nMean completion time / errors by age band (HC only):\n")
hc <- participants[participants$group == "HC", ]
print(aggregate(cbind(A_time_s, B_time_s, B_div_A, A_errors, B_errors, B_rt_cv) ~ age_band,
                data = hc, FUN = function(x) round(mean(x), 2)))
cat("\nBy cognitive status:\n")
print(aggregate(cbind(A_time_s, B_time_s, B_div_A, B_errors, B_rt_cv) ~ group,
                data = participants, FUN = function(x) round(mean(x), 2)))
cat("\ncor(log A, log B) =", round(cor(log(participants$A_time_s), log(participants$B_time_s)), 2), "\n")
