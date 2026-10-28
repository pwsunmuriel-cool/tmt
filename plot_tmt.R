# =============================================================================
# Trail Making Test: suggested ggplot2 figures
#   Reads tmt_sim_participants.csv / tmt_sim_transitions.csv (or real exports
#   from the app, same columns) and writes PNGs to ./figures
#   Packages: ggplot2, dplyr
# =============================================================================

library(ggplot2)
library(dplyr)

dir.create("figures", showWarnings = FALSE)

p  <- read.csv("tmt_sim_participants.csv")
tr <- read.csv("tmt_sim_transitions.csv")

# Real app data has no `group` column; treat everyone as one group
if (!"group" %in% names(p)) p$group <- "All"
p$group <- factor(p$group, levels = intersect(c("HC", "MCI", "Dementia", "All"), unique(p$group)))
tr <- left_join(tr, p[, c("participant_id", "age", "group")], by = "participant_id")

# Colorblind-safe (Okabe-Ito)
group_cols <- c(HC = "#0072B2", MCI = "#E69F00", Dementia = "#D55E00", All = "#555555")
part_cols  <- c(A = "#56B4E9", B = "#CC79A7")

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(), legend.position = "top"))

save <- function(plot, name, w = 7, h = 4.5) {
  ggsave(file.path("figures", name), plot, width = w, height = h, dpi = 200, bg = "white")
}

# Participant-level data in long format (one row per participant x part)
p_long <- bind_rows(
  transmute(p, participant_id, age, group, part = "A", time_s = A_time_s, errors = A_errors, rt_cv = A_rt_cv),
  transmute(p, participant_id, age, group, part = "B", time_s = B_time_s, errors = B_errors, rt_cv = B_rt_cv)
)

# ---- Fig 1. Age trend: completion time vs age --------------------------------
f1 <- ggplot(p_long, aes(age, time_s, colour = group)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 1) +
  scale_y_log10() +
  scale_colour_manual(values = group_cols) +
  facet_wrap(~ part, labeller = labeller(part = c(A = "TMT-A", B = "TMT-B"))) +
  labs(x = "Age (years)", y = "Completion time (s, log scale)", colour = NULL,
       title = "Completion time rises with age, more steeply in TMT-B")
save(f1, "fig1_age_trend.png", w = 9)

# ---- Fig 2. Group differences: A, B and B/A ----------------------------------
p_idx <- bind_rows(
  transmute(p, group, index = "TMT-A time (s)", value = A_time_s),
  transmute(p, group, index = "TMT-B time (s)", value = B_time_s),
  transmute(p, group, index = "B / A ratio",     value = B_div_A)
)
p_idx$index <- factor(p_idx$index, levels = c("TMT-A time (s)", "TMT-B time (s)", "B / A ratio"))
f2 <- ggplot(p_idx, aes(group, value, fill = group)) +
  geom_violin(alpha = 0.35, colour = NA, scale = "width") +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
  scale_fill_manual(values = group_cols, guide = "none") +
  facet_wrap(~ index, scales = "free_y") +
  labs(x = NULL, y = NULL, title = "Group differences in TMT indices")
save(f2, "fig2_group_indices.png", w = 9)

# ---- Fig 3. Transition profile: RT across the 24 transitions -----------------
# RTs are right-skewed, so average on the log scale (geometric mean)
prof <- tr |>
  group_by(part, group, transition_index) |>
  summarise(m = mean(log(rt_ms)), se = sd(log(rt_ms)) / sqrt(n()), .groups = "drop")
f3 <- ggplot(prof, aes(transition_index, exp(m), colour = group, fill = group)) +
  geom_ribbon(aes(ymin = exp(m - 1.96 * se), ymax = exp(m + 1.96 * se)), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_y_log10() +
  scale_colour_manual(values = group_cols) +
  scale_fill_manual(values = group_cols) +
  facet_wrap(~ part, labeller = labeller(part = c(A = "TMT-A", B = "TMT-B"))) +
  labs(x = "Transition (1 = first circle to second)", y = "Geometric mean RT (ms, log scale, 95% CI)",
       colour = NULL, fill = NULL, title = "Inter-node RT profile across the trail")
save(f3, "fig3_transition_profile.png", w = 9)

# ---- Fig 4. Distance effect: RT vs distance between circles ------------------
dist_df <- tr |>
  filter(errors == 0) |>                         # error transitions carry correction time
  mutate(dist_bin = cut(distance_norm, quantile(distance_norm, 0:10 / 10), include.lowest = TRUE)) |>
  group_by(part, group, dist_bin) |>
  summarise(dist = median(distance_norm), rt = median(rt_ms), .groups = "drop")
f4 <- ggplot(dist_df, aes(dist, rt, colour = group, linetype = part)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = group_cols) +
  scale_linetype_discrete(labels = c(A = "TMT-A", B = "TMT-B")) +
  labs(x = "Distance between circles (proportion of board diagonal, decile medians)",
       y = "Median inter-node RT (ms)", colour = NULL, linetype = NULL,
       title = "Longer jumps take longer (error-free transitions)")
save(f4, "fig4_distance_effect.png")

# ---- Fig 5. Intra-individual variability (RT CV) -----------------------------
f5 <- ggplot(p_long, aes(age, rt_cv, colour = group)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 1) +
  scale_colour_manual(values = group_cols) +
  facet_wrap(~ part, labeller = labeller(part = c(A = "TMT-A", B = "TMT-B"))) +
  labs(x = "Age (years)", y = "RT coefficient of variation (SD / mean)", colour = NULL,
       title = "Response-time variability within a person")
save(f5, "fig5_rt_variability.png", w = 9)

# ---- Fig 6. Error distribution ------------------------------------------------
err_df <- p_long |>
  mutate(err_cat = factor(pmin(errors, 4), levels = 0:4, labels = c("0", "1", "2", "3", "4+"))) |>
  count(part, group, err_cat) |>
  group_by(part, group) |>
  mutate(prop = n / sum(n))
f6 <- ggplot(err_df, aes(err_cat, prop, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  scale_fill_manual(values = group_cols) +
  facet_wrap(~ part, labeller = labeller(part = c(A = "TMT-A", B = "TMT-B"))) +
  labs(x = "Number of errors", y = "Share of participants", fill = NULL,
       title = "Errors are rare in TMT-A and common in impaired groups on TMT-B")
save(f6, "fig6_errors.png", w = 9)

cat("Saved", length(list.files("figures", pattern = "png$")), "figures to ./figures\n")
