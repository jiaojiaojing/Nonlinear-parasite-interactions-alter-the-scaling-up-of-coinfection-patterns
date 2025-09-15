# ============================
# FAST sweep: theta -> min p
# ============================

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(ggplot2)
  # Optional parallel:
  # library(future); library(future.apply)
})

# ---- knobs ----
theta_grid <- c(0.001, 0.01, 0.02, 0.04, 0.1)  # larger theta = less noise
b1_common  <- 0.25
b2_linear  <- 0.0
b2_quadr   <- -0.20
n_hosts    <- 1200

# search sizes (FAST)
K_bins_coarse  <- 3
R_parts_coarse <- 60       # was 200
n_steps_coarse <- 4        # was 8
n_perm_coarse  <- 199      # was 999

# refine only the best theta
n_perm_refine  <- 1999

min_n_set      <- 5

# ---- optional parallel ----
# plan(multisession, workers = max(1, parallel::detectCores() - 1))

run_one_theta_fast <- function(th, seed_base = 1L) {
  set.seed(seed_base + as.integer(th*1e6) %% 1e6)
  
  # Simulate independent datasets
  seed_lin <- sample.int(1e7, 1)
  seed_nl  <- sample.int(1e7, 1)
  set.seed(seed_lin)
  dat_lin  <- sim.func(n = n_hosts, b1 = b1_common, b2 = b2_linear, theta = th)
  set.seed(seed_nl)
  dat_quad <- sim.func(n = n_hosts, b1 = b1_common, b2 = b2_quadr,  theta = th)
  
  # Coarse random-bin search
  rb <- random_bin_search(
    dat_lin, dat_quad,
    x = "Mite_load",
    K = K_bins_coarse, R = R_parts_coarse,
    n_perm = n_perm_coarse,
    min_n_per_set = min_n_set,
    seed = sample.int(1e7, 1)
  )
  
  if (is.null(rb$best)) {
    return(list(theta = th, min_p_coarse = NA_real_, best_row = NULL, trials = NULL))
  }
  
  # Early stop: if the best random-bin p is very large, skip shrinking
  if (is.na(rb$best$p) || rb$best$p > 0.5) {
    return(list(theta = th, min_p_coarse = rb$best$p, best_row = rb$best, trials = NULL))
  }
  
  # Local shrink with coarse permutations
  ls <- local_shrink_search(
    best_row = rb$best,
    dat_lin  = dat_lin, dat_nonlin = dat_quad,
    x = "Mite_load",
    step = NULL,        # width/20
    n_steps = n_steps_coarse,
    n_perm  = n_perm_coarse,
    min_n_per_set = min_n_set,
    seed = sample.int(1e7, 1)
  )
  
  trials <- ls$trials
  min_p  <- suppressWarnings(min(trials$p, na.rm = TRUE))
  
  list(theta = th, min_p_coarse = min_p, best_row = rb$best, trials = trials)
}

# ---- run coarse sweep (parallel-ready) ----
# coarse_res <- future_lapply(theta_grid, run_one_theta_fast)   # if using future.apply
coarse_res <- lapply(theta_grid, run_one_theta_fast)

coarse_tbl <- bind_rows(lapply(coarse_res, function(x) {
  tibble(theta = x$theta, min_p_coarse = x$min_p_coarse)
}))
print(coarse_tbl)

# ---- pick the best theta to refine ----
best_idx <- which.min(coarse_tbl$min_p_coarse)
best_theta <- theta_grid[best_idx]
message(sprintf("Refining theta = %s", best_theta))

# ---- refine the best theta with larger n_perm only ----
ref <- run_one_theta_fast(best_theta)  # rerun to get datasets & best row
if (!is.null(ref$best_row)) {
  # Re-simulate to avoid reusing the same seeds as coarse run
  set.seed(12345 + best_idx)
  dat_lin_ref  <- sim.func(n = n_hosts, b1 = b1_common, b2 = b2_linear, theta = best_theta)
  dat_quad_ref <- sim.func(n = n_hosts, b1 = b1_common, b2 = b2_quadr,  theta = best_theta)
  
  ls_ref <- local_shrink_search(
    best_row = ref$best_row,
    dat_lin  = dat_lin_ref, dat_nonlin = dat_quad_ref,
    x = "Mite_load",
    step = NULL,
    n_steps = n_steps_coarse,
    n_perm  = n_perm_refine,       # high precision only here
    min_n_per_set = min_n_set,
    seed = 8675309
  )
  
  trials_ref <- ls_ref$trials
  min_p_ref  <- suppressWarnings(min(trials_ref$p, na.rm = TRUE))
} else {
  trials_ref <- NULL
  min_p_ref  <- NA_real_
}

# ---- final table for plotting ----
final_tbl <- coarse_tbl %>%
  mutate(min_p_final = ifelse(theta == best_theta, min_p_ref, min_p_coarse))

print(final_tbl)

# ---- plot: min p vs theta ----
ggplot(final_tbl, aes(theta, min_p_final)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = theta_grid) +
  labs(
    title    = "Minimum p-value vs. noise (theta) — fast, two-stage search",
    subtitle = "Coarse n_perm for all θ; refined n_perm only for best θ",
    x = expression(theta~"(NB size; larger = less noise)"),
    y = "Minimum p-value across bins"
  ) +
  theme_minimal(base_size = 12)
