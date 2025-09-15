####generate simulated data for nonlinear and linear system, which would be used for later comparsion

sim.func <- function(mite_lambda = exp(0.4488), n = 1e3, 
                     b0 = 0.9210, b1 = 0.1068, b2 = -0.002088, 
                     theta = 0.1) {
  
  set.seed(1)
  
  # Step 1: Simulate MITE load as the independent variable
  Mite_load <- rnbinom(n = n, mu = mite_lambda, size = theta)
  
  # Step 2: Make GREGARINE load a function of Mite_load
  log_lambda <- b0 + b1 * Mite_load + b2 * Mite_load^2
  Greg_load <- rnbinom(n = n, mu = exp(log_lambda), size = theta)
  
  return(data.frame(
    Mite_load = Mite_load,
    Greg_load = Greg_load,
    Mite_Infected = as.numeric(Mite_load > 0),
    Greg_Infected = as.numeric(Greg_load > 0)
  ))
}


## ---- No association scenario (match style of your LinCoi block) ----
sim.dat.NoCoi <- sim.func(b1 = 0, b2 = 0)

sim.means.NoCoi <- replicate(1000, {
  sim.dat.NoCoi %>%
    mutate(., 
           Mite_Infected_sim = sample(Mite_Infected, replace = FALSE),
           Greg_Infected_sim = sample(Greg_Infected, replace = FALSE)) %>%
    dplyr::select(., c("Mite_Infected_sim", "Greg_Infected_sim")) %>%
    C.score()
}) %>%
  as.data.frame() %>%
  summarise(
    mean  = mean(.),
    sd    = sd(.),
    n     = n(),
    upper = quantile(., prob = 0.975),
    lower = quantile(., prob = 1 - 0.975)
  ) %>%
  mutate(se = sd / sqrt(n), Species = "No association")

sim.means.NoCoi$Observed <- sim.dat.NoCoi %>%
  dplyr::select(., c("Mite_Infected", "Greg_Infected")) %>%
  C.score()

print(sim.means.NoCoi)

#######

sim.dat.LinCoi <- sim.func(b1=0.25, b2=0)

sim.means.LinCoi <- replicate(1000, {
  sim.dat.LinCoi %>%
    mutate(., Mite_Infected_sim = sample(Mite_Infected, replace=FALSE),
           Greg_Infected_sim = sample(Greg_Infected, replace=FALSE)) %>%
    dplyr::select(., c("Mite_Infected_sim","Greg_Infected_sim")) %>% 
    C.score() }) %>% 
  as.data.frame() %>%
  summarise(mean=mean(.), sd=sd(.), n=n(), upper=quantile(.,prob=0.975), lower=quantile(., prob=1-0.975)) %>% 
  mutate(se=sd/sqrt(n), Species="Linear association")

sim.means.LinCoi$Observed <- sim.dat.LinCoi %>% dplyr::select(., c("Mite_Infected","Greg_Infected")) %>% 
  C.score()

print(sim.means.LinCoi)


#######weak quadratic

sim.dat.WeakQuadCoi <- sim.func(b1=0.25, b2=-0.0005)

sim.means.WeakQuadCoi <- replicate(1000, {
  sim.dat.WeakQuadCoi %>%
    mutate(., Mite_Infected_sim = sample(Mite_Infected, replace=FALSE),
           Greg_Infected_sim = sample(Greg_Infected, replace=FALSE)) %>%
    dplyr::select(., c("Mite_Infected_sim","Greg_Infected_sim")) %>% 
    C.score() }) %>% 
  as.data.frame() %>%
  summarise(mean=mean(.), sd=sd(.), n=n(), upper=quantile(.,prob=0.975), lower=quantile(., prob=1-0.975)) %>% 
  mutate(se=sd/sqrt(n), Species="Weak quadratic association")

sim.means.WeakQuadCoi$Observed <- sim.dat.WeakQuadCoi %>% dplyr::select(., c("Mite_Infected","Greg_Infected")) %>% 
  C.score()

print(sim.means.WeakQuadCoi)

#######strong quadratic

sim.dat.QuadCoi <- sim.func(b1=0.25, b2=-0.1)

attach(sim.dat.QuadCoi)

plot(Mite_load, Greg_load, type = "o")

sim.means.QuadCoi <- replicate(1000, {
  sim.dat.QuadCoi %>%
    mutate(., Mite_Infected_sim = sample(Mite_Infected, replace=FALSE),
           Greg_Infected_sim = sample(Greg_Infected, replace=FALSE)) %>%
    dplyr::select(., c("Mite_Infected_sim","Greg_Infected_sim")) %>% 
    C.score() }) %>% 
  as.data.frame() %>%
  summarise(mean=mean(.), sd=sd(.), n=n(), upper=quantile(.,prob=0.975), lower=quantile(., prob=1-0.975)) %>% 
  mutate(se=sd/sqrt(n), Species="Strong quadratic association")

sim.means.QuadCoi$Observed <- sim.dat.QuadCoi %>% dplyr::select(., c("Mite_Infected","Greg_Infected")) %>% 
  C.score()

print(sim.means.QuadCoi)

sim.dat <- rbind(sim.means.NoCoi, sim.means.LinCoi, sim.means.WeakQuadCoi, sim.means.QuadCoi)

sim.dat$Species <- factor(sim.dat$Species, levels=c("No association", "Linear association", "Weak quadratic association", "Strong quadratic association"))



suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
})

# ------- Core C-score helpers (pair: 2 columns = two parasites) -------
cscore_pair <- function(pa) {
  stopifnot(is.matrix(pa), ncol(pa) == 2)
  A <- pa[,1] == 1; B <- pa[,2] == 1
  b <- sum(A & !B)   # only A
  c <- sum(!A & B)   # only B
  b * c              # checkerboard units
}

cscore_in_bin <- function(df, spA = "Mite_Infected", spB = "Greg_Infected") {
  pa <- as.matrix(df[, c(spA, spB)])
  tibble(n = nrow(df),
         cscore = cscore_pair(pa),
         prevA = mean(pa[,1]),
         prevB = mean(pa[,2]))
}

# Label-permutation test (two-sided) for difference between datasets within a bin
permtest_diff <- function(dfA, dfB, spA = "Mite_Infected", spB = "Greg_Infected",
                          n_perm = 999, seed = 1) {
  set.seed(seed)
  nA <- nrow(dfA); nB <- nrow(dfB)
  dat <- rbind(
    data.frame(dfA[, c(spA, spB)]),
    data.frame(dfB[, c(spA, spB)])
  )
  obsA <- cscore_pair(as.matrix(dfA[, c(spA, spB)]))
  obsB <- cscore_pair(as.matrix(dfB[, c(spA, spB)]))
  obs_diff <- obsA - obsB
  
  null <- replicate(n_perm, {
    idxA <- sample.int(nA + nB, nA, replace = FALSE)
    A2 <- as.matrix(dat[idxA, ])
    B2 <- as.matrix(dat[-idxA, ])
    cscore_pair(A2) - cscore_pair(B2)
  })
  
  p_two <- (sum(abs(null) >= abs(obs_diff)) + 1) / (n_perm + 1)
  z_diff <- (obs_diff - mean(null)) / sd(null)
  list(diff = obs_diff, p = p_two, z = z_diff, nA = nA, nB = nB)
}

# ------- Validity checks for a candidate bin -------
has_variation <- function(df, spA = "Mite_Infected", spB = "Greg_Infected") {
  length(unique(df[[spA]])) > 1 && length(unique(df[[spB]])) > 1
}

bin_ok <- function(dfA, dfB, min_n_per_set = 5, spA = "Mite_Infected", spB = "Greg_Infected") {
  nrow(dfA) >= min_n_per_set && nrow(dfB) >= min_n_per_set &&
    has_variation(dfA, spA, spB) && has_variation(dfB, spA, spB)
}

# ------- Evaluate one bin [L, U] -------
evaluate_bin <- function(L, U, dat_lin, dat_nonlin,
                         x = "Mite_load", spA = "Mite_Infected", spB = "Greg_Infected",
                         n_perm = 999, min_n_per_set = 5, seed = 1) {
  A <- dat_lin    %>% filter(.data[[x]] >= L, .data[[x]] < U)
  B <- dat_nonlin %>% filter(.data[[x]] >= L, .data[[x]] < U)
  if (!bin_ok(A, B, min_n_per_set, spA, spB)) {
    return(tibble(L = L, U = U, mid = (L+U)/2,
                  nA = nrow(A), nB = nrow(B),
                  cA = NA_real_, cB = NA_real_,
                  diff = NA_real_, p = NA_real_, z = NA_real_,
                  reason = "too_small_or_no_variation"))
  }
  csA <- cscore_in_bin(A, spA, spB)
  csB <- cscore_in_bin(B, spA, spB)
  pt  <- permtest_diff(A, B, spA, spB, n_perm = n_perm, seed = seed)
  tibble(L = L, U = U, mid = (L+U)/2,
         nA = pt$nA, nB = pt$nB,
         cA = csA$cscore, cB = csB$cscore,
         diff = pt$diff, p = pt$p, z = pt$z,
         reason = NA_character_)
}

# ------- Step 1: Random partitions to find ONE promising bin -------
# Make R random partitions of K bins each over the pooled range; score all bins; pick best by p, then |diff|
random_bin_search <- function(dat_lin, dat_nonlin,
                              x = "Mite_load",
                              K = 3,     # bins per partition
                              R = 100,   # number of random partitions
                              n_perm = 999,
                              min_n_per_set = 5,
                              seed = 42) {
  set.seed(seed)
  x_all <- c(dat_lin[[x]], dat_nonlin[[x]])
  xmin <- min(x_all, na.rm = TRUE)
  xmax <- max(x_all, na.rm = TRUE)
  
  results <- vector("list", R)
  for (r in seq_len(R)) {
    # sample K-1 interior cutpoints uniformly, sort, add endpoints
    cuts <- sort(runif(K - 1, min = xmin, max = xmax))
    br   <- c(xmin, cuts, xmax + 1e-9)  # open upper bound
    bins <- map_dfr(seq_len(K), function(i) {
      evaluate_bin(br[i], br[i+1], dat_lin, dat_nonlin,
                   x = x, n_perm = n_perm, min_n_per_set = min_n_per_set, seed = seed + r + i)
    }) %>% mutate(partition = r, bin_id = row_number())
    results[[r]] <- bins
  }
  all_bins <- bind_rows(results)
  
  # pick best bin: min p (non-NA), tiebreak by largest |diff|
  cand <- all_bins %>% filter(!is.na(p))
  if (nrow(cand) == 0) {
    return(list(all_bins = all_bins, best = NULL))
  }
  ord <- cand %>% arrange(p, desc(abs(diff)))
  list(all_bins = all_bins, best = ord[1, ])
}

# ------- Step 2: Local shrink around the chosen bin -------
# Given best [L,U], try shrinking/expanding by small steps while keeping >=5 per set
local_shrink_search <- function(best_row, dat_lin, dat_nonlin,
                                x = "Mite_load",
                                step = NULL,         # step size in original units; if NULL, use 1/20 of width
                                n_steps = 6,         # how many steps inward/outward to try on each side
                                n_perm = 999,
                                min_n_per_set = 5,
                                seed = 99) {
  stopifnot(!is.null(best_row), nrow(best_row) == 1)
  L0 <- best_row$L; U0 <- best_row$U
  w0 <- U0 - L0
  if (is.null(step)) step <- w0 / 20  # gentle granularity
  
  # grid of (dL, dU) offsets to shrink/expand, but enforce L' < U'
  dLs <- seq(-n_steps, n_steps) * step   # negative = expand left, positive = shrink rightwards
  dUs <- seq(-n_steps, n_steps) * step   # negative = shrink leftwards (move U left), positive = expand right
  grid <- expand.grid(dL = dLs, dU = dUs)
  grid <- grid %>% filter(L0 + dL < U0 + dU)
  
  trials <- pmap_dfr(grid, function(dL, dU) {
    L <- L0 + dL
    U <- U0 + dU
    evaluate_bin(L, U, dat_lin, dat_nonlin,
                 x = x, n_perm = n_perm, min_n_per_set = min_n_per_set, seed = seed)
  }) %>% arrange(p, desc(abs(diff)))
  
  # best by p, then |diff|
  best2 <- trials %>% filter(!is.na(p)) %>% slice_min(order_by = p, n = 1, with_ties = FALSE)
  list(trials = trials, best = if (nrow(best2)) best2 else NULL)
}

# ============================
# ====== RUN THE PIPELINE =====
# ============================

# Your two datasets
dat_lin    <- sim.dat.LinCoi
dat_nonlin <- sim.dat.QuadCoi

# 1) Random partitions to find one promising bin
rb <- random_bin_search(
  dat_lin, dat_nonlin,
  x = "Mite_load",
  K = 3,             # random 3-bin partitions (adjustable)
  R = 200,           # number of random partitions to try
  n_perm = 999,
  min_n_per_set = 5, # as requested
  seed = 2025
)

cat("\n=== Random bin search: best candidate ===\n")
print(rb$best)

# 2) Local shrink around that bin (skip if none found)
if (!is.null(rb$best)) {
  ls <- local_shrink_search(
    best_row = rb$best,
    dat_lin, dat_nonlin,
    x = "Mite_load",
    step = NULL,       # defaults to width/20
    n_steps = 8,       # explores +/- 8 steps on each side
    n_perm = 999,
    min_n_per_set = 5,
    seed = 303
  )
  
  cat("\n=== Local shrink: best refined bin ===\n")
  print(ls$best)
  
  # ---- Plot the refinement landscape with exact p-values ----
  trials <- ls$trials %>%
    mutate(width = U - L,
           label = ifelse(is.na(p), "NA", sprintf("p=%.3f", p)))
  
  p_land <- ggplot(trials, aes(x = (L+U)/2, y = width)) +
    geom_point(aes(fill = p, size = nA + nB), shape = 21, color = "black", na.rm = TRUE) +
    geom_text(aes(label = label), vjust = -0.8, size = 3, na.rm = TRUE) +
    scale_fill_viridis_c(option = "C", direction = -1, na.value = "grey85") +
    labs(title = "Local shrink/expand around the promising bin",
         subtitle = "Point fill = exact p-value; size = total n per bin",
         x = "Bin midpoint (Mite_load)",
         y = "Bin width") +
    theme_minimal() +
    guides(fill = guide_colorbar(title = "p-value"),
           size = guide_legend(title = "Total n"))
  
  print(p_land)
  
  # ---- Simple “best vs original” comparison with exact p ----
  if (nrow(ls$best)) {
    best_lab <- sprintf("Refined bin: [%.2f, %.2f], p = %.3f, ΔC = %.0f",
                        ls$best$L, ls$best$U, ls$best$p, ls$best$diff)
    orig_lab <- sprintf("Original bin: [%.2f, %.2f], p = %.3f, ΔC = %.0f",
                        rb$best$L, rb$best$U, rb$best$p, rb$best$diff)
    message(best_lab, "\n", orig_lab)
  }
} else {
  message("No usable bin found in random search (all failed min_n/variation). Consider increasing R, lowering min_n_per_set, or broadening K.")
}
