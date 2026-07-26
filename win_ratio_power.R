# =============================================================================
#  Hierarchical Win Ratio Simulation
#  Matches the interactive HTML calculator exactly.
#
#  Endpoint hierarchy:
#    Level 1 — All-Cause Mortality (ACM)   : Weibull, shape = nu
#    Level 2 — CV Hospitalisation (CVH)    : NB first-event (Lomax), kappa
#  ACM & CVH correlated via Gaussian copula (rho)
#  LTFU: uniform dropout (total % over study period), independent per arm
#  Randomisation ratio: nT / nC = ratio
#
#  Variance of CVH count = mean + mean^2 * kappa  (NB parameterisation)
# =============================================================================

# ---- USER-ADJUSTABLE ASSUMPTIONS (mirrors HTML sliders) --------------------
fu          <- 30       # Follow-up duration (months)
nC          <- 300      # Control arm size
ratio       <- 1        # Randomisation ratio trt:ctrl  (e.g. 2 → 2:1)
ltfuc_pct   <- 5        # % lost to follow-up, control  (over whole study)
ltfut_pct   <- 5        # % lost to follow-up, treatment
acm_rate    <- 0.25     # Placebo ACM event rate (proportion by end of fu)
cvh_rate    <- 0.30     # Placebo CVH rate (events per year)
hr          <- 0.77     # HR for ACM  (treatment vs control)
rr          <- 0.50     # RR for CVH  (treatment vs control)
nu          <- 1.50     # Weibull shape for ACM
kappa       <- 1.50     # NB overdispersion for CVH  (Var = mean + mean^2*kappa)
rho         <- 0.50     # Gaussian copula correlation (ACM-CVH)
n_sim       <- 150      # Number of Monte Carlo replicates

# Derived
nT    <- round(nC * ratio)
ltfuc <- ltfuc_pct / 100
ltfut <- ltfut_pct / 100

cat("=== Hierarchical Win Ratio Simulation ===\n")
cat(sprintf("Study:   follow-up=%d mo | nC=%d | nT=%d (ratio %.2f:1)\n", fu, nC, nT, ratio))
cat(sprintf("LTFU:    ctrl=%.0f%% | trt=%.0f%% over study\n", ltfuc_pct, ltfut_pct))
cat(sprintf("Placebo: ACM=%.0f%% | CVH=%.2f/yr\n", acm_rate*100, cvh_rate))
cat(sprintf("Effect:  HR=%.2f (ACM) | RR=%.2f (CVH)\n", hr, rr))
cat(sprintf("Model:   Weibull nu=%.2f | NB kappa=%.2f | copula rho=%.2f\n", nu, kappa, rho))
cat(sprintf("Sims:    %d replicates x %d x %d = %s pairs\n\n",
            n_sim, nT, nC, format(n_sim*nT*nC, big.mark=",")))

# ---- HELPER FUNCTIONS -------------------------------------------------------

# Standard normal CDF (base R pnorm is fine)
# Gaussian copula: generate correlated (U_death, U_cvh) pair for one patient
# Z_death ~ N(0,1); Z_cvh = rho*Z_death + sqrt(1-rho^2)*Z2, Z2~N(0,1)
gen_copula_uniforms <- function(n, rho) {
  z1  <- rnorm(n)
  z2  <- rnorm(n)
  z_cvh <- rho * z1 + sqrt(1 - rho^2) * z2
  list(u_death = pnorm(z1), u_cvh = pnorm(z_cvh))
}

# Weibull inverse CDF:  t = scale * (-ln(1-u))^(1/shape)
# scale back-calculated so P(T <= fu) = acm_rate for control arm
weibull_inv <- function(u, scale, shape) {
  scale * (-log(1 - pmin(u, 1 - 1e-9)))^(1 / shape)
}

# Lomax (NB first-event) inverse CDF
# Var = mean + mean^2 * kappa  =>  t = 1/(lambda*kappa) * ((1-u)^{-kappa} - 1)
lomax_inv <- function(u, lambda, kappa) {
  (1 / (lambda * kappa)) * ((1 - pmin(u, 1 - 1e-9))^(-kappa) - 1)
}

# LTFU: uniform dropout. Patient lost if U_ltfu < ltfu_frac;
# dropout time = U_ltfu / ltfu_frac * fu
gen_ltfu_cens <- function(n, ltfu_frac, fu) {
  if (ltfu_frac <= 0) return(rep(Inf, n))
  u <- runif(n)
  ifelse(u < ltfu_frac, (u / ltfu_frac) * fu, Inf)
}

# ---- SINGLE REPLICATE -------------------------------------------------------
run_once <- function(fu, nC, nT, acm_rate, cvh_rate, hr, rr,
                     nu, kappa, rho, ltfuc, ltfut) {

  # Weibull scale: P(T <= fu | control) = acm_rate
  scale_ctrl <- fu / (-log(1 - acm_rate))^(1 / nu)
  scale_trt  <- scale_ctrl * hr^(-1 / nu)

  lambda_ctrl <- cvh_rate / 12   # per month
  lambda_trt  <- lambda_ctrl * rr

  # -- Control arm --
  cu_c   <- gen_copula_uniforms(nC, rho)
  t_dc   <- weibull_inv(cu_c$u_death, scale_ctrl, nu)
  t_ltfu_c <- gen_ltfu_cens(nC, ltfuc, fu)
  cens_c <- pmin(fu, t_ltfu_c)
  obs_dc <- pmin(t_dc, cens_c);  ev_dc <- as.integer(t_dc <= cens_c)
  ltfu_ctrl_n <- sum(t_ltfu_c < fu & t_dc > t_ltfu_c)

  t_hc   <- lomax_inv(cu_c$u_cvh, lambda_ctrl, kappa)
  obs_hc <- pmin(t_hc, obs_dc);  ev_hc <- as.integer(t_hc <= obs_dc)

  # -- Treatment arm --
  cu_t   <- gen_copula_uniforms(nT, rho)
  t_dt   <- weibull_inv(cu_t$u_death, scale_trt, nu)
  t_ltfu_t <- gen_ltfu_cens(nT, ltfut, fu)
  cens_t <- pmin(fu, t_ltfu_t)
  obs_dt <- pmin(t_dt, cens_t);  ev_dt <- as.integer(t_dt <= cens_t)
  ltfu_trt_n <- sum(t_ltfu_t < fu & t_dt > t_ltfu_t)

  t_ht   <- lomax_inv(cu_t$u_cvh, lambda_trt, kappa)
  obs_ht <- pmin(t_ht, obs_dt);  ev_ht <- as.integer(t_ht <= obs_dt)

  # -- Pairwise comparisons (nT x nC) --
  wd <- ld <- wh <- lh <- tf <- 0L

  for (i in seq_len(nT)) {
    for (j in seq_len(nC)) {
      # Level 1: ACM
      if      (ev_dc[j] == 1 && obs_dt[i] > obs_dc[j]) { wd <- wd + 1L }
      else if (ev_dt[i] == 1 && obs_dc[j] > obs_dt[i]) { ld <- ld + 1L }
      else {
        # Level 2: CVH
        if      (ev_hc[j] == 1 && obs_ht[i] > obs_hc[j]) { wh <- wh + 1L }
        else if (ev_ht[i] == 1 && obs_hc[j] > obs_ht[i]) { lh <- lh + 1L }
        else                                               { tf <- tf + 1L }
      }
    }
  }

  tot <- nT * nC
  list(
    wd = wd, ld = ld, wh = wh, lh = lh, tf = tf, tot = tot,
    cvh_ctrl  = mean(ev_hc),    cvh_trt  = mean(ev_ht),
    acm_ctrl  = mean(ev_dc),    acm_trt  = mean(ev_dt),
    ltfu_ctrl = ltfu_ctrl_n/nC, ltfu_trt = ltfu_trt_n/nT
  )
}

# ---- MAIN SIMULATION LOOP ---------------------------------------------------
set.seed(42)
res_list <- vector("list", n_sim)

cat("Running replicates")
for (s in seq_len(n_sim)) {
  res_list[[s]] <- run_once(fu, nC, nT, acm_rate, cvh_rate,
                             hr, rr, nu, kappa, rho, ltfuc, ltfut)
  if (s %% 30 == 0) cat(sprintf(" %d/%d", s, n_sim))
}
cat(" done.\n\n")

# ---- AGGREGATE RESULTS -------------------------------------------------------
acc <- list(
  wd  = sum(sapply(res_list, `[[`, "wd")),
  ld  = sum(sapply(res_list, `[[`, "ld")),
  wh  = sum(sapply(res_list, `[[`, "wh")),
  lh  = sum(sapply(res_list, `[[`, "lh")),
  tf  = sum(sapply(res_list, `[[`, "tf")),
  tot = sum(sapply(res_list, `[[`, "tot"))
)

pct  <- function(x, tot) round(100 * x / tot, 1)
WR   <- function(w, l)   if (l > 0) round(w / l, 3) else NA

tot  <- acc$tot
wdP  <- pct(acc$wd, tot);  ldP  <- pct(acc$ld, tot)
whP  <- pct(acc$wh, tot);  lhP  <- pct(acc$lh, tot)
tfP  <- pct(acc$tf, tot)
tdP  <- pct(acc$tf + acc$wh + acc$lh, tot)  # death-level ties (escalated)
wTot <- pct(acc$wd + acc$wh, tot)
lTot <- pct(acc$ld + acc$lh, tot)
dWR  <- WR(acc$wd, acc$ld)
hWR  <- WR(acc$wh, acc$lh)
oWR  <- WR(acc$wd + acc$wh, acc$ld + acc$lh)

obs_acm_ctrl  <- round(mean(sapply(res_list, `[[`, "acm_ctrl"))  * 100, 1)
obs_acm_trt   <- round(mean(sapply(res_list, `[[`, "acm_trt"))   * 100, 1)
obs_cvh_ctrl  <- round(mean(sapply(res_list, `[[`, "cvh_ctrl"))  * 100, 1)
obs_cvh_trt   <- round(mean(sapply(res_list, `[[`, "cvh_trt"))   * 100, 1)
obs_ltfu_ctrl <- round(mean(sapply(res_list, `[[`, "ltfu_ctrl")) * 100, 1)
obs_ltfu_trt  <- round(mean(sapply(res_list, `[[`, "ltfu_trt"))  * 100, 1)

# ---- PRINT RESULTS ----------------------------------------------------------
cat("=== Observed event rates ===\n")
cat(sprintf("  ACM:  placebo=%.1f%%  treatment=%.1f%%\n", obs_acm_ctrl, obs_acm_trt))
cat(sprintf("  CVH:  placebo=%.1f%%  treatment=%.1f%%\n", obs_cvh_ctrl, obs_cvh_trt))
cat(sprintf("  LTFU: ctrl=%.1f%%     trt=%.1f%%\n\n", obs_ltfu_ctrl, obs_ltfu_trt))

result_tbl <- data.frame(
  Level     = c("1. ACM", "2. CVH (among death ties)", "Overall (hierarchical)"),
  Wins_pct  = c(wdP,  whP,  wTot),
  Losses_pct= c(ldP,  lhP,  lTot),
  Ties_pct  = c(tdP,  tfP,  tfP),
  Win_Ratio = c(dWR,  hWR,  oWR)
)

cat("=== Win / Loss / Tie Summary ===\n")
print(result_tbl, row.names = FALSE)

cat(sprintf("\nOverall Win Ratio = %.3f\n", oWR))
cat(sprintf("  Death-level WR   = %.3f  (%.1f%% wins, %.1f%% losses)\n", dWR, wdP, ldP))
cat(sprintf("  CVH-level WR     = %.3f  (%.1f%% wins, %.1f%% losses)\n", hWR, whP, lhP))
cat(sprintf("  Final ties       = %.1f%%\n", tfP))

# =============================================================================
#  FINKELSTEIN-SCHOENFELD (F-S) STATISTIC & STUDY POWER
# =============================================================================
#
#  The F-S statistic generalises the log-rank test to hierarchical endpoints.
#  For each treatment-control pair (i,j), define the pairwise score:
#
#    D_ij = +1  if treatment wins  (i beats j in hierarchy)
#           -1  if treatment loses
#            0  if tie
#
#  The F-S test statistic is:
#    U   = sum_{i,j} D_ij          (= W - L, net wins)
#    V   = variance of U under H0
#    Z   = U / sqrt(V)  ~  N(0,1) under H0
#
#  Variance under H0 (Dong et al. 2020, Stat Med):
#    V = [nT*nC / (nT+nC)] * [mean(pi_i^2) + mean(pi_j^2)]  ... (complex)
#
#  Practical simulation-based approach used here:
#    - Per replicate, compute U_s = W_s - L_s and record it
#    - Under H0 (null simulation with HR=1, RR=1), compute V_0 = Var(U)
#    - Z_s = U_s / sqrt(V_0)   [standardised using null variance]
#    - Power = proportion of replicates where |Z_s| > z_alpha/2
#
#  We also compute the analytical null variance:
#    Under H0, E[D_ij] = 0, and V(U) = nT*nC*(p_conc)
#    where p_conc = P(pair is concordant) = 1 - P(tie) = (W+L) / (nT*nC)
#    Simplified: V_analytic = (W + L) per replicate (variance of Bernoulli scores)
# =============================================================================

alpha     <- 0.05     # two-sided significance level
z_alpha   <- qnorm(1 - alpha/2)

cat("\n=== Finkelstein-Schoenfeld Statistics & Power ===\n\n")

# ---- Per-replicate net score U_s = W_s - L_s --------------------------------
U_vec <- sapply(res_list, function(r) (r$wd + r$wh) - (r$ld + r$lh))
W_vec <- sapply(res_list, function(r)  r$wd + r$wh)
L_vec <- sapply(res_list, function(r)  r$ld + r$lh)
T_vec <- sapply(res_list, function(r)  r$tf)
N_pairs <- nT * nC

# ---- Z score per replicate using analytical null variance -------------------
# Under H0, concordant pairs contribute: Var(D_ij) = P(win)*1 + P(loss)*1
# => V_analytic_s = (W_s + L_s)  [sum of squared non-zero scores]
# Z_s = U_s / sqrt(W_s + L_s)
Z_analytic <- U_vec / sqrt(pmax(W_vec + L_vec, 1))

# ---- Z score using simulation-based null variance ---------------------------
# Run n_sim null replicates (HR=1, RR=1) to estimate V under H0
cat("Running null simulations (HR=1, RR=1) for empirical variance...\n")
set.seed(123)
null_list <- vector("list", n_sim)
for (s in seq_len(n_sim)) {
  null_list[[s]] <- run_once(fu, nC, nT, acm_rate, cvh_rate,
                              1.0, 1.0,   # HR=1, RR=1 (null)
                              nu, kappa, rho, ltfuc, ltfut)
}
U_null <- sapply(null_list, function(r) (r$wd + r$wh) - (r$ld + r$lh))
V_null <- var(U_null)       # empirical null variance
sd_null <- sqrt(V_null)

Z_empirical <- U_vec / sd_null

# ---- Power estimates --------------------------------------------------------
power_analytic  <- mean(abs(Z_analytic)  > z_alpha)
power_empirical <- mean(abs(Z_empirical) > z_alpha)

# ---- Z score summary --------------------------------------------------------
z_summary <- function(Z, label) {
  data.frame(
    Estimator  = label,
    Mean_Z     = round(mean(Z), 3),
    Median_Z   = round(median(Z), 3),
    SD_Z       = round(sd(Z), 3),
    P2.5       = round(quantile(Z, 0.025), 3),
    P97.5      = round(quantile(Z, 0.975), 3),
    Pct_sig    = round(100 * mean(abs(Z) > z_alpha), 1)
  )
}

z_tbl <- rbind(
  z_summary(Z_analytic,  "Analytical V (W+L)"),
  z_summary(Z_empirical, "Empirical V (null sim)")
)

cat("\n--- Z score summary across", n_sim, "replicates ---\n")
print(z_tbl, row.names = FALSE)

cat(sprintf("\n--- Power at alpha=%.2f (two-sided z > %.3f) ---\n", alpha, z_alpha))
cat(sprintf("  Power [analytical null var]:  %.1f%%\n", power_analytic  * 100))
cat(sprintf("  Power [empirical null var]:   %.1f%%\n", power_empirical * 100))

cat("\n--- Mean U (net wins) and null variance ---\n")
cat(sprintf("  Mean U (W-L) under alternative: %.1f  (SD=%.1f)\n", mean(U_vec), sd(U_vec)))
cat(sprintf("  Mean U (W-L) under null:        %.1f  (SD=%.1f)\n", mean(U_null), sd(U_null)))
cat(sprintf("  Empirical null SD:              %.1f\n", sd_null))
cat(sprintf("  Non-centrality parameter:       %.3f\n", mean(U_vec) / sd_null))

# ---- Distribution plot -------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  df_z <- data.frame(
    Z     = c(Z_empirical, U_null / sd_null),
    Group = rep(c("Alternative (simulated)", "Null (HR=RR=1)"), each = n_sim)
  )
  p_z <- ggplot(df_z, aes(x = Z, fill = Group, colour = Group)) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    geom_vline(xintercept = c(-z_alpha, z_alpha),
               linetype = "dashed", colour = "#e74c3c", linewidth = 0.7) +
    annotate("text", x = z_alpha + 0.15, y = Inf, vjust = 1.5, hjust = 0,
             label = sprintf("±%.2f", z_alpha), size = 3.2, colour = "#e74c3c") +
    scale_fill_manual(values = c("#3b82f6", "#94a3b8")) +
    scale_colour_manual(values = c("#1d4ed8", "#475569")) +
    labs(
      title    = "F-S Statistic: Z score distribution under null vs alternative",
      subtitle = sprintf("Power = %.1f%%  |  WR = %.3f  |  n=%d ctrl, %d trt  |  fu=%d mo",
                         power_empirical * 100, oWR, nC, nT, fu),
      x = "Z = U / SD(U_null)",
      y = "Density",
      fill = NULL, colour = NULL,
      caption = sprintf(
        "HR=%.2f (ACM)  RR=%.2f (CVH)  Weibull ν=%.2f  NB κ=%.2f  copula ρ=%.2f\nLTFU ctrl=%.0f%%  trt=%.0f%%  Placebo ACM=%.0f%%  CVH=%.2f/yr  %d replicates",
        hr, rr, nu, kappa, rho, ltfuc_pct, ltfut_pct, acm_rate*100, cvh_rate, n_sim)
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 9, colour = "grey40"),
          plot.caption  = element_text(size = 8, colour = "grey50", hjust = 0))
  print(p_z)
  ggsave("fs_power_zdist.png", plot = p_z, width = 9, height = 5.5, dpi = 150)
  message("Z distribution plot saved to fs_power_zdist.png")
} else {
  cat("\n(Install ggplot2 for the Z distribution plot)\n")
}
