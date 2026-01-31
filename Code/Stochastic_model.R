# ============================================================
# Stochastic decision-urgency inference from ecological shocks
# Real data version: State x Month with Rain, NDVI, Suitability, War
# Outputs:
#   - urgency_time_series_by_state.csv
#   - Figure_Urgency_Inference_by_state.png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
})

set.seed(10)

# -----------------------------
# USER SETTINGS
# -----------------------------
Tmax <- 1            # 1 year horizon
month_levels <- c("june","july","august","september")  # your data months
target_fraction <- 0.20  # policy: must keep >=50% of initial stock
#M0 initial stock index
df <- df %>%
  mutate(Biomass = ndvi_to_biomass(NDVI))

M0_by_state <- df %>%
  group_by(State) %>%
  summarise(M0 = mean(Biomass, na.rm = TRUE))
# shock index weights (sum to 1)
w_ops  <- 1/3
w_rain <- 1/3
w_T    <- 1/3

# sigma bounds
sigma_min <- 0.10
sigma_max <- 0.40

# gamma uncertainty (from your table / assumptions)
Qa_mean <- 0.80; Qa_sd <- 0.13
Qi_mean <- 0.50; Qi_sd <- 0.15
nsim_gamma <- 2000

# rain extreme thresholds
rain_low_q  <- 0.10
rain_high_q <- 0.95

# inversion numeric bounds
k_lower <- 1e-8
k_upper <- 200

# guardrails for y
y_eps <- 1e-6  # keep y in (eps, 1-eps)

# -----------------------------
# 1) Helpers: NDVI -> Biomass -> gamma(t)
# -----------------------------
ndvi_to_biomass <- function(ndvi) {
  R <- 11750 * ndvi - 3120
  pmax(R, 0)
}

gamma_holling <- function(R, a, h) {
  1 / (1 + a * h * R)
}

simulate_gamma_time <- function(ndvi_vec,
                                Mpred = 1, a0 = 1, I0 = 1,
                                Qa_mean = 0.80, Qa_sd = 0.13,
                                Qi_mean = 0.50, Qi_sd = 0.15,
                                nsim = 2000) {
  R <- ndvi_to_biomass(ndvi_vec)
  
  Qa <- rnorm(nsim, Qa_mean, Qa_sd)
  Qi <- rnorm(nsim, Qi_mean, Qi_sd)
  
  a    <- a0 * (Mpred ^ Qa)
  Imax <- I0 * (Mpred ^ Qi)
  h    <- 1 / Imax
  
  gamma_mat <- sapply(seq_along(R), function(i) gamma_holling(R[i], a, h))
  gamma_mat <- t(gamma_mat)
  
  tibble(
    step = seq_along(ndvi_vec),
    ndvi = ndvi_vec,
    biomass = R,
    gamma_mean = rowMeans(gamma_mat),
    gamma_p025 = apply(gamma_mat, 1, quantile, 0.025),
    gamma_p50  = apply(gamma_mat, 1, quantile, 0.50),
    gamma_p975 = apply(gamma_mat, 1, quantile, 0.975)
  )
}

# -----------------------------
# 2) Helpers: sigma(t) from ops + rain extremes + suitability
# -----------------------------
minmax01 <- function(x) {
  if (all(is.na(x))) return(x)
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

compute_sigma <- function(rain_mm, rain_percentile, suitability, war,
                          w_ops, w_rain, w_T, sigma_min, sigma_max,
                          low_q = 0.10, high_q = 0.95) {
  
  S_ops  <- ifelse(war == 1, 1, 0)
  S_rain <- ifelse(rain_percentile <= low_q | rain_percentile >= high_q, 1, 0)
  S_T    <- suitability  # expected in [0,1]
  
  S_index <- w_ops*S_ops + w_rain*S_rain + w_T*S_T
  
  sigma_min + (sigma_max - sigma_min) * S_index
}

# -----------------------------
# 3) Brownian multiplier Z(t) with time-varying sigma(t)
# -----------------------------
simulate_Zt <- function(sigma_vec, dt, seed = 10) {
  set.seed(seed)
  n <- length(sigma_vec)
  dW <- sqrt(dt) * rnorm(n)
  logZ <- numeric(n)
  acc <- 0
  for (i in 1:n) {
    acc <- acc - 0.5 * (sigma_vec[i]^2) * dt + sigma_vec[i] * dW[i]
    logZ[i] <- acc
  }
  exp(logZ)
}

# -----------------------------
# 4) Kernel + inversion for k(t)
# -----------------------------
G_kernel <- function(t, k, Tmax) {
  num <- 1 - exp(-k * (Tmax - t))
  den <- exp(k * t) - exp(-k * (Tmax - t))
  num / den
}

solve_k_from_y <- function(y, t, Tmax, k_lower = 1e-8, k_upper = 200) {
  if (!is.finite(y) || y <= 0 || y >= 1) return(NA_real_)
  f <- function(k) G_kernel(t, k, Tmax) - y
  
  fl <- f(k_lower); fu <- f(k_upper)
  if (!is.finite(fl) || !is.finite(fu)) return(NA_real_)
  
  if (fl * fu > 0) {
    for (ku in c(500, 1000, 2000)) {
      fu2 <- f(ku)
      if (is.finite(fu2) && fl * fu2 <= 0) {
        k_upper <- ku
        fu <- fu2
        break
      }
    }
    if (fl * fu > 0) return(NA_real_)
  }
  
  uniroot(f, lower = k_lower, upper = k_upper, tol = 1e-10)$root
}

# -----------------------------
# 5) Main state runner: infer r(t)
# -----------------------------
run_state_pipeline <- function(df_state,
                               Tmax, M0, target_fraction,
                               w_ops, w_rain, w_T,
                               sigma_min, sigma_max,
                               Qa_mean, Qa_sd, Qi_mean, Qi_sd, nsim_gamma,
                               rain_low_q, rain_high_q,
                               k_lower, k_upper, y_eps,
                               seed = 10) {
  
  n_steps <- nrow(df_state)
  dt <- Tmax / n_steps
  time <- seq(dt, Tmax, by = dt)
  
  # gamma(t) with uncertainty from NDVI
  gamma_sum <- simulate_gamma_time(
    ndvi_vec = df_state$NDVI,
    Qa_mean = Qa_mean, Qa_sd = Qa_sd,
    Qi_mean = Qi_mean, Qi_sd = Qi_sd,
    nsim = nsim_gamma
  )
  
  # rain percentile within state (already computed outside, but keep here if needed)
  # sigma(t)
  sigma_t <- compute_sigma(
    rain_mm = df_state$Rain_mm,
    rain_percentile = df_state$rain_pct,
    suitability = df_state$Suitability_01,
    war = df_state$War,
    w_ops = w_ops, w_rain = w_rain, w_T = w_T,
    sigma_min = sigma_min, sigma_max = sigma_max,
    low_q = rain_low_q, high_q = rain_high_q
  )
  
  # Z(t)
  Z_t <- simulate_Zt(sigma_t, dt, seed = seed)
  
  # policy target
  Mt_target <- rep(target_fraction * M0, n_steps)
  
  # inversion
  gamma_t <- gamma_sum$gamma_mean
  y_raw <- Mt_target / (M0 * Z_t)
  
  # clamp y into (eps, 1-eps) so inversion is defined
  y <- pmin(pmax(y_raw, y_eps), 1 - y_eps)
  
  k_hat <- sapply(seq_len(n_steps), function(i) {
    solve_k_from_y(y[i], time[i], Tmax, k_lower, k_upper)
  })
  
  r_hat <- ifelse(
    is.finite(k_hat),
    (1 - gamma_t) * k_hat + 0.5 * (sigma_t^2) * gamma_t * (gamma_t - 1),
    NA_real_
  )
  
  out <- tibble(
    State = df_state$State[1],
    Month = df_state$Month,
    step = seq_len(n_steps),
    time = time,
    Rain_mm = df_state$Rain_mm,
    NDVI = df_state$NDVI,
    Suitability = df_state$Suitability,
    War = df_state$War,
    rain_pct = df_state$rain_pct,
    gamma_mean = gamma_sum$gamma_mean,
    gamma_p025 = gamma_sum$gamma_p025,
    gamma_p975 = gamma_sum$gamma_p975,
    sigma = sigma_t,
    Z = Z_t,
    Mt_target = Mt_target,
    y_raw = y_raw,
    y = y,
    k = k_hat,
    r = r_hat,
    # diagnostic stock implied by local kernel
    G = ifelse(is.finite(k_hat), G_kernel(time, k_hat, Tmax), NA_real_),
    M_star = M0 * ifelse(is.finite(k_hat), G_kernel(time, k_hat, Tmax), NA_real_) * Z_t
  )
  
  out
}

# ============================================================
# 6) READ YOUR REAL DATA
# ============================================================
# Put your data into a CSV with columns:
# State, Month, Rain_mm, NDVI, Suitability, War
# Example path (change):
df <- read_csv('/Users/kagboka/Desktop/Method_paper_stocchastic DL/Clean_data-code_paper.shp/Climate_Rain_NDVI_Suitability_Admin1_2019_FINAL.csv')


# ============================================================
# 7) PREP + RUN ALL STATES
# ============================================================

# --- Example: assume you already loaded it as df ---
# df <- read_csv("...")

# Minimal cleaning + month order + suitability scaling + rain percentiles per state
prepare_inputs <- function(df, month_levels) {
  df %>%
    mutate(
      Month = tolower(str_trim(Month)),
      Month = factor(Month, levels = month_levels, ordered = TRUE),
      State = as.character(State)
    ) %>%
    arrange(State, Month) %>%
    group_by(State) %>%
    mutate(
      rain_pct = percent_rank(Rain_mm),          # within-state percentile (0..1)
      Suitability_01 = pmin(pmax(Suitability, 0), 1) # assume already 0..1; clamp anyway
    ) %>%
    ungroup()
}

# ============================================================
# 8) DRIVER (call this after df is loaded)
# ============================================================
run_all_states <- function(df) {
  df2 <- prepare_inputs(df, month_levels)
  
  out <- df2 %>%
    group_by(State) %>%
    group_split() %>%
    lapply(function(dfs) {
      run_state_pipeline(
        df_state = dfs,
        Tmax = Tmax, M0 = M0, target_fraction = target_fraction,
        w_ops = w_ops, w_rain = w_rain, w_T = w_T,
        sigma_min = sigma_min, sigma_max = sigma_max,
        Qa_mean = Qa_mean, Qa_sd = Qa_sd, Qi_mean = Qi_mean, Qi_sd = Qi_sd,
        nsim_gamma = nsim_gamma,
        rain_low_q = rain_low_q, rain_high_q = rain_high_q,
        k_lower = k_lower, k_upper = k_upper, y_eps = y_eps,
        seed = 10
      )
    }) %>%
    bind_rows()
  
  out
}

# ============================================================
# 9) PLOTTING (SEPARATE, PUBLICATION-READY FIGURES)
# ============================================================

plot_results_separate <- function(df_out) {
  
  theme_pub <- theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 9)  # reduce label size to avoid overlap
    )
  
  month_labs <- c(
    june = "Jun",
    july = "Jul",
    august = "Aug",
    september = "Sep"
  )
  
  # ------------------------------------------------------------
  # (a) Biological saturation gamma(t)
  # ------------------------------------------------------------
  p_gamma <- ggplot(df_out, aes(x = Month, y = gamma_mean, group = 1)) +
    geom_ribbon(aes(ymin = gamma_p025, ymax = gamma_p975), alpha = 0.2) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~State, scales = "free_y") +
    scale_x_discrete(labels = month_labs) +
    labs(
      x = NULL,
      y = expression(gamma(t)),
      title = ""
    ) +
    theme_pub
  
  # ------------------------------------------------------------
  # (b) Volatility sigma(t)
  # ------------------------------------------------------------
  p_sigma <- ggplot(df_out, aes(x = Month, y = sigma, group = 1)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    facet_wrap(~State, scales = "free_y") +
    scale_x_discrete(labels = month_labs) +
    labs(
      x = NULL,
      y = expression(sigma(t)),
      title = ""
    ) +
    theme_pub
  
  # ------------------------------------------------------------
  # (c) Inferred urgency r(t)
  # ------------------------------------------------------------
  p_r <- ggplot(df_out, aes(x = Month, y = r, group = 1)) +
    geom_line(linewidth = 0.9, na.rm = TRUE) +
    geom_point(size = 1.8, na.rm = TRUE) +
    facet_wrap(~State, scales = "free_y") +
    scale_x_discrete(labels = month_labs) +
    labs(
      x = NULL,
      y = "Inferred urgency r(t)",
      title = ""
    ) +
    theme_pub
  
  # ------------------------------------------------------------
  # (d) Diagnostic: implied stock vs policy target
  # ------------------------------------------------------------
  p_M <- ggplot(df_out, aes(x = Month, group = 1)) +
    geom_line(aes(y = M_star), linewidth = 0.9, na.rm = TRUE) +
    geom_point(aes(y = M_star), size = 1.8, na.rm = TRUE) +
    geom_line(aes(y = Mt_target), linetype = "dashed", linewidth = 0.9) +
    facet_wrap(~State, scales = "free_y") +
    scale_x_discrete(labels = month_labs) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      x = NULL,
      y = expression(M[t]),
      title = ""
    ) +
    theme_pub
  
  # ---- Explicit return (critical) ----
  return(list(
    gamma   = p_gamma,
    sigma   = p_sigma,
    urgency = p_r,
    stock   = p_M
  ))
}

# ============================================================
# 10) RUN + EXPORT FIGURES
# ============================================================

df_out <- run_all_states(df)
write_csv(df_out, "urgency_time_series_by_state.csv")

plots <- plot_results_separate(df_out)

# ---- Print to screen ----
print(plots$gamma)
print(plots$sigma)
print(plots$urgency)
print(plots$stock)

# ---- Save each figure separately ----
ggsave("Figure_gamma_by_state.png",   plots$gamma,   width = 12, height = 8, dpi = 600)
ggsave("Figure_sigma_by_state.png",   plots$sigma,   width = 12, height = 8, dpi = 600)
ggsave("Figure_urgency_r_by_state.png", plots$urgency, width = 12, height = 8, dpi = 600)
ggsave("Figure_stock_vs_target_by_state.png", plots$stock, width = 12, height = 8, dpi = 600)

# ============================================================
# 11) SPATIO–TEMPORAL URGENCY MAPPING (ADM1 / ADM2)
# ============================================================

library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(scales)

# -----------------------------
# 11.1 Load administrative boundaries
# -----------------------------
adm1 <- st_read('/Users/kagboka/Desktop/Method_paper_stocchastic DL/sdn_admin_boundaries.shp/sdn_admin1.shp')   # State-level
adm2 <- st_read('/Users/kagboka/Desktop/Method_paper_stocchastic DL/sdn_admin_boundaries.shp/sdn_admin2.shp')   # District-level (optional)

# Harmonize state naming
adm1 <- adm1 %>% mutate(State = adm1_name)
adm2 <- adm2 %>% mutate(State = adm1_name)

# -----------------------------
# 11.2 Aggregate urgency metrics by State
# -----------------------------
map_data <- df_out %>%
  group_by(State) %>%
  summarise(
    r_mean     = mean(r, na.rm = TRUE),
    r_max      = max(r, na.rm = TRUE),
    sigma_mean = mean(sigma, na.rm = TRUE),
    gamma_mean = mean(gamma_mean, na.rm = TRUE),
    .groups = "drop"
  )

# Join to geometries
adm1_map <- adm1 %>% left_join(map_data, by = "State")
adm2_map <- adm2 %>% left_join(map_data, by = "State")

# -----------------------------
# 11.3 Generic mapping function
# -----------------------------
plot_map <- function(sf_obj, var, title, fill_lab) {
  ggplot(sf_obj) +
    geom_sf(aes(fill = .data[[var]]),
            color = "black", linewidth = 0.25) +
    scale_fill_viridis_c(
      option = "plasma",
      na.value = "grey90"
    ) +
    labs(title = title, fill = fill_lab) +
    theme_void(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}

# -----------------------------
# 11.4 ADM1-level maps
# -----------------------------
p_r_mean <- plot_map(
  adm1_map, "r_mean",
  "Mean Decision Urgency (r̄)",
  "r̄"
)

p_r_max <- plot_map(
  adm1_map, "r_max",
  "Peak Decision Urgency (max r)",
  "max r"
)

p_sigma <- plot_map(
  adm1_map, "sigma_mean",
  "Mean Shock Volatility (σ̄)",
  "σ̄"
)

p_gamma <- plot_map(
  adm1_map, "gamma_mean",
  "Mean Biological Saturation (γ̄)",
  "γ̄"
)

# -----------------------------
# 11.5 Combined spatial atlas (ADM1)
# -----------------------------
urgency_atlas_adm1 <- (p_r_mean | p_r_max) /
  (p_sigma  | p_gamma)

print(urgency_atlas_adm1)

ggsave(
  "Figure_Spatial_Urgency_Atlas_ADM1.png",
  urgency_atlas_adm1,
  width = 12,
  height = 10,
  dpi = 600
)


# ============================================================
# 11.5 ADM1 visualization (DECISION MAP – urgency classes)
# ============================================================

# --- Compute quantile thresholds ONLY for mapping ---
R_vals <- df_out$r
R_vals <- R_vals[is.finite(R_vals)]

r_Q25 <- quantile(R_vals, 0.25)
r_Q50 <- quantile(R_vals, 0.50)
r_Q75 <- quantile(R_vals, 0.75)

# --- Aggregate urgency by ADM1 (State) ---
map_data_adm1 <- df_out %>%
  group_by(State) %>%
  summarise(
    r_mean = mean(r, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Urgency_class = case_when(
      r_mean <= r_Q25                 ~ "Low",
      r_mean > r_Q25 & r_mean <= r_Q50 ~ "Moderate",
      r_mean > r_Q50 & r_mean <= r_Q75 ~ "High",
      r_mean > r_Q75                  ~ "Critical"
    )
  )

# --- Join to ADM1 geometry ---
adm1_map_class <- adm1 %>% 
  left_join(map_data_adm1, by = "State")

# --- Plot categorical decision map (ADM1) ---
p_adm1_r_class <- ggplot(adm1_map_class) +
  geom_sf(
    aes(fill = Urgency_class),
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "Low"      = "#2c7bb6",
      "Moderate" = "#abd9e9",
      "High"     = "#fdae61",
      "Critical" = "#d7191c"
    ),
    drop = FALSE
  ) +
  labs(
    title = "",
    fill  = "Urgency level"
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(p_adm1_r_class)

ggsave(
  "Figure_Spatial_Urgency_ADM1_DecisionClasses.png",
  p_adm1_r_class,
  width = 12,
  height = 8,
  dpi = 600
)
