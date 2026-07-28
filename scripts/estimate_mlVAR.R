# ============================================================================
#
#  Bivariate VAR(1) with ADHD (time-invariant), Social Context (time-varying),
#  and Weekend (time-varying) — Bayesian multilevel via brms
#
# ============================================================================

# ============================================================================
# 1. LOAD LIBRARIES
# ============================================================================
library(readxl)
library(lubridate)
library(tidyverse)
library(dplyr)
library(tidyr)
library(purrr)
library(brms)
library(ggplot2)
library(patchwork)
library(ggdist)
library(bayesplot)
library(ggpubr)

# ============================================================================
# 2. LOAD RAW DATA
# ============================================================================
in_dat <- read_xlsx(
  "./data/SPIG EMA updated data.xlsx"
)

in_dat$Response[in_dat$Response == "#skipped#"] <- NA
names(in_dat) <- make.names(names(in_dat))

# ============================================================================
# 3. PIVOT TO WIDE FORMAT (PER PROMPT)
# ============================================================================
all_data <- map_dfr(unique(in_dat$Participant.ID), function(i) {
  iso_dat <- in_dat %>% filter(Participant.ID == i)
  
  wide_iso_dat <- pivot_wider(
    iso_dat,
    id_cols     = c(Notification.Time, Session.Instance.No, ADHD),
    names_from  = Prompt.Label,
    values_from = Response
  ) %>%
    distinct() %>%
    mutate(ID = i) %>%
    select(
      Notification.Time, Session.Instance.No, ID, ADHD,
      adhdig_ema_1, adhdig_ema_2, adhdig_ema_3, adhdig_ema_4,
      adhdig_ema_5, adhdig_ema_6, adhdig_ema_7, adhdig_ema_8,
      `adhdig_ema_9_I am alone`,
      `adhdig_ema_9_Through social media`,
      `adhdig_ema_9_In-person`
    )
  
  wide_iso_dat
})


# ============================================================================
# 3b. PARSE DATETIME, CREATE WEEKEND FLAG, ORDER & CREATE OBS COUNTER
# ============================================================================

all_data <- all_data %>%
  filter(!is.na(Session.Instance.No)) %>%
  mutate(
    Notification.Time = ymd_hms(Notification.Time, tz = "UTC"),
    day_of_week       = wday(Notification.Time, label = TRUE),
    weekend           = as.integer(wday(Notification.Time) %in% c(1, 7))
  )

# Check for parsing failures
n_bad <- sum(is.na(all_data$Notification.Time))
if (n_bad > 0) {
  cat(sprintf("  WARNING: %d rows failed datetime parsing\n", n_bad))
}


# Convert to numeric AFTER extracting weekend
all_data <- all_data %>%
  mutate(Notification.Time = as.numeric(Notification.Time))

# Sort by ID then time, create observation counter
all_data <- all_data %>%
  arrange(ID, Notification.Time) %>%
  group_by(ID) %>%
  mutate(obs = row_number()) %>%
  ungroup()

# ============================================================================
# 5. RENAME COLUMNS
# ============================================================================

all_data <- all_data %>%
  rename(
    time_col     = Notification.Time,
    SessionNumb  = Session.Instance.No,
    id           = ID,
    item1        = adhdig_ema_1,
    item2        = adhdig_ema_2,
    item3        = adhdig_ema_3,
    item4        = adhdig_ema_4,
    item5        = adhdig_ema_5,
    item6        = adhdig_ema_6,
    item7        = adhdig_ema_7,
    item8        = adhdig_ema_8,
    Alone        = `adhdig_ema_9_I am alone`,
    SocialMedia  = `adhdig_ema_9_Through social media`,
    InPerson     = `adhdig_ema_9_In-person`
  )

# ============================================================================
# 5b. Create response pattern diff data
# ============================================================================
rep_diff_long <- all_data %>%
  rename(ID = id) %>%
  mutate(
    item1_present = as.numeric(!is.na(item1)),
    item2_present = as.numeric(!is.na(item2)),
    item3_present = as.numeric(!is.na(item3)),
    item4_present = as.numeric(!is.na(item4)),
    item5_present = as.numeric(!is.na(item5)),
    item6_present = as.numeric(!is.na(item6)),
    item7_present = as.numeric(!is.na(item7)),
    item8_present = as.numeric(!is.na(item8))
  ) %>% 
  group_by(ID) %>%
  arrange(ID, time_col) %>%
  filter(rowSums(!is.na(across(item1:item8))) > 0) %>%
  mutate(
    RawHours    = as.numeric(difftime(time_col, first(time_col), units = "hours")),
    gap_hours   = RawHours - lag(RawHours),
    new_session = is.na(gap_hours) | gap_hours > SESSION_GAP_HOURS,
    session     = cumsum(new_session)
  ) %>%
  ungroup() %>%
  pivot_longer(
    cols = starts_with("item") & ends_with("_present"),
    names_to = "item",
    values_to = "responded",
    names_pattern = "item(\\d+)_present"
  ) %>%
  mutate(item = as.numeric(item))

## Now run a glmer testing for temporal response pattern differences
#mod_rep1 = brm(responded ~ item + RawHours * ADHD + (1 + RawHours|ID), data = rep_diff_long, cores = 4)
#mod_rep2 = brm(responded ~ item + s(RawHours) + ADHD + (1 + RawHours|ID), data = rep_diff_long)
mod_rep1P = lme4::glmer(responded ~ item + RawHours * ADHD + (1 + RawHours|ID), data = rep_diff_long)
mod_rep1P = mgcv::gamm(responded ~ item + s(RawHours) + s(RawHours, by = ADHD), data = rep_diff_long, family="binomial")

# ============================================================================
# 6. COMPUTE TIME, SESSIONS (30-HOUR GAP RULE), AND BASIC VARIABLES
# ============================================================================
SESSION_GAP_HOURS <- 30

ema_data <- all_data %>%
  rename(ID = id) %>%
  mutate(across(item1:item8, as.numeric),
         across(c(Alone, SocialMedia, InPerson), as.numeric)) %>%
  group_by(ID) %>%
  arrange(ID, time_col) %>%
  filter(rowSums(!is.na(across(item1:item8))) > 0) %>%
  mutate(
    RawHours    = as.numeric(difftime(time_col, first(time_col), units = "hours")),
    gap_hours   = RawHours - lag(RawHours),
    new_session = is.na(gap_hours) | gap_hours > SESSION_GAP_HOURS,
    session     = cumsum(new_session)
  ) %>%
  ungroup()

# Session diagnostics
session_summary <- ema_data %>%
  group_by(ID, session) %>%
  summarise(n_obs = n(), .groups = "drop")

gap_summary <- ema_data %>%
  filter(!is.na(gap_hours)) %>%
  summarise(
    mean_gap   = mean(gap_hours),
    median_gap = median(gap_hours),
    max_gap    = max(gap_hours),
    n_large    = sum(gap_hours > SESSION_GAP_HOURS)
  )

# Within-session time
ema_data <- ema_data %>%
  group_by(ID, session) %>%
  mutate(
    TimeHours     = as.numeric(difftime(time_col, first(time_col), units = "hours")),
    time_interval = TimeHours - lag(TimeHours, default = 0)
  ) %>%
  ungroup()

# ============================================================================
# 7. EXTRACT TIME-INVARIANT COVARIATES (ADHD STATUS)
# ============================================================================
time_invar <- ema_data %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup() %>%
  select(ID, ADHD) %>%
  mutate(ADHD_numeric = ifelse(ADHD == "Y", 1, 0))

time_invar$ADHD_numeric[is.na(time_invar$ADHD_numeric)] <- 0
## Now add the other assessments of pathology here too
add_covar = in_dat <- read_xlsx(
  "./data/pennStateEMAAnalyses/ADHD analyses R21/possible covariates .xlsx"
)
spig_key = read_csv("./data/pennStateEMAAnalyses/ADHD analyses R21/spig_id_key.csv")

## Now merge
colnames(add_covar)[1] = "SPIG_ID"
colnames(spig_key) = c("ID", "SPIG_ID")
add_covar = merge(spig_key, add_covar)

time_invar = merge(time_invar, add_covar, by = "ID")
## Now fix the names
names(time_invar) <- gsub(" ", "_", names(time_invar))

cat(sprintf("  ADHD: %d,  Control: %d\n",
            sum(time_invar$ADHD_numeric == 1),
            sum(time_invar$ADHD_numeric == 0)))
psych::describe(time_invar)

# ============================================================================
# 8. FILTER: MINIMUM OBSERVATIONS + VARIANCE CHECK
# ============================================================================
sufficient_persons <- ema_data %>%
  count(ID) %>%
  filter(n >= 0) %>%
  pull(ID)

ema_data <- ema_data %>% filter(ID %in% sufficient_persons)

cat(sprintf("  After >=15 obs filter: %d obs, %d persons\n",
            nrow(ema_data), length(unique(ema_data$ID))))

variance_flags <- ema_data %>%
  group_by(ID) %>%
  summarise(
    var_items1_4  = var(c(item1, item2, item3, item4), na.rm = TRUE),
    var_items5_8  = var(c(item5, item6, item7, item8), na.rm = TRUE),
    flag_any      = (var_items1_4 < 0 | is.na(var_items1_4)) |
      (var_items5_8 < 0 | is.na(var_items5_8)),
    .groups       = "drop"
  )

ema_data <- ema_data %>%
  left_join(variance_flags %>% select(ID, flag_any), by = "ID") %>%
  filter(!flag_any) %>%
  select(-flag_any)

cat(sprintf("  After variance filter: %d obs, %d persons\n",
            nrow(ema_data), length(unique(ema_data$ID))))

session_summary <- ema_data %>%
  group_by(ID, session) %>%
  summarise(n_obs = n(), .groups = "drop")

# ============================================================================
# 9. COMPUTE MEAN SCORE
# ============================================================================

ema_data <- ema_data %>%
  mutate(
    anti_mean = rowMeans(cbind(item1, item2, item3, item4), na.rm = TRUE),
    pro_mean  = rowMeans(cbind(item5, item6, item7, item8), na.rm = TRUE),
    anti_mean = ifelse(is.nan(anti_mean), NA_real_, anti_mean),
    pro_mean  = ifelse(is.nan(pro_mean),  NA_real_, pro_mean)
  )

# Grand-mean center and standard
anti_gm        <- mean(ema_data$anti_mean, na.rm = TRUE)
pro_gm         <- mean(ema_data$pro_mean,  na.rm = TRUE)
anti_sd_global <- sd(ema_data$anti_mean,   na.rm = TRUE)
pro_sd_global  <- sd(ema_data$pro_mean,    na.rm = TRUE)

ema_data <- ema_data %>%
  mutate(
    Y1 = (anti_mean - anti_gm) / anti_sd_global,
    Y2 = (pro_mean  - pro_gm)  / pro_sd_global
  )

cat(sprintf("  Grand means: anti=%.3f, pro=%.3f\n", anti_gm, pro_gm))


# ============================================================================
# 10. PROCESS BINARY SOCIAL CONTEXT INDICATORS
# ============================================================================
context_vars <- c("Alone", "SocialMedia", "InPerson")

for (vname in context_vars) {
  n_na <- sum(is.na(ema_data[[vname]]))
  if (n_na > 0) cat(sprintf("  %s: %d NAs → 0\n", vname, n_na))
  ema_data[[vname]][is.na(ema_data[[vname]])] <- 0
  ema_data[[vname]] <- as.integer(ema_data[[vname]] != 0)
}

# ---- Context frequency table ----
ema_data %>%
  count(Alone, SocialMedia, InPerson, name = "n") %>%
  arrange(desc(n)) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

# ---- Create ADHD × Context interaction terms ----
## First pick the variables that best represent our intern and extern dims
test_mod = psych::fa(time_invar[,-c(1:5)], nfactors = 2)

if (!"ADHD_numeric" %in% names(ema_data)) {
  ema_data <- ema_data %>%
    left_join(time_invar %>% select(ID, ADHD_numeric, Age, C_SCARED_Score, DBD_Cp_Scr), by = "ID")
  ema_data$ADHD_numeric[is.na(ema_data$ADHD_numeric)] <- 0
}

interaction_td_names <- paste0("ADHD_x_", base_td_names)

for (k in seq_along(base_td_names)) {
  ema_data[[interaction_td_names[k]]] <-
    as.integer(ema_data$ADHD_numeric * ema_data[[base_td_names[k]]])
}

# Also create ADHD × weekend interaction
ema_data <- ema_data %>%
  mutate(ADHD_x_weekend = as.integer(ADHD_numeric * weekend))

td_var_names <- c(base_td_names, interaction_td_names)
n_TDpred     <- length(td_var_names)

# ---- Descriptives by ADHD group ----
desc_vars <- c(td_var_names, "weekend", "ADHD_x_weekend")

ema_data %>%
  group_by(ADHD_numeric) %>%
  summarise(
    n_obs = n(),
    across(all_of(desc_vars), ~ round(100 * mean(.x, na.rm = TRUE), 1),
           .names = "pct_{.col}"),
    .groups = "drop"
  ) %>%
  print()

# ---- Within-person variation check ----
person_variation <- ema_data %>%
  group_by(ID) %>%
  summarise(
    across(all_of(c(td_var_names, "weekend")),
           list(sd = ~ sd(.x, na.rm = TRUE)),
           .names = "{.col}_sd"),
    .groups = "drop"
  )

# ============================================================================
# 11. CREATE SEQUENTIAL PERSON IDs
# ============================================================================
id_mapping <- ema_data %>%
  distinct(ID) %>%
  arrange(ID) %>%
  mutate(person_id = row_number())

ema_data <- ema_data %>%
  select(-any_of("person_id")) %>%
  left_join(id_mapping, by = "ID")

# ============================================================================
# 12. HANDLE SESSION BREAKS
# ============================================================================
MIN_OBS_PER_SESSION <- 3

ema_data <- ema_data %>%
  mutate(ctsem_id = person_id * 1000L + session) %>%
  group_by(ctsem_id) %>%
  arrange(ctsem_id, time_col) %>%
  mutate(
    ctsem_time = as.numeric(
      difftime(time_col, first(time_col), units = "hours")
    )
  ) %>%
  ungroup()

session_obs <- ema_data %>% count(ctsem_id)

n_before  <- n_distinct(ema_data$ctsem_id)
valid_ids <- session_obs %>% filter(n >= MIN_OBS_PER_SESSION) %>% pull(ctsem_id)
ema_data  <- ema_data %>% filter(ctsem_id %in% valid_ids)

# Verify no within-session gaps exceed threshold
within_check <- ema_data %>%
  group_by(ctsem_id) %>%
  arrange(ctsem_time) %>%
  mutate(dt = ctsem_time - lag(ctsem_time)) %>%
  filter(!is.na(dt)) %>%
  summarise(max_dt = max(dt), .groups = "drop")


# ============================================================================
# 12b. HANDLE ZERO TIME INTERVALS
# ============================================================================
ema_data <- ema_data %>%
  arrange(ctsem_id, ctsem_time) %>%
  group_by(ctsem_id) %>%
  mutate(
    dt_check   = ctsem_time - lag(ctsem_time),
    is_zero_dt = !is.na(dt_check) & dt_check == 0
  ) %>%
  ungroup()

n_zero <- sum(ema_data$is_zero_dt)

if (n_zero > 0) {
  # Remove exact duplicates
  ema_data <- ema_data %>%
    arrange(ctsem_id, ctsem_time, time_col) %>%
    group_by(ctsem_id) %>%
    mutate(
      is_exact_dup = is_zero_dt &
        (is.na(Y1) & is.na(lag(Y1)) |
           (!is.na(Y1) & !is.na(lag(Y1)) & Y1 == lag(Y1))) &
        (is.na(Y2) & is.na(lag(Y2)) |
           (!is.na(Y2) & !is.na(lag(Y2)) & Y2 == lag(Y2)))
    ) %>%
    ungroup()
  
  n_exact_dup <- sum(ema_data$is_exact_dup, na.rm = TRUE)
  cat(sprintf("  Exact duplicates removed: %d\n", n_exact_dup))
  if (n_exact_dup > 0) {
    ema_data <- ema_data %>% filter(!is_exact_dup | is.na(is_exact_dup))
  }
}

ema_data <- ema_data %>%
  select(-any_of(c("dt_check", "is_zero_dt", "is_exact_dup", "zero_run")))

# Final verification
final_dt_check <- ema_data %>%
  arrange(ctsem_id, ctsem_time) %>%
  group_by(ctsem_id) %>%
  mutate(dt = ctsem_time - lag(ctsem_time)) %>%
  filter(!is.na(dt)) %>%
  ungroup()


# ============================================================================
# 13. PREPARE brms DATA
# ============================================================================

brms_dat <- ema_data %>%
  group_by(ID) %>%
  mutate(
    Y1_pm = mean(Y1, na.rm = TRUE),
    Y2_pm = mean(Y2, na.rm = TRUE),
    Y1_wc = Y1 - Y1_pm,
    Y2_wc = Y2 - Y2_pm
  ) %>%
  ungroup()

  arrange(ctsem_id, ctsem_time) %>%
  group_by(ctsem_id) %>%
  mutate(
    # ── Lags of within-centered scores (pure within-person dynamics) ─────
    Y1_lag_wc = lag(Y1_wc),
    Y2_lag_wc = lag(Y2_wc),
    
    # ── Person-mean lags (between-person level control; constant within
    #    person but needs to exist on lag rows) ────────────────────────────
    Y1_pm_lag = lag(Y1_pm),
    Y2_pm_lag = lag(Y2_pm),
    
    # ── Time interval ────────────────────────────────────────────────────
    dt = ctsem_time - lag(ctsem_time)
  ) %>%
  ungroup() %>%
  filter(
    !is.na(Y1_lag_wc),
    !is.na(Y2_lag_wc),
    !is.na(dt),
    dt > 0
  ) %>%
  mutate(
    log_dt    = log(dt),
    person_id = factor(ID),
    ADHD_x_Y1_lag  = ADHD_numeric * Y1_lag_wc,
    ADHD_x_Y2_lag  = ADHD_numeric * Y2_lag_wc,
    ADHD_x_SocialMedia = ADHD_numeric * SocialMedia,
    ADHD_x_InPerson    = ADHD_numeric * InPerson,
    ADHD_x_Alone       = ADHD_numeric * Alone,
    ADHD_x_weekend     = ADHD_numeric * weekend,
    ADHD_x_log_dt = ADHD_numeric * log_dt
  )

# ============================================================================
# 13b. AR(2)-COMPATIBLE DATA (used for BOTH models so comparison is fair)
# ============================================================================

brms_dat_ar2 <- ema_data %>%
  arrange(ctsem_id, ctsem_time) %>%
  group_by(ctsem_id) %>%
  mutate(
    Y1_lag  = lag(Y1),
    Y1_lag2 = lag(Y1, 2),
    Y2_lag  = lag(Y2),
    Y2_lag2 = lag(Y2, 2),
    dt      = ctsem_time - lag(ctsem_time)
  ) %>%
  ungroup() %>%
  filter(
    !is.na(Y1_lag),  !is.na(Y2_lag),
    !is.na(Y1_lag2), !is.na(Y2_lag2),   # ← KEY: require 2nd lag
    !is.na(Y1),      !is.na(Y2),
    !is.na(dt),      dt > 0
  ) %>%
  mutate(
    log_dt    = log(dt),
    person_id = factor(person_id)
  )

# ============================================================================
# 14. DEFINE brms MULTIVARIATE VAR(1) MODEL (with weekend)
# ============================================================================

# ============================================================================
# 14a. AR(1) MODEL — refit on AR(2)-compatible data for fair comparison
# ============================================================================

bf_Y1_ar1 <- bf(
  Y1 ~ 1 +
    Y1_lag + Y2_lag +
    (1 + Y1_lag + log_dt | p | person_id)
) + gaussian()

bf_Y2_ar1 <- bf(
  Y2 ~ 1 +
    Y1_lag + Y2_lag +
    (1 + Y2_lag + log_dt | p | person_id)
) + gaussian()

fit_ar1 <- brm(
  bf_Y1_ar1 + bf_Y2_ar1 + set_rescor(TRUE),
  data    = brms_dat_ar2,   # ← same data as AR(2)
  cores   = 4,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 42,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)


# ============================================================================
# 14b. AR(2) MODEL
# ============================================================================

bf_Y1_ar2 <- bf(
  Y1 ~ 1 +
    # ── AR(1) + AR(2) lags ──
    Y1_lag + Y1_lag2 + Y2_lag + Y2_lag2 +
    (1 + Y1_lag + Y1_lag2 + Y2_lag + Y2_lag2 + log_dt | p | person_id)
) + gaussian()

bf_Y2_ar2 <- bf(
  Y2 ~ 1 +
    Y1_lag + Y1_lag2 + Y2_lag + Y2_lag2 +
    (1 +Y1_lag + Y1_lag2 + Y2_lag + Y2_lag2 + log_dt | p | person_id)
) + gaussian()

fit_ar2 <- brm(
  bf_Y1_ar2 + bf_Y2_ar2 + set_rescor(TRUE),
  data    = brms_dat_ar2,
  cores   = 4,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 42,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

# ============================================================================
# 14c. MODEL COMPARISON: AR(1) vs AR(2)
# ============================================================================

# Add LOO criterion to each model
loo_ar1 <- add_criterion(fit_ar1, "loo")
loo_ar2 <- add_criterion(fit_ar2, "loo")

# Compare
loo_ar1 = loo(fit_ar1)
loo_ar2 = loo(fit_ar2)
loo_comp <- loo_compare(loo_ar1, loo_ar2)
print(loo_comp)

best_model <- rownames(loo_comp)[1]
elpd_diff  <- loo_comp[2, "elpd_diff"]
se_diff    <- loo_comp[2, "se_diff"]

# ============================================================================
# 15. FIT THE MODEL
# ============================================================================
bf_Y1 <- bf(
  Y1 ~ 1 +
    # Temporal dynamics
    Y1_lag + Y2_lag +
    # Between-person
    ADHD_numeric +
    # ADHD moderates dynamics
    Y1_lag:ADHD_numeric + Y2_lag:ADHD_numeric +
    # Social context
    SocialMedia + InPerson + Alone +
    # ADHD × Context
    ADHD_x_SocialMedia + ADHD_x_InPerson + ADHD_x_Alone +
    # ── Weekend ──
    weekend +
    # Internalizing
    C_SCARED_Score +
    # Extern
    DBD_Cp_Scr +
    # Age 
    Age +
    # Spacing correction
    log_dt + log_dt:ADHD_numeric +
    # Random effects
    (1 + Y2_lag + Y1_lag + log_dt | p | person_id)
) + gaussian()

bf_Y2 <- bf(
  Y2 ~ 1 +
    Y1_lag + Y2_lag +
    ADHD_numeric +
    Y1_lag:ADHD_numeric + Y2_lag:ADHD_numeric +
    SocialMedia + InPerson + Alone +
    ADHD_x_SocialMedia + ADHD_x_InPerson + ADHD_x_Alone +
    weekend  +
    C_SCARED_Score +
    DBD_Cp_Scr +
    Age +
    log_dt + log_dt:ADHD_numeric +
    (1 + Y2_lag + Y1_lag + log_dt | p | person_id)
) + gaussian()

fit_brms <- brm(
  bf_Y1 + bf_Y2 + set_rescor(TRUE),
  data    = brms_dat,
  prior   = model_priors,
  cores   = 4,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 42,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)