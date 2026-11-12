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
SESSION_GAP_HOURS = 30
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
## Scale before estimation
rep_diff_long$RawHours_s = scale(rep_diff_long$RawHours)[,]
rep_diff_long$item = factor(rep_diff_long$item)
mod_rep1 = brm(responded ~ item + RawHours_s * ADHD + weekend + (1 + RawHours_s|ID), 
               data = rep_diff_long, cores = 4, backend = "cmdstanr", family = "bernoulli",
               file = "./data/missingLogitModel.rds")
## Write summary table
write.csv(round(summary(mod_rep1)$fixed,2), file = "./data/out_endorse_table.csv", quote = F, row.names = T)

#mod_rep2 = brm(responded ~ item + s(RawHours) + ADHD + (1 + RawHours|ID), data = rep_diff_long)
mod_rep1P = lme4::glmer(responded ~ item + RawHours_s * ADHD + weekend + (1 + RawHours_s|ID),
                        data = rep_diff_long, family = "binomial")
#mod_rep1P = gamm(responded ~ item + s(RawHours) + ADHD + s(RawHours, by = "item"), data = rep_diff_long, family="binomial")

# ============================================================================
# 6. COMPUTE TIME, SESSIONS (30-HOUR GAP RULE), AND BASIC VARIABLES
# ============================================================================
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
  "./data/possible covariates .xlsx"
)
spig_key = read_csv("./data/spig_id_key.csv")

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

ema_data_with_counts <- ema_data %>%
  mutate(
    # Count how many items are missing in the anti scale (item1 to item4)
    missing_anti = rowSums(is.na(select(., item1:item4))),
    
    # Count how many items are missing in the pro scale (item5 to item8)
    missing_pro  = rowSums(is.na(select(., item5:item8))),
    
    # Count total missing items across all 8 items combined
    missing_total = rowSums(is.na(select(., item1:item8)))
  )

# To see a breakdown/frequency table of how many observations have 0, 1, 2, etc. missing items:
table(ema_data_with_counts$missing_total)

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

interaction_td_names <- paste0("ADHD_x_", context_vars)
base_td_names = context_vars
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
    Y2_wc = Y2 - Y2_pm,
    ## Now make the person centered binary social indicator variables
    alone_pm = mean(Alone, na.rm = TRUE),
    socia_pm = mean(SocialMedia, na.rm = TRUE),
    inper_pm = mean(InPerson, na.rm = TRUE),
    alone_wc = Alone - alone_pm,
    socia_wc = SocialMedia - socia_pm,
    inper_wc = InPerson - inper_pm
  ) %>%
  ungroup() |>
  arrange(ctsem_id, ctsem_time) %>%
  group_by(ctsem_id) %>%
  mutate(
    # ── Lags of within-centered scores (pure within-person dynamics) ─────
    Y1_lag_wc = lag(Y1_wc),
    Y2_lag_wc = lag(Y2_wc),
    # ── Person-mean lags (between-person level control; constant within
    # person but needs to exist on lag rows) ────────────────────────────
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
    log_dt = log(dt),
    person_id = factor(ID),
    ADHD_x_Y1_lag = ADHD_numeric * Y1_lag_wc,
    ADHD_x_Y2_lag = ADHD_numeric * Y2_lag_wc,
    ADHD_x_SocialMedia = ADHD_numeric * SocialMedia,
    ADHD_x_InPerson = ADHD_numeric * InPerson,
    ADHD_x_Alone = ADHD_numeric * Alone,
    ADHD_x_weekend = ADHD_numeric * weekend,
    ADHD_x_SocialMediawc = ADHD_numeric * socia_wc,
    ADHD_x_InPersonwc = ADHD_numeric * inper_wc,
    ADHD_x_Alonewc = ADHD_numeric * alone_wc,
    ADHD_x_weekend = ADHD_numeric * weekend,
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
    !is.na(Y1_lag2), !is.na(Y2_lag2),   
    !is.na(Y1),      !is.na(Y2),
    !is.na(dt),      dt > 0
  ) %>%
  mutate(
    log_dt    = log(dt),
    person_id = factor(person_id)
  )

# ============================================================================
# 14. DEFINE brms MULTIVARIATE VAR(1) MODEL 
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
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  backend = "cmdstanr",
  file = "./data/ar1.RDS"
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
  seed    = 16,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  backend = "cmdstanr",
  file = "./data/ar2.RDS"
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
    Y1_lag_wc + Y2_lag_wc +
    # Between-person
    ADHD_numeric + Y1_pm + Y2_pm +
    # ADHD moderates dynamics
    Y1_lag_wc:ADHD_numeric + Y2_lag_wc:ADHD_numeric +
    # Social context
    SocialMedia + InPerson + Alone +
    # ADHD × Context
    ADHD_x_SocialMedia + ADHD_x_InPerson + ADHD_x_Alone +
    # ── Weekend ──
    weekend + ADHD_x_weekend +
    # Internalizing
    C_SCARED_Score +
    # Extern
    DBD_Cp_Scr +
    # Age 
    Age +
    # Spacing correction
    log_dt + log_dt:ADHD_numeric +
    # Random effects
    (1 + Y2_lag_wc + Y1_lag_wc + log_dt | p | person_id)
) + gaussian()

bf_Y2 <- bf(
  Y2 ~ 1 +
    Y1_lag_wc + Y2_lag_wc +
    ADHD_numeric + Y1_pm + Y2_pm +
    Y1_lag_wc:ADHD_numeric + Y2_lag_wc:ADHD_numeric +
    SocialMedia + InPerson + Alone +
    ADHD_x_SocialMedia + ADHD_x_InPerson + ADHD_x_Alone +
    weekend  + ADHD_x_weekend +
    C_SCARED_Score +
    DBD_Cp_Scr +
    Age +
    log_dt + log_dt:ADHD_numeric +
    (1 + Y2_lag_wc + Y1_lag_wc + log_dt | p | person_id)
) + gaussian()

fit_brms <- brm(
  bf_Y1 + bf_Y2 + set_rescor(TRUE),
  data    = brms_dat,
  cores   = 4,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 16,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  file = "./data/rawBinaryVars.RDS"
)

# ============================================================================
# 15b. FIT THE MODEL -- Centered Social Shocks
# ============================================================================
bf_Y1 <- bf(
  Y1 ~ 1 +
    # Temporal dynamics
    Y1_lag_wc + Y2_lag_wc +
    alone_pm + inper_pm + socia_pm +
    # Between-person
    ADHD_numeric + Y1_pm + Y2_pm +
    # ADHD moderates dynamics
    Y1_lag_wc:ADHD_numeric + Y2_lag_wc:ADHD_numeric +
    # Social context
    socia_wc + inper_wc + alone_wc +
    # ADHD × Context
    ADHD_x_SocialMediawc + ADHD_x_InPersonwc + ADHD_x_Alonewc +
    # ── Weekend ──
    weekend + ADHD_x_weekend +
    # Internalizing
    C_SCARED_Score +
    # Extern
    DBD_Cp_Scr +
    # Age
    Age +
    # Spacing correction
    log_dt + log_dt:ADHD_numeric +
    log_dt:Y1_lag_wc + log_dt:Y2_lag_wc + 
    log_dt:ADHD_numeric:Y1_lag_wc + log_dt:ADHD_numeric:Y2_lag_wc +
    # Random effects
    (1 + Y2_lag_wc + Y1_lag_wc + log_dt | p | person_id)
) + gaussian()
bf_Y2 <- bf(
  Y2 ~ 1 +
    Y1_lag_wc + Y2_lag_wc +
    alone_pm + inper_pm + socia_pm +
    ADHD_numeric + Y1_pm + Y2_pm +
    Y1_lag_wc:ADHD_numeric + Y2_lag_wc:ADHD_numeric +
    socia_wc + inper_wc + alone_wc +
    ADHD_x_SocialMediawc + ADHD_x_InPersonwc + ADHD_x_Alonewc +
    weekend + ADHD_x_weekend +
    C_SCARED_Score +
    DBD_Cp_Scr +
    Age +
    log_dt + log_dt:ADHD_numeric +
    log_dt:Y1_lag_wc + log_dt:Y2_lag_wc + 
    log_dt:ADHD_numeric:Y1_lag_wc + log_dt:ADHD_numeric:Y2_lag_wc +
    (1 + Y2_lag_wc + Y1_lag_wc + log_dt | p | person_id)
) + gaussian()
fit_brms2 <- brm(
  bf_Y1 + bf_Y2 + set_rescor(TRUE),
  data = brms_dat,
  cores   = 4,
  chains  = 4,
  iter    = 4000,
  warmup  = 2000,
  seed    = 16,
  control = list(adapt_delta = 0.97, max_treedepth = 18),
  backend = "cmdstanr",
  file = "./data/scaledBinaryVars22-ARCRLagSpace2.RDS"
)

## Prepare the output table
sum_vals <- summary(fit_brms2)$fixed

# Flag rows where the CI excludes 0 (lower and upper have the same sign)
lower <- sum_vals[["l-95% CI"]]
upper <- sum_vals[["u-95% CI"]]
sig   <- lower > 0 | upper < 0

# Format helper (change digits as needed)
fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

# Build the output table
out_tab <- data.frame(
  Parameter = rownames(sum_vals),
  `Estimate {95% CI}` = paste0(
    fmt(sum_vals$Estimate), ifelse(sig, "*", ""),
    " [", fmt(lower), ", ", fmt(upper), "]"
  ),
  `SE` = fmt(sum_vals$Est.Error),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## NOw split so Y1 variables on the left and Y2 on the right
out_tab_Y1 = out_tab[grep(x = out_tab$Parameter, pattern = "^Y1_"),]
out_tab_Y2 = out_tab[grep(x = out_tab$Parameter, pattern = "^Y2_"),]
## NOw cbind these?
to_write = bind_cols(out_tab_Y1, out_tab_Y2)
write.table(x = to_write, file = "./data/outputPretty.csv", 
            quote = F, row.names = F, sep ="\t")


# ============================================================================
# 16. PLOT EMM FROM SOCIAL VARS
# ============================================================================
library(emmeans)

# Outcome Y1 (Hostile)
inper_y1_0     <- emmeans(fit_brms2, ~ inper_wc, at = list(inper_wc = c(-1, 0, 1), ADHD_x_InPersonwc = 0, ADHD_numeric = 0), resp = "Y1")
inper_y1_1_neg <- emmeans(fit_brms2, ~ inper_wc, at = list(inper_wc = -1, ADHD_x_InPersonwc = -1, ADHD_numeric = 1), resp = "Y1")
inper_y1_1_zer <- emmeans(fit_brms2, ~ inper_wc, at = list(inper_wc = 0,  ADHD_x_InPersonwc = 0,  ADHD_numeric = 1), resp = "Y1")
inper_y1_1_pos <- emmeans(fit_brms2, ~ inper_wc, at = list(inper_wc = 1,  ADHD_x_InPersonwc = 1,  ADHD_numeric = 1), resp = "Y1")

# Outcome Y2 (Prosocial)
inper_y2_0     <- emmeans(fit_brms2, ~ inper_wc, at = list(inper_wc = c(-1, 0, 1), ADHD_x_InPersonwc = 0, ADHD_numeric = 0), resp = "Y2")
inper_y2_1_neg <- emmeans(fit_brms2, ~ inper_wc, at = list(inper_wc = -1, ADHD_x_InPersonwc = -1, ADHD_numeric = 1), resp = "Y2")
inper_y2_1_zer <- emmeans(fit_brms2, ~ inper_wc, at = list(inper_wc = 0,  ADHD_x_InPersonwc = 0,  ADHD_numeric = 1), resp = "Y2")
inper_y2_1_pos <- emmeans(fit_brms2, ~ inper_wc, at = list(inper_wc = 1,  ADHD_x_InPersonwc = 1,  ADHD_numeric = 1), resp = "Y2")

# ==========================================
# 2. SOCIAL MEDIA SETTING (Values: -1, 0, 1)
# ==========================================
# Outcome Y1 (Hostile)
socia_y1_0     <- emmeans(fit_brms2, ~ socia_wc, at = list(socia_wc = c(-1, 0, 1), ADHD_x_SocialMediawc = 0, ADHD_numeric = 0), resp = "Y1")
socia_y1_1_neg <- emmeans(fit_brms2, ~ socia_wc, at = list(socia_wc = -1, ADHD_x_SocialMediawc = -1, ADHD_numeric = 1), resp = "Y1")
socia_y1_1_zer <- emmeans(fit_brms2, ~ socia_wc, at = list(socia_wc = 0,  ADHD_x_SocialMediawc = 0,  ADHD_numeric = 1), resp = "Y1")
socia_y1_1_pos <- emmeans(fit_brms2, ~ socia_wc, at = list(socia_wc = 1,  ADHD_x_SocialMediawc = 1,  ADHD_numeric = 1), resp = "Y1")

# Outcome Y2 (Prosocial)
socia_y2_0     <- emmeans(fit_brms2, ~ socia_wc, at = list(socia_wc = c(-1, 0, 1), ADHD_x_SocialMediawc = 0, ADHD_numeric = 0), resp = "Y2")
socia_y2_1_neg <- emmeans(fit_brms2, ~ socia_wc, at = list(socia_wc = -1, ADHD_x_SocialMediawc = -1, ADHD_numeric = 1), resp = "Y2")
socia_y2_1_zer <- emmeans(fit_brms2, ~ socia_wc, at = list(socia_wc = 0,  ADHD_x_SocialMediawc = 0,  ADHD_numeric = 1), resp = "Y2")
socia_y2_1_pos <- emmeans(fit_brms2, ~ socia_wc, at = list(socia_wc = 1,  ADHD_x_SocialMediawc = 1,  ADHD_numeric = 1), resp = "Y2")

# ==========================================
# 3. ALONE SETTING (Values: -1, 0, 1)
# ==========================================
# Outcome Y1 (Hostile)
alone_y1_0     <- emmeans(fit_brms2, ~ alone_wc, at = list(alone_wc = c(-1, 0, 1), ADHD_x_Alonewc = 0, ADHD_numeric = 0), resp = "Y1")
alone_y1_1_neg <- emmeans(fit_brms2, ~ alone_wc, at = list(alone_wc = -1, ADHD_x_Alonewc = -1, ADHD_numeric = 1), resp = "Y1")
alone_y1_1_zer <- emmeans(fit_brms2, ~ alone_wc, at = list(alone_wc = 0,  ADHD_x_Alonewc = 0,  ADHD_numeric = 1), resp = "Y1")
alone_y1_1_pos <- emmeans(fit_brms2, ~ alone_wc, at = list(alone_wc = 1,  ADHD_x_Alonewc = 1,  ADHD_numeric = 1), resp = "Y1")

# Outcome Y2 (Prosocial)
alone_y2_0     <- emmeans(fit_brms2, ~ alone_wc, at = list(alone_wc = c(-1, 0, 1), ADHD_x_Alonewc = 0, ADHD_numeric = 0), resp = "Y2")
alone_y2_1_neg <- emmeans(fit_brms2, ~ alone_wc, at = list(alone_wc = -1, ADHD_x_Alonewc = -1, ADHD_numeric = 1), resp = "Y2")
alone_y2_1_zer <- emmeans(fit_brms2, ~ alone_wc, at = list(alone_wc = 0,  ADHD_x_Alonewc = 0,  ADHD_numeric = 1), resp = "Y2")
alone_y2_1_pos <- emmeans(fit_brms2, ~ alone_wc, at = list(alone_wc = 1,  ADHD_x_Alonewc = 1,  ADHD_numeric = 1), resp = "Y2")

# ==========================================
# 4. COMBINE 
# ==========================================
df_inper <- bind_rows(
  as.data.frame(inper_y1_0)     %>% mutate(Outcome = "Hostile", Group = "Control"),
  as.data.frame(inper_y1_1_neg) %>% mutate(Outcome = "Hostile", Group = "ADHD"),
  as.data.frame(inper_y1_1_zer) %>% mutate(Outcome = "Hostile", Group = "ADHD"),
  as.data.frame(inper_y1_1_pos) %>% mutate(Outcome = "Hostile", Group = "ADHD"),
  as.data.frame(inper_y2_0)     %>% mutate(Outcome = "Prosocial", Group = "Control"),
  as.data.frame(inper_y2_1_neg) %>% mutate(Outcome = "Prosocial", Group = "ADHD"),
  as.data.frame(inper_y2_1_zer) %>% mutate(Outcome = "Prosocial", Group = "ADHD"),
  as.data.frame(inper_y2_1_pos) %>% mutate(Outcome = "Prosocial", Group = "ADHD")
) %>% rename(Setting_Val = 1) %>% mutate(Context = "In-Person")

df_socia <- bind_rows(
  as.data.frame(socia_y1_0)     %>% mutate(Outcome = "Hostile", Group = "Control"),
  as.data.frame(socia_y1_1_neg) %>% mutate(Outcome = "Hostile", Group = "ADHD"),
  as.data.frame(socia_y1_1_zer) %>% mutate(Outcome = "Hostile", Group = "ADHD"),
  as.data.frame(socia_y1_1_pos) %>% mutate(Outcome = "Hostile", Group = "ADHD"),
  as.data.frame(socia_y2_0)     %>% mutate(Outcome = "Prosocial", Group = "Control"),
  as.data.frame(socia_y2_1_neg) %>% mutate(Outcome = "Prosocial", Group = "ADHD"),
  as.data.frame(socia_y2_1_zer) %>% mutate(Outcome = "Prosocial", Group = "ADHD"),
  as.data.frame(socia_y2_1_pos) %>% mutate(Outcome = "Prosocial", Group = "ADHD")
) %>% rename(Setting_Val = 1) %>% mutate(Context = "Social Media")

df_alone <- bind_rows(
  as.data.frame(alone_y1_0)     %>% mutate(Outcome = "Hostile", Group = "Control"),
  as.data.frame(alone_y1_1_neg) %>% mutate(Outcome = "Hostile", Group = "ADHD"),
  as.data.frame(alone_y1_1_zer) %>% mutate(Outcome = "Hostile", Group = "ADHD"),
  as.data.frame(alone_y1_1_pos) %>% mutate(Outcome = "Hostile", Group = "ADHD"),
  as.data.frame(alone_y2_0)     %>% mutate(Outcome = "Prosocial", Group = "Control"),
  as.data.frame(alone_y2_1_neg) %>% mutate(Outcome = "Prosocial", Group = "ADHD"),
  as.data.frame(alone_y2_1_zer) %>% mutate(Outcome = "Prosocial", Group = "ADHD"),
  as.data.frame(alone_y2_1_pos) %>% mutate(Outcome = "Prosocial", Group = "ADHD")
) %>% rename(Setting_Val = 1) %>% mutate(Context = "Alone")

plot_data <- bind_rows(df_inper, df_socia, df_alone) %>%
  mutate(Setting_Val = factor(Setting_Val),
         Context = factor(Context, levels = c("Social Media", "In-Person", "Alone")))

# ==========================================
# 5. PLOT
# ==========================================
fig_3 = ggplot(
  plot_data, 
  aes(x = Setting_Val, y = emmean, fill = Group)
) +
  geom_col(
    position = position_dodge(width = 0.8), 
    width = 0.7, 
    color = "black", 
    linewidth = 0.3
  ) +
  geom_errorbar(
    aes(ymin = lower.HPD, ymax = upper.HPD), 
    position = position_dodge(width = 0.8), 
    width = 0.2,
    linewidth = 0.5
  ) +
  facet_grid(Outcome ~ Context, scales = "free_y") +
  labs(
    #title = "Estimated Marginal Means Across Social Settings & ADHD Status",
    #subtitle = "Comparing Levels (-1, 0, 1) Within Each Social Category for Hostile and Prosocial Outcomes",
    x = "Person-mean Social Setting Level",
    y = "Predicted Within-person Deviations",
    fill = "ADHD Status"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1.5, "lines")
  ) +
  scale_fill_manual(values = c("Control" = "#4575b4", "ADHD" = "#d73027"))

ggsave(filename = "figure_3.png", plot = fig_3,
       dpi = 300, height = 8, width = 8)



# ============================================================================
# 17. MODEL CONVERGENCE
# ============================================================================
rhat_vals <- brms::rhat(fit_brms2)

# 2. Convert to a data frame for ggplot2
rhat_df <- data.frame(rhat = rhat_vals)

# 3. Plot a histogram of the R-hat values
plot = ggplot(rhat_df, aes(x = rhat)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 1.05, color = "red", linetype = "dashed") +
  labs(
    title = "Histogram of R-hat Convergence Diagnostics",
    x = "R-hat Value",
    y = "Count"
  ) +
  theme_minimal()

ggsave(filename = "./data/out_rhat_example.png", plot = plot, 
       width = 7, height = 4, dpi = 300)

plot = mcmc_trace(fit_brms2, pars = vars(starts_with("b_Y1_Y1_lag_wc")))
ggsave(filename = "./data/out_trace_example.png", plot = plot, 
       width = 7, height = 4, dpi = 300)

plot = mcmc_acf(fit_brms2, pars = vars(starts_with("b_Y1_Y1_lag_wc")))
ggsave(filename = "./data/out_acf_example.png", plot = plot, 
       width = 7, height = 4, dpi = 300)

# 1. Compute both Bulk and Tail ESS using bayestestR
ess_summary <- bayestestR::effective_sample(fit_brms2)

# 2. Reshape 
ess_df <- ess_summary %>%
  pivot_longer(
    cols = c(ESS, ESS_tail), 
    names_to = "ESS_Type", 
    values_to = "ESS_Value"
  ) %>%
  mutate(ESS_Type = gsub("ESS_", "", ESS_Type))

# 3. Plot facetted histograms
plot = ggplot(ess_df, aes(x = ESS_Value, fill = ESS_Type)) +
  geom_histogram(bins = 30, color = "white", show.legend = FALSE) +
  facet_wrap(~ ESS_Type, scales = "free_x") +
  scale_fill_manual(values = c("Bulk" = "steelblue", "Tail" = "darkorange")) +
  geom_vline(xintercept = 400, color = "red", linetype = "dashed") +
  labs(
    title = "Histograms of Bulk and Tail Effective Sample Size (ESS)",
    x = "Effective Sample Size Value",
    y = "Count"
  ) +
  theme_minimal()

ggsave(filename = "./data/out_ess_example.png", plot = plot, 
       width = 7, height = 4, dpi = 300)



# ============================================================================
# 17. RANDOM EFFECT VCOV MATRIX
# ============================================================================

# 1. EXTRACT VCOV
vc <- VarCorr(fit_brms2, summary = TRUE)$person_id

# Correct indexing: [row, statistic, col]
est_mat <- vc$cor[, "Estimate", ]
lo_mat  <- vc$cor[, "Q2.5",     ]
hi_mat  <- vc$cor[, "Q97.5",    ]

cat("Correlation matrix (posterior means):\n")
print(round(est_mat, 2))

# 2. ORGANIZE
param_order <- c(
  "Y1_Intercept",
  "Y1_Y1_lag_wc",
  "Y1_Y2_lag_wc",
  "Y1_log_dt",
  "Y2_Intercept",
  "Y2_Y1_lag_wc",
  "Y2_Y2_lag_wc",
  "Y2_log_dt"
)

clean_labels <- c(
  "Y1_Intercept"  = "HOS\nIntercept",
  "Y1_Y1_lag_wc"  = "HOS\nAR",
  "Y1_Y2_lag_wc"  = "HOS←PS\nCL",
  "Y1_log_dt"     = "HOS\nlog Δt",
  "Y2_Intercept"  = "PS\nIntercept",
  "Y2_Y1_lag_wc"  = "PS←HOS\nCL",
  "Y2_Y2_lag_wc"  = "PS\nAR",
  "Y2_log_dt"     = "PS\nlog Δt"
)

n <- length(param_order)

# Reorder matrices to match param_order
est_mat <- est_mat[param_order, param_order]
lo_mat  <- lo_mat[param_order,  param_order]
hi_mat  <- hi_mat[param_order,  param_order]

# EXTRACT SD
sd_vec  <- vc$sd[param_order, "Estimate"]
var_vec <- sd_vec

cat("\nVariances (SD²):\n")
print(round(var_vec, 4))

# ORGANIZE MAT
display_mat <- est_mat
diag(display_mat) <- var_vec

# Credibility: 95% BCI excludes zero (off-diagonal only)
cred_mat <- (lo_mat > 0 | hi_mat < 0)
diag(cred_mat) <- FALSE

# RESHAPE
label_vec <- clean_labels[param_order]

mat_long <- expand.grid(
  row_idx = seq_len(n),
  col_idx = seq_len(n)
) %>%
  mutate(
    row_lab  = factor(label_vec[row_idx], levels = rev(label_vec)),
    col_lab  = factor(label_vec[col_idx], levels = label_vec),
    value    = mapply(function(r, c) display_mat[r, c], row_idx, col_idx),
    is_diag  = row_idx == col_idx,
    credible = mapply(function(r, c) cred_mat[r, c], row_idx, col_idx),
    text_val = case_when(
      is_diag ~ sprintf("%.2f", value),
      TRUE    ~ sprintf("%.2f",     value)
    ),
    text_col = case_when(
      is_diag          ~ "grey20",
      abs(value) > 0.5 ~ "white",
      TRUE             ~ "grey15"
    ),
    text_face = ifelse(credible, "bold", "plain")
  )

# PLOT
p <- ggplot(mat_long, aes(x = col_lab, y = row_lab)) +
  
  # Off-diagonal: correlation fill
  geom_tile(
    data  = filter(mat_long, !is_diag),
    aes(fill = value),
    color = "white", linewidth = 0.7
  ) +
  geom_tile(
    data  = filter(mat_long, is_diag),
    aes(fill = value),
    #fill  = "#DCDCDC",
    color = "white",
    linewidth = 0.7
  ) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "#F7F7F7",
    high     = "#D6604D",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Posterior mean correlation",
    guide    = guide_colorbar(
      barwidth       = 10,
      barheight      = 0.8,
      title.position = "top",
      title.hjust    = 0.5
    )
  ) +
  
  geom_text(
    aes(label = text_val, colour = text_col),
    size       = 6.0,
    lineheight = 0.85,
    fontface   = mat_long$text_face
  ) +
  scale_colour_identity() +
  
  # Credibility star
  geom_text(
    data    = filter(mat_long, credible),
    aes(label = "✶"),
    colour  = "black",
    size    = 4.5,
    nudge_x =  0.37,
    nudge_y =  0.33
  ) +
  
  # Block dividers between Y1 and Y2 groups
  geom_hline(yintercept = 4.5, linewidth = 1.1, colour = "grey25") +
  geom_vline(xintercept = 4.5, linewidth = 1.1, colour = "grey25") +
  
  # Equation block labels
  annotate("text", x = 2.5, y = n + 0.9,
           label    = "Hostile Equation",
           fontface = "bold", size = 3.6, colour = "grey20") +
  annotate("text", x = 6.5, y = n + 0.9,
           label    = "Prosocial Equation",
           fontface = "bold", size = 3.6, colour = "grey20") +
  
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.x     = element_text(size = 12,  colour = "grey20",
                                   lineheight = 0.85),
    axis.text.y     = element_text(size = 12,  colour = "grey20",
                                   lineheight = 0.85),
    panel.grid      = element_blank(),
    legend.position = "bottom",
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 9, colour = "grey45"),
    plot.margin     = margin(28, 10, 10, 10)
  )

print(p)

ggsave("plot_random_effect_cormat_v3.png", p,
       width = 8, height = 8.5, dpi = 200, bg = "white")
