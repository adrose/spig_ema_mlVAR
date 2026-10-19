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
  mutate(
    Notification.TimeOrig = Notification.Time,
    time_of_day = format(as.POSIXct(all_data$Notification.Time), "%H:%M:%S"),
    Notification.Time = as.numeric(Notification.Time))

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
#mod_rep1 = brm(responded ~ item + RawHours * ADHD + (1 + RawHours|ID), data = rep_diff_long, cores = 4)
#mod_rep2 = brm(responded ~ item + s(RawHours) + ADHD + (1 + RawHours|ID), data = rep_diff_long)
mod_rep1P = lme4::glmer(responded ~ item + RawHours * ADHD + (1 + RawHours|ID), data = rep_diff_long)
#mod_rep1P = mgcv::gamm(responded ~ item + s(RawHours) + s(RawHours, by = ADHD), data = rep_diff_long, family="binomial")

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
all_endorse = ema_data %>%
  count(Alone, SocialMedia, InPerson, name = "n") %>%
  arrange(desc(n)) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  print()

ema_data %>%
  group_by(ADHD) |> 
  count(Alone, SocialMedia, InPerson, name = "n") %>%
  arrange(ADHD, desc(n)) %>%
  mutate(pct = round(100 * n / sum(n), 1))

plot_data <- ema_data %>%
  mutate(
    Endorsement = case_when(
      Alone == 0 & SocialMedia == 0 & InPerson == 0 ~ "None",
      Alone == 1 & SocialMedia == 0 & InPerson == 0 ~ "Alone only",
      Alone == 0 & SocialMedia == 1 & InPerson == 0 ~ "Social media only",
      Alone == 0 & SocialMedia == 0 & InPerson == 1 ~ "In person only",
      Alone == 1 & SocialMedia == 1 & InPerson == 0 ~ "Alone + social media",
      Alone == 1 & SocialMedia == 0 & InPerson == 1 ~ "Alone + in person",
      Alone == 0 & SocialMedia == 1 & InPerson == 1 ~ "Social media + in person",
      Alone == 1 & SocialMedia == 1 & InPerson == 1 ~ "All three",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Endorsement)) %>%
  count(ID, Endorsement, name = "n") %>%
  complete(
    ID,
    Endorsement = c(
      "None",
      "Alone only",
      "Social media only",
      "In person only",
      "Alone + social media",
      "Alone + in person",
      "Social media + in person",
      "All three"
    ),
    fill = list(n = 0)
  ) %>%
  group_by(ID) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

ggplot(plot_data, aes(x = Endorsement, y = pct, fill = Endorsement)) +
  #geom_violin(trim = FALSE, show.legend = FALSE) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white"
  ) +
  geom_jitter() +
  labs(
    x = NULL,
    y = "Percentage of observations within participant",
    title = "Distribution of EMA Endorsement Combinations"
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none"
  )

##Now try a glmer model here
## First isolate the data
iso_vals = ema_data |> 
  dplyr::select(ID, ADHD, Alone, SocialMedia, InPerson, day_of_week, time_of_day) |> 
  pivot_longer(
    cols = c(Alone, SocialMedia, InPerson),
    names_to = "endorse",
    values_to = "value"
  )
## Convert time to numeric
iso_vals$time_period = hms(iso_vals$time_of_day)
iso_vals$hours_numeric <- scale(as.numeric(iso_vals$time_period, "hours"))[,]
## Now run spline on time and interactions on other variables
mod_est = brm(value ~ endorse + day_of_week + ADHD + s(hours_numeric) + (1 + endorse | ID), 
              data = iso_vals, family = bernoulli(link = "logit"))
