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

# ============================================================================
# 8. FILTER: MINIMUM OBSERVATIONS + VARIANCE CHECK
# ============================================================================

sufficient_persons <- ema_data %>%
  count(ID) %>%
  filter(n >= 3) %>%
  pull(ID)

ema_data <- ema_data %>% filter(ID %in% sufficient_persons)

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

session_summary <- ema_data %>%
  group_by(ID, session) %>%
  summarise(n_obs = n(), .groups = "drop")

# ============================================================================
# 9. COMPUTE MEAN
# ============================================================================
ema_data <- ema_data %>%
  mutate(
    anti_mean = rowMeans(cbind(item1, item2, item3, item4), na.rm = TRUE),
    pro_mean  = rowMeans(cbind(item5, item6, item7, item8), na.rm = TRUE),
    anti_mean = ifelse(is.nan(anti_mean), NA_real_, anti_mean),
    pro_mean  = ifelse(is.nan(pro_mean),  NA_real_, pro_mean)
  )

# Grand-mean 
anti_gm        <- mean(ema_data$anti_mean, na.rm = TRUE)
pro_gm         <- mean(ema_data$pro_mean,  na.rm = TRUE)
anti_sd_global <- sd(ema_data$anti_mean,   na.rm = TRUE)
pro_sd_global  <- sd(ema_data$pro_mean,    na.rm = TRUE)

ema_data <- ema_data %>%
  mutate(
    Y1 = (anti_mean - anti_gm) / anti_sd_global,
    Y2 = (pro_mean  - pro_gm)  / pro_sd_global
  )

# ============================================================================
# 9b. COMPUTE WITHIN-PERSON CENTRED VARIABLES
# ============================================================================

ema_data <- ema_data %>%
  group_by(ID) %>%
  mutate(
    Y1_wpc = Y1 - mean(Y1, na.rm = TRUE),   # remove person mean from GMC score
    Y2_wpc = Y2 - mean(Y2, na.rm = TRUE)
  ) %>%
  ungroup()

# ============================================================================
# 10. PROCESS BINARY SOCIAL CONTEXT INDICATORS
# ============================================================================
if (!"ADHD_numeric" %in% names(ema_data)) {
  ema_data <- ema_data %>%
    left_join(time_invar %>% select(ID, ADHD_numeric), by = "ID")
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

P <- nrow(id_mapping)
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

violations <- within_check %>% filter(max_dt > SESSION_GAP_HOURS)
if (nrow(violations) > 0) {
  cat(sprintf("  WARNING: %d sessions with gap > %d hours\n",
              nrow(violations), SESSION_GAP_HOURS))
}
# ============================================================================
#  FIGURE: Individual time series
# ============================================================================
# ============================================================================
# 1. BASIC DATA STRUCTURE
# ============================================================================
person_summary <- ema_data %>%
  group_by(ID) %>%
  summarise(
    n_obs       = n(),
    n_sessions  = n_distinct(session),
    total_hours = max(RawHours) - min(RawHours),
    ADHD        = first(ADHD_numeric),
    anti_mean   = mean(Y1, na.rm = TRUE),
    anti_sd     = sd(Y1, na.rm = TRUE),
    pro_mean    = mean(Y2, na.rm = TRUE),
    pro_sd      = sd(Y2, na.rm = TRUE),
    anti_range  = max(Y1, na.rm = TRUE) - min(Y1, na.rm = TRUE),
    pro_range   = max(Y2, na.rm = TRUE) - min(Y2, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(Group = ifelse(ADHD == 1, "ADHD", "Control"))

# ============================================================================
# 2. OVERALL DISTRIBUTIONS (GMC and WPC side by side)
# ============================================================================

for (centering in c("GMC", "WPC")) {
  v1 <- if (centering == "GMC") "Y1" else "Y1_wpc"
  v2 <- if (centering == "GMC") "Y2" else "Y2_wpc"
  label <- if (centering == "GMC") "Grand-mean centred" else "Within-person centred"
}

dist_long <- bind_rows(
  ema_data %>%
    select(Y1, Y2) %>%
    pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
    mutate(
      Centering = "Grand-mean centred",
      Variable  = recode(Variable, Y1 = "Hostile", Y2 = "Prosocial")
    ),
  ema_data %>%
    select(Y1_wpc, Y2_wpc) %>%
    pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
    mutate(
      Centering = "Within-person centred",
      Variable  = recode(Variable, Y1_wpc = "Hostile", Y2_wpc = "Prosocial")
    )
) %>%
  mutate(Centering = factor(Centering,
                            levels = c("Grand-mean centred",
                                       "Within-person centred")))

p_hist <- ggplot(dist_long, aes(x = Value, fill = Centering)) +
  geom_histogram(bins = 30, colour = "white", alpha = 0.85) +
  scale_fill_manual(values = c("Grand-mean centred"    = "#2166ac",
                               "Within-person centred" = "#d6604d")) +
  facet_grid(Centering ~ Variable, scales = "free") +
  labs(
    x = "Average Endorsment", y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    strip.text    = element_text(face = "bold"),
    legend.position = "none"
  )

# ============================================================================
# 3. TIME INTERVAL DISTRIBUTION
# ============================================================================

intervals <- ema_data %>%
  arrange(ctsem_id, ctsem_time) %>%
  group_by(ctsem_id) %>%
  mutate(dt = ctsem_time - lag(ctsem_time)) %>%
  filter(!is.na(dt)) %>%
  ungroup()


p_dt <- ggplot(intervals, aes(x = dt)) +
  geom_histogram(bins = 50, fill = "black", colour = "white", alpha = 0.8) + #, fill = "#e66101"
  scale_x_continuous(breaks = seq(0, 30, by = 2)) +
  labs(x = "Δt (hours)", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))


# ============================================================================
# 4. INDIVIDUAL SPAGHETTI PLOTS — GMC and WPC
# ============================================================================
make_spaghetti <- function(data, y1_var, y2_var, centering_label) {
  data %>%
    mutate(Group = ifelse(ADHD_numeric == 1, "ADHD", "Control")) %>%
    select(ID, Group, ctsem_id, ctsem_time,
           Y1 = all_of(y1_var), Y2 = all_of(y2_var)) %>%
    pivot_longer(c(Y1, Y2), names_to = "Process", values_to = "Value") %>%
    mutate(Process = recode(Process, Y1 = "Hostile", Y2 = "Prosocial")) %>%
    ggplot(aes(x     = ctsem_time,
               y     = Value,
               group = interaction(ID, ctsem_id),
               colour = Group)) +
    geom_line(alpha = 0.15, linewidth = 0.3) +
    geom_smooth(aes(group = Group), method = "loess", se = TRUE,
                linewidth = 1.2, alpha = 0.3) +
    scale_colour_manual(values = c("Control" = "#4575b4", "ADHD" = "#d73027")) +
    facet_wrap(~ Process) +
    labs(
      x = "Time within session (hours)", y = centering_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom"
    )
}

p_spaghetti_gmc <- make_spaghetti(ema_data, "Y1",     "Y2",     "Grand-mean centred")
p_spaghetti_wpc <- make_spaghetti(ema_data, "Y1_wpc", "Y2_wpc", "Within-person centred")

library(ggpubr)

fig_2 <- ggarrange(
  p_hist,
  p_dt,
  p_spaghetti_gmc,
  p_spaghetti_wpc,
  ncol   = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D")
)

ggsave(filename = "figure_2.png", plot = fig_2,
       dpi = 300, height = 8, width = 8)
