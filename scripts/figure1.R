library(readxl)
library(tidyverse)
library(lubridate)
library(viridis)
library(lme4)
library(brms)
library(ggpubr)

# ============================================================================
# LOAD AND PROCESS RAW DATA
# ============================================================================

in_dat <- read_xlsx("./data/SPIG EMA updated data.xlsx")
in_dat$Response[which(in_dat$Response == "#skipped#")] <- NA

# Rename columns immediately to avoid backtick issues
names(in_dat) <- make.names(names(in_dat))

# ============================================================================
# PROCESS DATA BY PARTICIPANT
# ============================================================================

## Isolate the questions of interest
question_vars <- c(
  "adhdig_ema_1", "adhdig_ema_2", "adhdig_ema_3", "adhdig_ema_4",
  "adhdig_ema_5", "adhdig_ema_6", "adhdig_ema_7", "adhdig_ema_8",
  "adhdig_ema_9"
)

## Process each individual
all_data2 <- NULL
n_processed <- 0

for(i in unique(in_dat$Participant.ID)){
  n_processed <- n_processed + 1
  if (n_processed %% 10 == 0) {
    cat(sprintf("  Processed %d participants...\n", n_processed))
  }
  
  iso_dat <- in_dat[in_dat$Participant.ID == i, ]
  iso_dat$ID <- i
  
  ## Convert from long to wide
  wide_iso_dat <- tidyr::pivot_wider(
    iso_dat,
    id_cols = c(Notification.Time, Session.Instance.No, ADHD),
    names_from = Prompt.Label,
    values_from = Response
  )
  
  wide_iso_dat$ID <- i
  
  ## Select columns
  wide_iso_dat <- wide_iso_dat %>%
    distinct() %>%
    select(Notification.Time, Session.Instance.No, ID, ADHD,
           adhdig_ema_1, adhdig_ema_2, adhdig_ema_3, adhdig_ema_4,
           adhdig_ema_5, adhdig_ema_6, adhdig_ema_7, adhdig_ema_8,
           `adhdig_ema_9_I am alone`,
           `adhdig_ema_9_Through social media`,
           `adhdig_ema_9_In-person`)
  
  all_data2 <- dplyr::bind_rows(all_data2, wide_iso_dat)
}

all_data <- all_data2
rm(all_data2)

# Remove problematic participant observations
all_data <- all_data[!all_data$ID %in% c(76826), ]

## First remove all NA observation counts
all_data = all_data[-which(is.na(all_data$Session.Instance.No)),]


# ============================================================================
# RENAME COLUMNS FOR EASIER HANDLING
# ============================================================================

colnames(all_data) <- c("time_col", "SessionNumb", "id", "ADHD",
                        "item1", "item2", "item3", "item4",
                        "item5", "item6", "item7", "item8",
                        "Alone", "SocialMedia", "InPerson")

# ============================================================================
# CALCULATE TIME INTERVALS AND PROCESS
# ============================================================================

ema_data <- all_data %>%
  rename(ID = id) %>%
  group_by(ID) %>%
  arrange(ID, time_col) %>%
  filter(complete.cases(item1, item2, item3, item4, item5, item6, item7, item8)) %>%
  mutate(
    TimeHours = as.numeric(difftime(time_col, first(time_col), units = "hours")),
    time_interval = TimeHours - lag(TimeHours, default = 0),
    
    # Convert to numeric
    item1 = as.numeric(item1),
    item2 = as.numeric(item2),
    item3 = as.numeric(item3),
    item4 = as.numeric(item4),
    item5 = as.numeric(item5),
    item6 = as.numeric(item6),
    item7 = as.numeric(item7),
    item8 = as.numeric(item8),
    Alone = as.numeric(Alone),
    SocialMedia = as.numeric(SocialMedia),
    InPerson = as.numeric(InPerson)
  ) %>%
  ungroup() %>%
  select(ID, time_col, TimeHours, time_interval,
         item1, item2, item3, item4, item5, item6, item7, item8,
         Alone, SocialMedia, InPerson, ADHD, SessionNumb)

## Now calculate reliability following: Nezlek 2016
rel_dat <- ema_data %>%
  group_by(ID) %>%
  mutate(occasion = row_number()) %>%
  ungroup() %>%
  pivot_longer(
    cols = starts_with("item"),
    names_to = "item",
    values_to = "value"
  ) %>%
  select(ID, occasion, item, value)
## Now estimate the three level nested model with only intercepts
mod_rel = lmer(value ~ 1 + (1|ID) + (1|occasion) + (1|item), data = rel_dat)
item_rel = VarCorr(mod_rel)$occasion[1] / (VarCorr(mod_rel)$occasion[1] + (VarCorr(mod_rel)$item[1] / 8))
print(item_rel)

## Now do within scale
mod_rel1 = lmer(value ~ 1 + (1|ID) + (1|occasion) + (1|item), data = rel_dat[which(rel_dat$item %in% paste("item", 1:4, sep='')),])
item_rel1 = VarCorr(mod_rel1)$occasion[1] / (VarCorr(mod_rel1)$occasion[1] + (VarCorr(mod_rel1)$item[1] / 4))
print(item_rel1)
mod_rel2 = lmer(value ~ 1 + (1|ID) + (1|occasion) + (1|item), data = rel_dat[which(!rel_dat$item %in% paste("item", 1:4, sep='')),])
item_rel2 = VarCorr(mod_rel2)$occasion[1] / (VarCorr(mod_rel2)$occasion[1] + (VarCorr(mod_rel2)$item[1] / 4))
print(item_rel2)


# ============================================================================
# CALCULATE ITEM-LEVEL VARIANCE FOR EACH PARTICIPANT
# ============================================================================

person_item_variance <- ema_data %>%
  select(ID, ADHD, item1:item8) %>%
  group_by(ID, ADHD) %>%
  summarise(
    n_obs = n(),
    
    # Variance for each item
    var_item1 = var(item1, na.rm = TRUE),
    var_item2 = var(item2, na.rm = TRUE),
    var_item3 = var(item3, na.rm = TRUE),
    var_item4 = var(item4, na.rm = TRUE),
    var_item5 = var(item5, na.rm = TRUE),
    var_item6 = var(item6, na.rm = TRUE),
    var_item7 = var(item7, na.rm = TRUE),
    var_item8 = var(item8, na.rm = TRUE),
    
    # SD for each item
    sd_item1 = sd(item1, na.rm = TRUE),
    sd_item2 = sd(item2, na.rm = TRUE),
    sd_item3 = sd(item3, na.rm = TRUE),
    sd_item4 = sd(item4, na.rm = TRUE),
    sd_item5 = sd(item5, na.rm = TRUE),
    sd_item6 = sd(item6, na.rm = TRUE),
    sd_item7 = sd(item7, na.rm = TRUE),
    sd_item8 = sd(item8, na.rm = TRUE),
    
    # Means for each item
    mean_item1 = mean(item1, na.rm = TRUE),
    mean_item2 = mean(item2, na.rm = TRUE),
    mean_item3 = mean(item3, na.rm = TRUE),
    mean_item4 = mean(item4, na.rm = TRUE),
    mean_item5 = mean(item5, na.rm = TRUE),
    mean_item6 = mean(item6, na.rm = TRUE),
    mean_item7 = mean(item7, na.rm = TRUE),
    mean_item8 = mean(item8, na.rm = TRUE),
    
    # Overall statistics
    overall_var = var(c_across(item1:item8), na.rm = TRUE),
    overall_sd = sd(c_across(item1:item8), na.rm = TRUE),
    overall_mean = mean(c_across(item1:item8), na.rm = TRUE),
    
    # Antisocial and Prosocial
    antisocial_mean = mean(c_across(item1:item4), na.rm = TRUE),
    antisocial_sd = sd(c_across(item1:item4), na.rm = TRUE),
    prosocial_mean = mean(c_across(item5:item8), na.rm = TRUE),
    prosocial_sd = sd(c_across(item5:item8), na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    ADHD_label = ifelse(ADHD == 1, "ADHD", "Control"),
    
    # Count items with zero or near-zero variance
    n_zero_var = rowSums(select(., starts_with("var_item")) == 0 | 
                           is.na(select(., starts_with("var_item"))), na.rm = TRUE),
    n_low_var = rowSums(select(., starts_with("var_item")) < 0.1, na.rm = TRUE),
    
    # Count items with very low SD
    n_zero_sd = rowSums(select(., starts_with("sd_item")) == 0 | 
                          is.na(select(., starts_with("sd_item"))), na.rm = TRUE),
    n_low_sd = rowSums(select(., starts_with("sd_item")) < 0.3, na.rm = TRUE),
    
    # Flags
    any_zero_var = n_zero_var > 0,
    multiple_low_var = n_low_var >= 3,
    overall_low_var = overall_sd < 0.5
  )

# ============================================================================
# COMPARE RESPONSE RATES
# ============================================================================
# source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
# rr_data = for_mod[!duplicated(for_mod$ID),]
# t.test(n_obs ~ ADHD, data = rr_data)
# summarySE(data = rr_data, measurevar = "n_obs", groupvars = "ADHD")

# ============================================================================
# EXAMINE PREDICTORS OF VARAINCE
# ============================================================================
library(visreg)
# RESHAPE FOR LMER MODEL
for_mod = reshape2::melt(person_item_variance, id.vars = c("ID", "ADHD", "n_obs"))
for_mod = for_mod[which(for_mod$variable %in% c("var_item1", "var_item2", "var_item3", "var_item4",
                                                "var_item5", "var_item6", "var_item7", "var_item8")),]
for_mod$value = as.numeric(for_mod$value)
mod = lmer(value ~ (ADHD + n_obs + variable)^3 + (1 | ID), data = for_mod)
mod_b = brm(value ~ (ADHD + n_obs + variable)^3 + (1 | ID), data = for_mod, warmup = 2000, iter = 4000)
car::Anova(mod)
visreg(mod, "ADHD", "variable", gg=TRUE) + theme_minimal()
conditional_effects(mod_b, effects = "ADHD:variable")

## Now do means
for_mod = reshape2::melt(person_item_variance, id.vars = c("ID", "ADHD", "n_obs"))
for_mod = for_mod[which(for_mod$variable %in% c("mean_item1", "mean_item2", "mean_item3", "mean_item4",
                                                "mean_item5", "mean_item6", "mean_item7", "mean_item8")),]
for_mod$value = as.numeric(for_mod$value)
mod = lmer(value ~ (ADHD + n_obs + variable)^3 + (1 | ID), data = for_mod)
mod_b = brm(value ~ (ADHD + n_obs + variable)^3 + (1 | ID), data = for_mod, warmup = 2000, iter = 4000)
car::Anova(mod)

# ============================================================================
# FIGURE: GROUP DIFFERENCES IN ITEM-LEVEL MEAN AND VARIANCE
# ============================================================================

library(tidyverse)
library(patchwork)
library(ggpubr)

# ============================================================================
# PREPARE PLOT DATA
# ============================================================================

# --- Variance data ---
var_plot_dat <- person_item_variance %>%
  select(ID, ADHD, starts_with("var_item")) %>%
  reshape2::melt(id.vars = c("ID", "ADHD")) %>%
  mutate(
    value = as.numeric(value),
    ADHD_label = factor(ifelse(ADHD == "Y", "ADHD", "Control"), 
                        levels = c("Control", "ADHD")),
    Item = factor(variable,
                  levels = paste0("var_item", 1:8),
                  labels = paste0("Item ", 1:8)),
    Subscale = ifelse(variable %in% c("var_item1", "var_item2", 
                                      "var_item3", "var_item4"),
                      "Hostile", "Prosocial")
  )

# --- Mean data ---
mean_plot_dat <- person_item_variance %>%
  select(ID, ADHD, starts_with("mean_item")) %>%
  reshape2::melt(id.vars = c("ID", "ADHD")) %>%
  mutate(
    value = as.numeric(value),
    ADHD_label = factor(ifelse(ADHD == "Y", "ADHD", "Control"), 
                        levels = c("Control", "ADHD")),
    Item = factor(variable,
                  levels = paste0("mean_item", 1:8),
                  labels = paste0("Item ", 1:8)),
    Subscale = ifelse(variable %in% c("mean_item1", "mean_item2", 
                                      "mean_item3", "mean_item4"),
                      "Hostile", "Prosocial")
  )

# ============================================================================
# SUMMARY STATISTICS FOR PLOTTING
# ============================================================================

var_summary <- var_plot_dat %>%
  group_by(ADHD_label, Item, Subscale) %>%
  summarise(
    mean_val  = mean(value, na.rm = TRUE),
    se_val    = sd(value, na.rm = TRUE) / sqrt(n()),
    lo        = mean_val - 1.96 * se_val,
    hi        = mean_val + 1.96 * se_val,
    .groups   = "drop"
  )

mean_summary <- mean_plot_dat %>%
  group_by(ADHD_label, Item, Subscale) %>%
  summarise(
    mean_val  = mean(value, na.rm = TRUE),
    se_val    = sd(value, na.rm = TRUE) / sqrt(n()),
    lo        = mean_val - 1.96 * se_val,
    hi        = mean_val + 1.96 * se_val,
    .groups   = "drop"
  )

# ============================================================================
# SHARED THEME AND PALETTE
# ============================================================================

adhd_colors <- c("Control" = "#3B82C4", "ADHD" = "#E05C3A")

shared_theme <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text         = element_text(face = "bold", size = 12),
    strip.background   = element_rect(fill = NULL, color = NA),
    legend.position    = "bottom",
    legend.title       = element_blank(),
    axis.title         = element_text(size = 12),
    axis.text.x        = element_text(size = 10, angle = 30, hjust = 1),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 11, color = "grey40")
  )

# ============================================================================
# PANEL A: ITEM-LEVEL MEAN RESPONSES
# ============================================================================

p_mean <- ggplot(mean_summary, 
                 aes(x = Item, y = mean_val, 
                     color = ADHD_label, group = ADHD_label)) +
  
  # Raw individual points (jittered, semi-transparent)
  geom_jitter(data = mean_plot_dat,
              aes(x = Item, y = value, color = ADHD_label),
              width = 0.15, alpha = 0.08, size = 1.2,
              inherit.aes = FALSE) +
  
  # Connecting lines
  geom_line(linewidth = 0.8, position = position_dodge(0.3)) +
  
  # Error bars (95% CI)
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = 0.15, linewidth = 0.7,
                position = position_dodge(0.3)) +
  
  # Points on top
  geom_point(size = 3.5, position = position_dodge(0.3)) +
  
  # Subscale facets
  facet_wrap(~ Subscale, scales = "free_x") +
  
  scale_color_manual(values = adhd_colors) +
  scale_y_continuous(name = "Mean Item Response", 
                     limits = c(1, NA), 
                     expand = expansion(mult = c(0.02, 0.1))) +
  labs(
    #title    = "A. Item-Level Mean Responses",
    #subtitle = "Points represent group means ± 95% CI; faded points are individual participants",
    x        = NULL
  ) +
  theme_minimal(base_size = 12) +
  shared_theme

# ============================================================================
# PANEL B: ITEM-LEVEL INTRAINDIVIDUAL VARIANCE
# ============================================================================

p_var <- ggplot(var_summary,
                aes(x = Item, y = mean_val,
                    color = ADHD_label, group = ADHD_label)) +
  
  # Raw individual points (jittered, semi-transparent)
  geom_jitter(data = var_plot_dat,
              aes(x = Item, y = value, color = ADHD_label),
              width = 0.15, alpha = 0.08, size = 1.2,
              inherit.aes = FALSE) +
  
  # Connecting lines
  geom_line(linewidth = 0.8, position = position_dodge(0.3)) +
  
  # Error bars (95% CI)
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = 0.15, linewidth = 0.7,
                position = position_dodge(0.3)) +
  
  # Points on top
  geom_point(size = 3.5, position = position_dodge(0.3)) +
  
  # Subscale facets
  facet_wrap(~ Subscale, scales = "free_x") +
  
  scale_color_manual(values = adhd_colors) +
  scale_y_continuous(name = "Mean Item Variance",
                     limits = c(0, NA),
                     expand = expansion(mult = c(0.02, 0.1))) +
  labs(
    #title    = "B. Item-Level Intraindividual Variance",
    #subtitle = "Points represent group means ± 95% CI; faded points are individual participants",
    x        = NULL
  ) +
  theme_minimal(base_size = 12) +
  shared_theme

# ============================================================================
# COMBINE AND SAVE
# ============================================================================

combined_figure <- p_mean / p_var +
  plot_annotation(
    #title   = "Group Differences in EMA Item Responses: ADHD vs. Control",
    #caption = "Note. Error bars reflect 95% confidence intervals around group means.\nAntisocial items = Items 1–4; Prosocial items = Items 5–8.",
    theme   = theme(
      plot.title   = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.caption = element_text(size = 9, color = "grey40", hjust = 0)
    )
  )

print(combined_figure)

combined_figure = ggpubr::ggarrange(p_mean, p_var, labels = "AUTO",ncol = 2, legend = "bottom", common.legend = TRUE)

ggsave("./Figure1_ItemLevel_Mean_Variance.png",
       plot   = combined_figure,
       width  = 10,
       height = 8,
       dpi    = 300,
       bg     = "white")
