############################################################
# PB behavioural response re-analysis (ESR revision)
# Unified modelling workflow (binary primary model + optional type model)
# Data file: PBResponse_Data.csv OR PBResponse_Data.xlsx (in working dir)
############################################################


# ---- Packages ----

library(tidyverse)
library(janitor)
library(splines)
library(MuMIn)
library(performance)
library(DHARMa)
library(pROC)
library(scales)

set.seed(123)

# ---- 1) Read + compact modelling table ----
data_path <- "PBResponse_Data.csv"  

pb <- readr::read_csv(data_path, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  transmute(
    record_id  = unique_record_id,
    confirmed  = str_to_upper(str_trim(confirmed_pb)),
    distance_m = readr::parse_number(closest_distance_meter),
    response_raw = str_to_lower(str_trim(response)),
    response_type_raw = str_to_lower(str_trim(behavior_change_response)),
    
    # activity
    activity2_raw = str_to_lower(str_trim(activity2)),
    nearest_activity_raw = str_to_lower(str_trim(nearest_activity)),
    
    # covariates
    season_iow_raw = str_to_upper(str_trim(season_i_ow)),
    survey_method_raw = str_to_lower(str_trim(survey_method_air_land_sea)),
    
    # group composition
    cub_yn_raw = str_to_upper(str_trim(cub_y_n)),
    m_cub = suppressWarnings(as.integer(m_cub)),
    f_cub = suppressWarnings(as.integer(f_cub)),
    u_cub = suppressWarnings(as.integer(u_cub)),
    
    # optional clustering/context (keep if you later decide to use it)
    year = suppressWarnings(as.integer(year)),
    unit = unit,
    site = map_site,
    grid_id = grid_id
  ) %>%
  mutate(
    confirmed = factor(if_else(confirmed == "Y", "yes", "no"), levels = c("no","yes")),
    
    # response: Response / No Response / Not Analyzed
    response = case_when(
      response_raw == "response"      ~ "yes",
      response_raw == "no response"   ~ "no",
      TRUE                            ~ NA_character_   # "not analyzed", blanks, etc.
    ) %>% factor(levels = c("no","yes")),
    response_num = as.integer(response == "yes"),
    
    # response type: only meaningful if response == "yes"
    response_type = case_when(
      response == "yes" & str_detect(response_type_raw, "walk") ~ "walk",
      response == "yes" & str_detect(response_type_raw, "swim") ~ "swim",
      response == "yes" & str_detect(response_type_raw, "run")  ~ "run",
      TRUE                                                      ~ NA_character_
    ) %>% factor(levels = c("walk","swim","run")),
    
    # activity class: prefer Activity2, fill missing from Nearest Activity
    activity_class = case_when(
      activity2_raw == "stationary" ~ "stationary",
      activity2_raw == "mobile"     ~ "mobile",
      nearest_activity_raw %in% c("pad operations", "pad operations ", "pad operations") ~ "stationary",
      nearest_activity_raw %in% c("driving","vessel","aircraft") ~ "mobile",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("stationary","mobile")),
    
    # season: just I vs OW 
    season = case_when(
      season_iow_raw == "I"  ~ "ice",
      season_iow_raw == "OW" ~ "open_water",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("ice","open_water")),
    
    # cubs: use cub flag OR cub counts 
    cubs_count = coalesce(m_cub,0L) + coalesce(f_cub,0L) + coalesce(u_cub,0L),
    cubs_present = (cub_yn_raw == "Y") | (cubs_count > 0),
    
    group_comp = if_else(cubs_present, "female_with_cubs", "no_cubs") %>%
      factor(levels = c("no_cubs","female_with_cubs")),
    
    # survey method (mostly Land, but keep if you want to test)
    survey_method = case_when(
      survey_method_raw %in% c("land") ~ "Land",
      survey_method_raw %in% c("sea")  ~ "Sea",
      survey_method_raw %in% c("air")  ~ "Air",
      survey_method_raw %in% c("hovercraft") ~ "Hovercraft",
      TRUE ~ NA_character_
    ) %>% factor()
  )

# ---- 2) Analysis dataset (core binary model) ----
dat <- pb %>%
  filter(
    confirmed == "yes",
    !is.na(distance_m), distance_m > 0,
    !is.na(response),
    !is.na(activity_class),
    !is.na(season),
    !is.na(group_comp)
  )

dat %>%
  summarise(
    n = n(),
    reaction_rate = mean(response == "yes"),
    dist_min = min(distance_m),
    dist_median = median(distance_m),
    dist_max = max(distance_m)
  ) %>% print()

# ---- 3) Candidate model set + AICc selection ---

# AICc calculator (works for glm; uses logLik df for k)
AICc <- function(model) {
  aic <- AIC(model)
  k   <- attr(logLik(model), "df")
  n   <- stats::nobs(model)
  # Guard against tiny n
  if (is.na(n) || (n - k - 1) <= 0) return(aic)
  aic + (2 * k * (k + 1)) / (n - k - 1)
}

# Build a small, transparent candidate set 
fit_candidate_set <- function(dat, dist_term) {
  f_base   <- as.formula(paste0("response ~ ", dist_term))
  f_act    <- as.formula(paste0("response ~ ", dist_term, " * activity_class"))
  f_core   <- as.formula(paste0("response ~ ", dist_term, " * activity_class + season * group_comp"))
  f_survey <- as.formula(paste0("response ~ ", dist_term, " * activity_class + season * group_comp + survey_method"))
  
  mods <- list(
    base = glm(f_base,   data = dat, family = binomial()),
    act  = glm(f_act,    data = dat, family = binomial()),
    core = glm(f_core,   data = dat, family = binomial())
  )
  
  # Only include survey_method model if survey_method exists and has >1 level
  if ("survey_method" %in% names(dat) && nlevels(droplevels(dat$survey_method)) > 1) {
    mods$survey <- glm(f_survey, data = dat, family = binomial())
  }
  
  mods
}

# Fit distance-form sets
mods_lin <- fit_candidate_set(dat, "I(distance_m/100)")
mods_log <- fit_candidate_set(dat, "log(distance_m)")
mods_spl <- fit_candidate_set(dat, "splines::ns(distance_m, df=4)")

# Combine all models in one named list
mods_all <- c(
  purrr::imap(mods_lin, ~ setNames(list(.x), paste0("lin_", .y))),
  purrr::imap(mods_log, ~ setNames(list(.x), paste0("log_", .y))),
  purrr::imap(mods_spl, ~ setNames(list(.x), paste0("spl_", .y)))
) %>% unlist(recursive = FALSE)

# Model selection table
sel_tbl <- purrr::imap_dfr(mods_all, function(m, nm) {
  tibble(
    model = nm,
    k     = attr(logLik(m), "df"),
    n     = nobs(m),
    AIC   = AIC(m),
    AICc  = AICc(m)
  )
}) %>%
  arrange(AICc) %>%
  mutate(
    delta = AICc - min(AICc),
    weight = exp(-0.5 * delta) / sum(exp(-0.5 * delta))
  )

print(sel_tbl)

# Pick top model
final_model <- mods_all[[sel_tbl$model[1]]]
summary(final_model)

# Using the rule forn "ΔAICc < 2" set:
top_names <- sel_tbl %>% filter(delta < 2) %>% pull(model)
top_models <- mods_all[top_names]
top_weights <- sel_tbl %>% filter(delta < 2) %>% pull(weight)

cat("\nTop set (ΔAICc < 2):\n")
print(sel_tbl %>% filter(delta < 2))

# ---- 4) Model Diagnostics ----
print(performance::check_collinearity(final_model))
print(performance::check_overdispersion(final_model))
# Influence diagnostics 
infl <- influence.measures(final_model)
summary(infl)

# Indices of points flagged as influential by any measure - might need to hcta to Kate about this one to see whether these points warrant a further look. For now looks fine.
which(apply(infl$is.inf, 1, any))

sim <- DHARMa::simulateResiduals(final_model, n = 500)
plot(sim)
DHARMa::testDispersion(sim)
DHARMa::testZeroInflation(sim)
DHARMa::testOutliers(sim)

# AUC
dat <- dat %>%
  mutate(pred_prob = predict(final_model, newdata = dat, type = "response"))
auc <- pROC::roc(response = dat$response_num, predictor = dat$pred_prob, quiet = TRUE) %>% pROC::auc()
cat("\nAUC:", as.numeric(auc), "\n")

roc_obj <- pROC::roc(
  response  = dat$response_num,
  predictor = dat$pred_prob,
  quiet     = TRUE
)

# Build ROC curve
roc_df <- data.frame(
  fpr = 1 - roc_obj$specificities,
  tpr = roc_obj$sensitivities
)

auc_label <- paste0("AUC = ", round(as.numeric(pROC::auc(roc_obj)), 3))

ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_ribbon(aes(ymin = 0, ymax = tpr), fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50", linewidth = 0.7) +
  annotate("text", x = 0.75, y = 0.15, label = auc_label,
           size = 4.5, colour = "#2c7bb6", fontface = "bold") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0.01, 0.01)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0.01, 0.01)) +
  labs(
    x     = "False positive rate (1 − Specificity)",
    y     = "True positive rate (Sensitivity)",
    title = "ROC Curve — binary reaction model"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# ---- 5) Prediction curve (distance-response) ----
make_pred_grid <- function(dat, season_ref = "open_water", group_ref = "no_cubs", survey_ref = NULL) {
  dist_seq <- seq(quantile(dat$distance_m, 0.01), quantile(dat$distance_m, 0.99), length.out = 250)
  
  g <- tidyr::expand_grid(
    distance_m = dist_seq,
    activity_class = levels(dat$activity_class),
    season = factor(season_ref, levels = levels(dat$season)),
    group_comp = factor(group_ref, levels = levels(dat$group_comp))
  )
  
  # include survey_method only if the model uses it
  if ("survey_method" %in% all.vars(formula(final_model))) {
    if (is.null(survey_ref)) survey_ref <- levels(dat$survey_method)[1]
    g <- g %>% mutate(survey_method = factor(survey_ref, levels = levels(dat$survey_method)))
  }
  
  g
}

pred_grid <- make_pred_grid(dat)
pred_link <- predict(final_model, newdata = pred_grid, type = "link", se.fit = TRUE)

pred_df <- pred_grid %>%
  mutate(
    link = pred_link$fit,
    se   = pred_link$se.fit,
    prob = plogis(link),
    lo   = plogis(link - 1.96 * se),
    hi   = plogis(link + 1.96 * se)
  )

p_curve <- ggplot(pred_df, aes(distance_m, prob, colour = activity_class, fill = activity_class)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Distance (m)",
    y = "Predicted probability of behavioural reaction",
    title = "Distance–response curve (binary reaction model)",
    subtitle = "95% CI ribbons; other covariates held at reference levels"
  ) +
  theme_minimal()

print(p_curve)

# ---- 6) Poliucy relevant outputs table (805 m and 1610 m) ----

key_dists <- c(805, 1610)

pred_at_distance <- function(d) {
  nd <- make_pred_grid(dat) %>%
    slice(1) %>%                         
    mutate(distance_m = d)
  
  pr <- predict(final_model, newdata = nd, type = "link", se.fit = TRUE)
  tibble(
    distance_m = d,
    prob = plogis(pr$fit),
    lo   = plogis(pr$fit - 1.96 * pr$se.fit),
    hi   = plogis(pr$fit + 1.96 * pr$se.fit)
  )
}

mgmt_tbl <- map_dfr(key_dists, pred_at_distance) %>%
  mutate(across(c(prob, lo, hi), ~ scales::percent(.x, accuracy = 0.1)))

print(mgmt_tbl)

# ---- 7) Save outputs ----
dir.create("results", showWarnings = FALSE)
write_csv(mgmt_tbl, "results/management_probs_805_1610.csv")
write_csv(as_tibble(sel), "results/model_selection_candidate_set.csv")
ggsave("results/pred_curve.png", p_curve, width = 9, height = 5, dpi = 300)
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

cat("\nDone. Results saved in ./results\n")



