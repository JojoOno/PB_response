############################################################
# Polar Bear Behavioural Response Re-analysis (ESR revision)
# Workflow: unified (primary) binary model + optional reaction-type model
# Author: <your name>
# Data: PBResponse_Data (CSV or Excel) in working directory
############################################################

## 0) Housekeeping ------------------------------------------------------------
rm(list = ls())
set.seed(123)

## 1) Packages (install + load) -----------------------------------------------
pkgs <- c(
  "tidyverse",   # dplyr, ggplot2, tidyr, purrr, readr, tibble, stringr, forcats
  "janitor",     # clean_names()
  "lubridate",   # date handling
  "broom",       # tidy model outputs
  "splines",     # ns()
  "MuMIn",       # AICc model selection + model averaging
  "performance", # model checks
  "DHARMa",      # residual diagnostics
  "pROC",        # AUC
  "yardstick",   # metrics (tidy)
  "patchwork",   # combine plots
  "nnet"         # multinomial (optional)
)

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

invisible(lapply(pkgs, library, character.only = TRUE))

## 2) Read data ---------------------------------------------------------------
# Expect either PBResponse_Data.csv OR PBResponse_Data.xlsx (adjust if needed)
read_pb_data <- function(base_name = "PBResponse_Data") {
  csv_path  <- paste0(base_name, ".csv")
  xlsx_path <- paste0(base_name, ".xlsx")
  
  if (file.exists(csv_path)) {
    readr::read_csv(csv_path, show_col_types = FALSE)
  } else if (file.exists(xlsx_path)) {
    # readxl is not in pkgs by default; install/load if needed
    if (!"readxl" %in% installed.packages()[, "Package"]) install.packages("readxl")
    library(readxl)
    readxl::read_xlsx(xlsx_path)
  } else if (file.exists(base_name)) {
    # if user supplies filename with extension already
    if (grepl("\\.csv$", base_name, ignore.case = TRUE)) {
      readr::read_csv(base_name, show_col_types = FALSE)
    } else if (grepl("\\.xlsx$", base_name, ignore.case = TRUE)) {
      if (!"readxl" %in% installed.packages()[, "Package"]) install.packages("readxl")
      library(readxl)
      readxl::read_xlsx(base_name)
    } else {
      stop("Could not detect file type for PBResponse_Data. Use .csv or .xlsx.")
    }
  } else {
    stop("Could not find PBResponse_Data.csv or PBResponse_Data.xlsx in the working directory.")
  }
}

pb_raw <- read_pb_data("PBResponse_Data")

## 3) Clean column names + light type coercion --------------------------------
pb <- pb_raw %>%
  janitor::clean_names()

# Helper: pick the first column name that exists from a set of candidates
pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

# Candidate columns (based on your screenshot + manuscript terms)
col_confirmed <- pick_col(pb, c("confirmed_pb", "confirm_pb", "confirmed", "pb_confirmed"))
col_distance1 <- pick_col(pb, c("closest_distance_m", "closest_distance_in_meters", "closest_distance_meters",
                                "sighting_distance_m", "sighting_distance_meter", "sighting_distance_in_meters",
                                "distance_m", "distance"))
col_response  <- pick_col(pb, c("response_y_or_n", "response", "reaction", "behavior_change_response"))
col_rtype     <- pick_col(pb, c("response_type", "reaction_type", "response_category", "response_cat"))
col_activity  <- pick_col(pb, c("activity_2", "activity", "activity_type", "encounter_activity"))
col_season    <- pick_col(pb, c("itr_season_and_or_year", "itr_season", "season"))
col_year      <- pick_col(pb, c("year"))
col_groupcomp <- pick_col(pb, c("confirmed_combined_sighting_comp", "group_composition", "sighting_comp"))
col_platform  <- pick_col(pb, c("survey_method_air_land_sea", "survey_method", "platform", "air_land_sea"))
col_site      <- pick_col(pb, c("map_site", "site", "grid_id", "map_grid_id"))
col_operator  <- pick_col(pb, c("unit", "operator", "company", "facility"))

# If the dataset has individual count columns, these may exist:
col_m_cub <- pick_col(pb, c("m_cub", "m_cubs"))
col_f_cub <- pick_col(pb, c("f_cub", "f_cubs"))
col_u_cub <- pick_col(pb, c("u_cub", "u_cubs"))

# Create modelling-friendly columns (intuitive names)
pb2 <- pb %>%
  mutate(
    confirmed_pb = if (!is.na(col_confirmed)) as.character(.data[[col_confirmed]]) else NA_character_,
    distance_m   = if (!is.na(col_distance1)) suppressWarnings(as.numeric(.data[[col_distance1]])) else NA_real_,
    response_raw = if (!is.na(col_response))  as.character(.data[[col_response]])  else NA_character_,
    response_type_raw = if (!is.na(col_rtype)) as.character(.data[[col_rtype]]) else NA_character_,
    activity_raw = if (!is.na(col_activity))  as.character(.data[[col_activity]]) else NA_character_,
    season_raw   = if (!is.na(col_season))    as.character(.data[[col_season]])   else NA_character_,
    year         = if (!is.na(col_year))      suppressWarnings(as.integer(.data[[col_year]])) else NA_integer_,
    group_raw    = if (!is.na(col_groupcomp)) as.character(.data[[col_groupcomp]]) else NA_character_,
    platform_raw = if (!is.na(col_platform))  as.character(.data[[col_platform]]) else NA_character_,
    site_raw     = if (!is.na(col_site))      as.character(.data[[col_site]]) else NA_character_,
    operator_raw = if (!is.na(col_operator))  as.character(.data[[col_operator]]) else NA_character_
  )

# Derive: response (binary), response_type (categorical), activity_class (mobile/stationary), group_comp
pb3 <- pb2 %>%
  mutate(
    # Standardise confirmed flag (if present)
    confirmed_pb = str_to_lower(confirmed_pb),
    confirmed_pb = case_when(
      confirmed_pb %in% c("y", "yes", "true", "1") ~ "yes",
      confirmed_pb %in% c("n", "no", "false", "0") ~ "no",
      TRUE ~ confirmed_pb
    ),
    
    # Binary response: yes/no
    response = str_to_lower(response_raw),
    response = case_when(
      response %in% c("y", "yes", "1", "true", "response", "reacted") ~ "yes",
      response %in% c("n", "no", "0", "false", "no response", "no_response") ~ "no",
      TRUE ~ response
    ),
    response = factor(response, levels = c("no", "yes")),
    
    # Response type (when response == yes)
    response_type = str_to_lower(response_type_raw),
    response_type = str_replace_all(response_type, "\\s+", "_"),
    response_type = case_when(
      response_type %in% c("walk", "walk_away") ~ "walk",
      response_type %in% c("swim", "swim_away") ~ "swim",
      response_type %in% c("run", "run_away")   ~ "run",
      TRUE ~ response_type
    ),
    
    # Activity class (mobile vs stationary) from free text
    activity = str_to_lower(activity_raw),
    activity = str_replace_all(activity, "\\s+", "_"),
    activity_class = case_when(
      str_detect(activity, "air|aircraft|helicopter|plane") ~ "mobile",
      str_detect(activity, "vehicle|truck|car|atv|snow")    ~ "mobile",
      str_detect(activity, "vessel|boat|ship")              ~ "mobile",
      str_detect(activity, "pad|facility|camp|stationary")  ~ "stationary",
      TRUE ~ NA_character_
    ),
    activity_class = factor(activity_class, levels = c("stationary", "mobile")),
    
    # Season (try to standardise)
    season = str_to_lower(season_raw),
    season = case_when(
      str_detect(season, "spring") ~ "spring",
      str_detect(season, "fall|autumn") ~ "fall",
      str_detect(season, "summer") ~ "summer",
      str_detect(season, "winter") ~ "winter",
      TRUE ~ season
    ),
    season = factor(season),
    
    # Group composition:
    # Prefer explicit cub count columns if available, else fallback to group_raw text
    cub_count = {
      cc <- rep(NA_real_, nrow(.))
      if (!is.na(col_m_cub) || !is.na(col_f_cub) || !is.na(col_u_cub)) {
        m <- if (!is.na(col_m_cub)) suppressWarnings(as.numeric(.data[[col_m_cub]])) else 0
        f <- if (!is.na(col_f_cub)) suppressWarnings(as.numeric(.data[[col_f_cub]])) else 0
        u <- if (!is.na(col_u_cub)) suppressWarnings(as.numeric(.data[[col_u_cub]])) else 0
        cc <- m + f + u
      }
      cc
    },
    group_comp = case_when(
      !is.na(cub_count) & cub_count > 0 ~ "female_with_cubs",
      !is.na(cub_count) & cub_count == 0 ~ "no_cubs",
      TRUE ~ NA_character_
    ),
    group_comp = if_else(
      is.na(group_comp) & !is.na(group_raw) & str_detect(str_to_lower(group_raw), "cub"),
      "female_with_cubs",
      group_comp
    ),
    group_comp = if_else(
      is.na(group_comp) & !is.na(group_raw),
      "no_cubs",
      group_comp
    ),
    group_comp = factor(group_comp, levels = c("no_cubs", "female_with_cubs")),
    
    # Basic IDs for possible clustering
    site     = factor(site_raw),
    operator = factor(operator_raw),
    platform = factor(platform_raw)
  )

## 4) Define analysis dataset --------------------------------------------------
dat <- pb3 %>%
  # If confirmed_pb exists, keep confirmed only
  { if (!all(is.na(.$confirmed_pb))) filter(., confirmed_pb == "yes") else . } %>%
  filter(!is.na(distance_m), distance_m > 0) %>%
  filter(!is.na(response)) %>%
  mutate(
    dist100 = distance_m / 100,
    log_dist = log(distance_m),
    # Optional: truncate very large distances if desired (uncomment if needed)
    # distance_m = pmin(distance_m, 10000),
    # dist100 = distance_m / 100,
    # log_dist = log(distance_m),
    response_num = as.integer(response == "yes")
  )

# Quick sanity checks
dat %>%
  summarise(
    n = n(),
    reaction_rate = mean(response == "yes"),
    dist_min = min(distance_m),
    dist_median = median(distance_m),
    dist_max = max(distance_m),
    missing_activity_class = mean(is.na(activity_class)),
    missing_season = mean(is.na(season)),
    missing_group_comp = mean(is.na(group_comp))
  ) %>% print()

## 5) Decide covariate set (keep realistic for 2-month revision) --------------
# Recommended "core" covariates (align with reviewer):
# - distance (main driver)
# - activity_class (mobile vs stationary)
# - season
# - group_comp (cubs vs none)
#
# Add platform/operator/site only if needed (and if sample sizes support it).

dat_core <- dat %>%
  filter(!is.na(activity_class), !is.na(season), !is.na(group_comp))

## 6) Model set: Binary response (primary, unified inference) -----------------
# We'll fit three distance functional forms:
#   A) linear per 100m
#   B) log(distance)
#   C) spline (natural cubic) on distance
#
# We'll allow a small number of biologically meaningful interactions:
#   - distance × activity_class
#   - season × group_comp  (optional)
#
# Use AICc model selection + (optional) model averaging.

options(na.action = "na.fail") # required by MuMIn::dredge

# Global models (keep modest—reviewer wants "candidate set", not a fishing expedition)
m_lin_global <- glm(
  response ~ dist100 * activity_class + season + group_comp + season:group_comp,
  data = dat_core,
  family = binomial()
)

m_log_global <- glm(
  response ~ log_dist * activity_class + season + group_comp + season:group_comp,
  data = dat_core,
  family = binomial()
)

m_spline_global <- glm(
  response ~ ns(distance_m, df = 4) * activity_class + season + group_comp + season:group_comp,
  data = dat_core,
  family = binomial()
)

# Dredge each global model
dd_lin    <- MuMIn::dredge(m_lin_global, trace = FALSE)
dd_log    <- MuMIn::dredge(m_log_global, trace = FALSE)
dd_spline <- MuMIn::dredge(m_spline_global, trace = FALSE)

# Compare best models across distance forms
best_tbl <- bind_rows(
  tibble(form = "linear",  best_aicc = min(dd_lin$AICc),    best_df = dd_lin$df[which.min(dd_lin$AICc)]),
  tibble(form = "log",     best_aicc = min(dd_log$AICc),    best_df = dd_log$df[which.min(dd_log$AICc)]),
  tibble(form = "spline",  best_aicc = min(dd_spline$AICc), best_df = dd_spline$df[which.min(dd_spline$AICc)])
) %>% arrange(best_aicc)

print(best_tbl)

# Pick the top-ranked distance form and models within ΔAICc < 2
pick_models <- function(dd) {
  mods <- MuMIn::get.models(dd, subset = delta < 2)
  list(
    dredge = dd,
    models = mods,
    avg = if (length(mods) > 1) MuMIn::model.avg(mods) else NULL,
    top = mods[[1]]
  )
}

res_lin    <- pick_models(dd_lin)
res_log    <- pick_models(dd_log)
res_spline <- pick_models(dd_spline)

# Decide final by the best AICc form
final_res <- switch(best_tbl$form[1],
                    "linear" = res_lin,
                    "log"    = res_log,
                    "spline" = res_spline
)

final_model <- final_res$top
summary(final_model)

# If there are multiple near-tied models, you can use model averaging
if (!is.null(final_res$avg)) {
  cat("\nModel averaging across ΔAICc < 2 models\n")
  print(summary(final_res$avg))
}

## 7) Model validation / robustness checks ------------------------------------
# Collinearity, overdispersion, influential points, and DHARMa residual diagnostics
performance::check_collinearity(final_model) %>% print()
performance::check_overdispersion(final_model) %>% print()
performance::check_outliers(final_model) %>% print()
performance::check_influential(final_model) %>% print()

# DHARMa residual simulation
sim <- DHARMa::simulateResiduals(final_model, n = 500)
plot(sim)
DHARMa::testDispersion(sim)
DHARMa::testZeroInflation(sim)   # for binomial, flags issues if excess zeros beyond expectation
DHARMa::testOutliers(sim)

# Predictive discrimination: AUC
dat_core <- dat_core %>%
  mutate(
    pred_prob = predict(final_model, type = "response")
  )

auc <- pROC::roc(response = dat_core$response_num, predictor = dat_core$pred_prob, quiet = TRUE) %>%
  pROC::auc()

cat("\nAUC:", as.numeric(auc), "\n")

## 8) Prediction plots (ggplot2 / tidyverse) ----------------------------------
# Create a smooth distance grid and predict by activity_class (and optionally season/group)
make_pred_grid <- function(dat, dist_seq = NULL) {
  if (is.null(dist_seq)) {
    dist_seq <- seq(
      from = quantile(dat$distance_m, 0.01, na.rm = TRUE),
      to   = quantile(dat$distance_m, 0.99, na.rm = TRUE),
      length.out = 250
    )
  }
  
  expand_grid(
    distance_m = dist_seq,
    activity_class = levels(dat$activity_class),
    season = levels(dat$season)[1],         # choose a reference level (can expand if you want)
    group_comp = levels(dat$group_comp)[1]  # choose a reference level
  ) %>%
    mutate(
      dist100 = distance_m / 100,
      log_dist = log(distance_m)
    )
}

pred_grid <- make_pred_grid(dat_core)

# Predict on link scale for confidence intervals
pred_link <- predict(final_model, newdata = pred_grid, type = "link", se.fit = TRUE)

pred_df <- pred_grid %>%
  mutate(
    link = pred_link$fit,
    se   = pred_link$se.fit,
    prob = plogis(link),
    lo   = plogis(link - 1.96 * se),
    hi   = plogis(link + 1.96 * se)
  )

p_curve <- ggplot(pred_df, aes(x = distance_m, y = prob, colour = activity_class, fill = activity_class)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Distance (m)",
    y = "Predicted probability of behavioural reaction",
    colour = "Activity class",
    fill = "Activity class",
    title = "Distance–response curve (binary reaction model)",
    subtitle = "Shaded bands: 95% confidence intervals (on link scale)"
  ) +
  theme_minimal()

p_curve

# Optional: add points (binned) for visual calibration
binned <- dat_core %>%
  mutate(dist_bin = cut(distance_m, breaks = quantile(distance_m, probs = seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE)) %>%
  group_by(activity_class, dist_bin) %>%
  summarise(
    dist_mid = median(distance_m, na.rm = TRUE),
    obs_rate = mean(response == "yes"),
    n = n(),
    .groups = "drop"
  )

p_curve_binned <- p_curve +
  geom_point(
    data = binned,
    aes(x = dist_mid, y = obs_rate, size = n, colour = activity_class),
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  scale_size_continuous(name = "N (bin)")

p_curve_binned

## 9) Management-style thresholds (distance at target probabilities) ----------
# Works for any model form by using prediction + root finding.
# Note: requires monotonic-ish relationship; if spline gets wiggly, keep df small.
predict_prob_at <- function(model, activity_class, season, group_comp, distance_m) {
  nd <- tibble(
    distance_m = distance_m,
    dist100 = distance_m / 100,
    log_dist = log(distance_m),
    activity_class = factor(activity_class, levels = levels(dat_core$activity_class)),
    season = factor(season, levels = levels(dat_core$season)),
    group_comp = factor(group_comp, levels = levels(dat_core$group_comp))
  )
  as.numeric(predict(model, newdata = nd, type = "response"))
}

find_distance_for_p <- function(model, target_p, activity_class, season, group_comp,
                                lower = 50, upper = 10000) {
  f <- function(d) predict_prob_at(model, activity_class, season, group_comp, d) - target_p
  # Ensure the function crosses 0 in [lower, upper]
  if (sign(f(lower)) == sign(f(upper))) return(NA_real_)
  uniroot(f, lower = lower, upper = upper)$root
}

targets <- c(0.01, 0.05, 0.10) # 1%, 5%, 10% reaction probabilities (edit as needed)
ref_season <- levels(dat_core$season)[1]
ref_group  <- levels(dat_core$group_comp)[1]

thresholds <- expand_grid(
  activity_class = levels(dat_core$activity_class),
  target_p = targets
) %>%
  mutate(
    distance_m = map2_dbl(activity_class, target_p,
                          ~ find_distance_for_p(final_model, .y, .x, ref_season, ref_group)
    )
  )

print(thresholds)

## 10) OPTIONAL: Reaction type model (if counts support it) -------------------
# If you want to address the reviewer’s “consider response type”:
# You can model categories (none/walk/swim/run) with multinomial regression.
# This can get fragile when some categories are rare or interactions are added.
#
# Minimal approach: use a multinomial with a small covariate set and no huge interactions.
#
# Uncomment to run.

# dat_multi <- dat_core %>%
#   mutate(
#     response_cat = case_when(
#       response == "no" ~ "none",
#       response == "yes" & response_type %in% c("walk", "swim", "run") ~ response_type,
#       TRUE ~ "other"
#     ),
#     response_cat = factor(response_cat)
#   ) %>%
#   filter(response_cat %in% c("none", "walk", "swim", "run")) %>%
#   droplevels()
#
# # Simple multinomial model (nnet::multinom)
# m_multinom <- nnet::multinom(
#   response_cat ~ dist100 + activity_class + season + group_comp,
#   data = dat_multi,
#   trace = FALSE
# )
#
# summary(m_multinom)
# AIC(m_multinom)
#
# # Predicted class probabilities vs distance (for one season/group reference)
# pred_grid2 <- make_pred_grid(dat_multi) %>%
#   mutate(dist100 = distance_m / 100)
#
# probs <- predict(m_multinom, newdata = pred_grid2, type = "probs")
# probs_df <- bind_cols(pred_grid2, as_tibble(probs)) %>%
#   pivot_longer(cols = -c(distance_m, activity_class, season, group_comp, dist100, log_dist),
#                names_to = "class", values_to = "prob")
#
# p_multi <- ggplot(probs_df, aes(distance_m, prob, colour = class)) +
#   geom_line(linewidth = 1) +
#   facet_wrap(~ activity_class) +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#   labs(
#     x = "Distance (m)",
#     y = "Predicted probability",
#     colour = "Response category",
#     title = "Multinomial predicted response type vs distance"
#   ) +
#   theme_minimal()
#
# p_multi

## 11) Outputs (tables + plots) ----------------------------------------------
dir.create("results", showWarnings = FALSE)

# Model selection tables
write_csv(as_tibble(dd_lin),    "results/model_selection_linear.csv")
write_csv(as_tibble(dd_log),    "results/model_selection_log.csv")
write_csv(as_tibble(dd_spline), "results/model_selection_spline.csv")
write_csv(best_tbl,             "results/best_distance_form_compare.csv")
write_csv(thresholds,           "results/threshold_distances.csv")

# Save plots
ggsave("results/pred_curve.png", p_curve, width = 9, height = 5, dpi = 300)
ggsave("results/pred_curve_binned.png", p_curve_binned, width = 9, height = 5, dpi = 300)

# Session info for reproducibility
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

cat("\nDone. Results saved in ./results\n")
