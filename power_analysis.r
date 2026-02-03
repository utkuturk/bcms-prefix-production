# =============================================================================
# Power Analysis for Morphological Priming Study (BCMS Prefixed Verbs)
# Based on effect sizes from Creemers et al. (2020, J Mem Lang)
# Following Solomon Kurz's Bayesian power analysis workflow:
# https://solomonkurz.netlify.app/blog/bayesian-power-analysis-part-i/
# =============================================================================
# The inefficient version is given after the code as a commented out. 

library(tidyverse)
library(brms)
library(cmdstanr)

# Set today's seed for reproducibility
set.seed(02022026)

# parameters from creemers study 

# Prior means (point estimates from Creemers)
prior_means <- list(
  intercept = 6.81,
  beta_MS = -0.12,
  beta_M  = -0.10,
  beta_Ph = -0.02,
  beta_pseudo = -0.04,
  sd_subject = 0.076,
  sd_item = 0.030,
  sd_residual = 0.126
)

# Prior SDs (uncertainty around those estimates)
prior_sds <- list(
  intercept = 0.5,
  beta_MS = 0.05,
  beta_M  = 0.05,
  beta_Ph = 0.05,
  beta_pseudo = 0.05,
  sd_subject = 0.03,
  sd_item = 0.02,
  sd_residual = 0.05
)

# simulate ta data

sim_data <- function(seed, n_subs, n_items_per_cond = 6, prior_means, prior_sds) {
  
  set.seed(seed)
  
  sampled <- list(
    intercept    = rnorm(1, prior_means$intercept, prior_sds$intercept),
    beta_MS      = rnorm(1, prior_means$beta_MS, prior_sds$beta_MS),
    beta_M       = rnorm(1, prior_means$beta_M, prior_sds$beta_M),
    beta_Ph      = rnorm(1, prior_means$beta_Ph, prior_sds$beta_Ph),
    beta_pseudo  = rnorm(1, prior_means$beta_pseudo, prior_sds$beta_pseudo),
    # SDs must be positive - use abs() or truncnorm
    sd_subject   = abs(rnorm(1, prior_means$sd_subject, prior_sds$sd_subject)),
    sd_item      = abs(rnorm(1, prior_means$sd_item, prior_sds$sd_item)),
    sd_residual  = abs(rnorm(1, prior_means$sd_residual, prior_sds$sd_residual))
  )
  
  conditions <- c("MS", "M", "Pseudo", "Ph")
  n_items <- n_items_per_cond * length(conditions)
  
  item_bank <- tibble(
    item = paste0("i_", 1:n_items),
    target_type = rep(conditions, each = n_items_per_cond)
  )
  
  d <- expand_grid(
    subject = paste0("s_", 1:n_subs),
    item = item_bank$item
  ) %>%
    left_join(item_bank, by = "item") %>%
    mutate(
      subj_num = as.numeric(gsub("s_", "", subject)),
      item_num = as.numeric(gsub("i_", "", item)),
      is_target = (subj_num %% 2) == (item_num %% 2),
      condition = if_else(is_target, target_type, "Control")
    )
  
  subj_re <- tibble(
    subject = paste0("s_", 1:n_subs),
    re_subj = rnorm(n_subs, 0, sampled$sd_subject)
  )
  
  item_re <- tibble(
    item = paste0("i_", 1:n_items),
    re_item = rnorm(n_items, 0, sampled$sd_item)
  )
  
  
  d <- d %>%
    left_join(subj_re, by = "subject") %>%
    left_join(item_re, by = "item") %>%
    mutate(
      cond_MS = as.numeric(condition == "MS"),
      cond_M = as.numeric(condition == "M"),
      cond_Pseudo = as.numeric(condition == "Pseudo"),
      cond_Ph = as.numeric(condition == "Ph"),
      
      mu = sampled$intercept +
        sampled$beta_MS * cond_MS +
        sampled$beta_M * cond_M +
        sampled$beta_pseudo * cond_Pseudo +
        sampled$beta_Ph * cond_Ph +
        re_subj + re_item,
      
      log_rt = rnorm(n(), mu, sampled$sd_residual)
    ) %>%
    select(subject, item, condition, cond_MS, cond_M, cond_Pseudo, cond_Ph, log_rt)
  
  return(list(data = d, sampled_params = sampled))
}


# Create priors based on Creemers et al. parameters
fit_priors <- c(
  # Intercept: centered on empirical estimate with moderate uncertainty
  set_prior("normal(6.81, 0.5)", class = "Intercept"),
  
  # Fixed effects: centered on expected values with moderate uncertainty
  set_prior("normal(-0.12, 0.05)", class = "b", coef = "cond_MS"),
  set_prior("normal(-0.10, 0.05)", class = "b", coef = "cond_M"),
  set_prior("normal(-0.04, 0.05)", class = "b", coef = "cond_Pseudo"),
  set_prior("normal(-0.02, 0.05)", class = "b", coef = "cond_Ph"),
  
  # Random effects SDs: weakly informative around empirical estimates
  # Subject SD
  set_prior("normal(0.076, 0.05)", class = "sd", coef = "Intercept", group = "subject"),
  # Item SD
  set_prior("normal(0.030, 0.02)", class = "sd", coef = "Intercept", group = "item"),
  
  # Residual SD
  set_prior("normal(0.126, 0.05)", class = "sigma")
)

# compile the model once


# Initial data for compilation
init_result <- sim_data(
  seed = 02022026, n_subs = 50, 
  prior_means = prior_means, prior_sds = prior_sds
)
d_init <- init_result$data

# Formula
f <- log_rt ~ cond_MS + cond_M + cond_Pseudo + cond_Ph + (1 | subject) + (1 | item)

# Fit initial model - THIS COMPILES THE STAN CODE ONCE
fit_init <- brm(
  formula = f,
  data = d_init,
  family = gaussian(),
  prior = fit_priors,
  chains = 4,
  cores = 4,
  threads = threading(8),
  iter = 4000,
  warmup = 2000,
  seed = 02022026,
  backend = "cmdstanr"
)
# simulate function

sim_d_and_fit <- function(seed, n_subs, fit, prior_means, prior_sds) {  # Simulate new data
  
  result <- sim_data(
    seed = seed, 
    n_subs = n_subs,
    prior_means = prior_means, 
    prior_sds = prior_sds
  )
  
  d_new <- result$data
  true_params <- result$sampled_params
  
  # Update fit with new data (uses pre-compiled Stan code!)
  fit_new <- update(
    fit,
    newdata = d_new,
    seed = seed,
    chains = 4,
    cores = 32,
    iter = 4000,
    warmup = 2000,
    refresh = 0  
  )
  
  # Extract fixed effects summary
  fe <- fixef(fit_new)
  
  tibble(
    condition = c("MS", "M", "Pseudo", "Ph"),
    estimate = fe[c("cond_MS", "cond_M", "cond_Pseudo", "cond_Ph"), "Estimate"],
    Q2.5 = fe[c("cond_MS", "cond_M", "cond_Pseudo", "cond_Ph"), "Q2.5"],
    Q97.5 = fe[c("cond_MS", "cond_M", "cond_Pseudo", "cond_Ph"), "Q97.5"],
    true_value = c(
      true_params$beta_MS, 
      true_params$beta_M,
      true_params$beta_pseudo, 
      true_params$beta_Ph
    ),
    significant = (Q2.5 > 0 & Q97.5 > 0) | (Q2.5 < 0 & Q97.5 < 0),
    ci_width = Q97.5 - Q2.5,
    covers_true = (Q2.5 <= true_value) & (true_value <= Q97.5)
  )
}

# power calculation 


calculate_power <- function(n_subs, n_sim, fit, prior_means, prior_sds) {
  
  cat(sprintf("Running %d simulations for n = %d...\n", n_sim, n_subs))
  
  sim_results <- map_dfr(
    1:n_sim,
    ~ {
      if (.x %% 10 == 0) cat(sprintf("  Simulation %d/%d\n", .x, n_sim))
      sim_d_and_fit(
        seed = .x, 
        n_subs = n_subs, 
        fit = fit,
        prior_means = prior_means, 
        prior_sds = prior_sds
      )
    },
    .id = "sim"
  )
  
  power_summary <- sim_results %>%
    group_by(condition) %>%
    summarize(
      power = mean(significant),
      mean_estimate = mean(estimate),
      mean_true = mean(true_value),
      sd_true = sd(true_value),  # Shows variation in sampled parameters
      bias = mean(estimate - true_value),
      coverage = mean(covers_true),  # Should be ~0.95
      mean_ci_width = mean(ci_width),
      n_sims = n()
    ) %>%
    mutate(n_participants = n_subs)
  
  return(power_summary)
}

# power across sample sizes 


sample_sizes <- c(40, 80, 120, 160)
n_simulations <- 1000  # Use 500-1000 for final analysis

power_results <- map_dfr(sample_sizes, function(n) {
  calculate_power(
    n_subs = n,
    n_sim = n_simulations,
    fit = fit_init,
    prior_means = prior_means,
    prior_sds = prior_sds
  )
})


# results

print(power_results, n = Inf)

write_csv(power_results, "power_results_prior_predictive.csv")

# Power curve
power_curve <- ggplot(power_results, aes(x = n_participants, y = power, color = condition)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0.90, linetype = "dotted", color = "darkred", linewidth = 0.8) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = sample_sizes) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Bayesian Power Analysis (Prior Predictive)",
    subtitle = sprintf("Parameters sampled from priors | %d simulations per N", n_simulations),
    x = "Number of Participants",
    y = "Power (95% CI excludes 0)",
    color = "Condition",
    caption = "Red dashed = 80% | Red dotted = 90%"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave("power_curve.png", power_curve, width = 10, height = 6, dpi = 300)



# Coverage plot (calibration check)
coverage_plot <- ggplot(power_results, aes(x = n_participants, y = coverage, color = condition)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  scale_y_continuous(labels = scales::percent, limits = c(0.8, 1)) +
  scale_x_continuous(breaks = sample_sizes) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Coverage: Does 95% CI contain true value?",
    subtitle = "Should be ~95% if model is well-calibrated",
    x = "Number of Participants",
    y = "Coverage",
    color = "Condition"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("coverage_plot.png", coverage_plot, width = 10, height = 6, dpi = 300)

# summary 

cat("\n\n========================================\n")
cat("Summary\n")
cat("========================================\n\n")


# Minimum N for 80% power
min_n_80 <- power_results %>%
  filter(power >= 0.80) %>%
  group_by(condition) %>%
  slice_min(n_participants) %>%
  select(condition, n_participants, power)

cat("Minimum N for 80% power:\n")
print(min_n_80)

# Power at each sample size
for (n in sample_sizes) {
  cat(sprintf("\n\nPower at n = %d:\n", n))
  power_results %>% filter(n_participants == n) %>% print()
}


######## INEFFICIENT OLD VERSION
# =============================================================================
# Power Analysis for Morphological Priming Study (BCMS Prefixed Verbs)
# Based on effect sizes from Creemers et al. (2020, J Mem Lang)
# =============================================================================


# # Load required packages
# library(tidyverse)
# library(brms)
# library(lme4)
# library(cmdstanr)


# # part 1: get values from their exp1
# params <- list(
#   # Fixed effects (relative to Control baseline)
#   intercept = 6.81,
#   beta_MS = -0.12, # Morphologically + Semantically related
#   beta_M  = -0.10, # Morphologically related only (opaque)
#   beta_Ph = -0.02, # Phonologically related (not significant)
#   # additionall we have: pseudo-morphological condition
#   # let's assume that it is somewhere between phonological and morphological.
#   beta_pseudo = -0.04, # conservative estimate for pseudo-morphological
  
#   # random effects from creemers
#   sd_subject = 0.076,
#   sd_item = 0.030,
#   sd_residual = 0.126
# )


# # Create priors based on Creemers et al. parameters
# mypriors <- c(
#   # Intercept: centered on empirical estimate with moderate uncertainty
#   set_prior("normal(6.81, 0.5)", class = "Intercept"),
  
#   # Fixed effects: centered on expected values with moderate uncertainty
#   set_prior("normal(-0.12, 0.05)", class = "b", coef = "cond_MS"),
#   set_prior("normal(-0.10, 0.05)", class = "b", coef = "cond_M"),
#   set_prior("normal(-0.04, 0.05)", class = "b", coef = "cond_Pseudo"),
#   set_prior("normal(-0.02, 0.05)", class = "b", coef = "cond_Ph"),
  
#   # Random effects SDs: weakly informative around empirical estimates
#   # Subject SD
#   set_prior("normal(0.076, 0.05)", class = "sd", coef = "Intercept", group = "subject"),
#   # Item SD
#   set_prior("normal(0.030, 0.02)", class = "sd", coef = "Intercept", group = "item"),
  
#   # Residual SD
#   set_prior("normal(0.126, 0.05)", class = "sigma")
# )


# # design parameters
# design <- list(
#   n_items = 30, 
#   n_conditions = 5,
#   items_per_cond = 6
#   # (5 per priming condition + 10 controls)
# )

# # generating the data 
# create_experiment_data <- function(n_subs, n_items_per_cond = 6, params) {
#   # 1. Define which items belong to which "Priming Type"
#   conditions <- c("MS", "M", "Pseudo", "Ph")
  
#   # Create the "Item Bank"
#   # Each target condition gets a set of items
#   item_bank <- data.frame(
#     item = paste0("i_", 1:(n_items_per_cond * length(conditions))),
#     target_type = rep(conditions, each = n_items_per_cond)
#   )
  
#   # 2. Create the Subject design
#   # Assign subjects to "List A" or "List B" (the counterbalancing)
#   data <- expand.grid(
#     subject = paste0("s_", 1:n_subs),
#     item = item_bank$item
#   ) %>%
#     left_join(item_bank, by = "item") %>%
#     group_by(subject) %>%
#     mutate(
#       # Subj 1, 3, 5... see Items 1-3 as Target, 4-6 as Control
#       # Subj 2, 4, 6... see Items 1-3 as Control, 4-6 as Target
#       subj_num = as.numeric(gsub("s_", "", subject)),
#       item_num = as.numeric(gsub("i_", "", item)),
      
#       # Logic: if (subj is even AND item is even) OR (subj is odd AND item is odd) -> Target
#       # This ensures a 50/50 split of Target vs Control for every item
#       is_target = (subj_num %% 2) == (item_num %% 2),
#       condition = if_else(is_target, target_type, "Control")
#     ) %>%
#     ungroup()
  
#   # 3. Apply the treatment coding for the model
#   data <- data %>%
#     mutate(
#       cond_MS = as.numeric(condition == "MS"),
#       cond_M = as.numeric(condition == "M"),
#       cond_Pseudo = as.numeric(condition == "Pseudo"),
#       cond_Ph = as.numeric(condition == "Ph")
#     )
  
#   return(data)
# }

# # simulaton function
# simulate_and_fit <- function(n_subs, n_items_per_cond = 6, mypriors, use_brms_fit = FALSE) {
  
#   # Step 1: Create dummy data for brms prior sampling
#   # create_experiment_data expects params, but we don't need it for structure
#   dummy_data <- data.frame(
#     subject = rep(paste0("s_", 1:2), each = 10),
#     item = rep(paste0("i_", 1:10), 2),
#     cond_MS = 0,
#     cond_M = 0,
#     cond_Pseudo = 0,
#     cond_Ph = 0,
#     log_rt = 0
#   )
  
#   # Step 2: Sample from priors only
#   prior_model <- brm(
#     log_rt ~ cond_MS + cond_M + cond_Pseudo + cond_Ph + 
#       (1 | subject) + (1 | item),
#     data = dummy_data,
#     family = gaussian(),
#     prior = mypriors,
#     cores = 32, threads = threading(4),
#     sample_prior = "only",
#     chains = 4, iter = 4000,
#     silent = 2, refresh = 0,
#     backend = "cmdstanr"
#   )
  
#   # Step 3: Create new experimental data structure
#   # Use a dummy params object since create_experiment_data needs it
#   dummy_params <- list(intercept = 0, beta_MS = 0, beta_M = 0, 
#                        beta_Ph = 0, beta_pseudo = 0,
#                        sd_subject = 0, sd_item = 0, sd_residual = 0)
                       
#   new_experiment <- create_experiment_data(n_subs, n_items_per_cond, dummy_params)
  
#   # Step 4: Generate data from prior predictive distribution
#   yrep <- posterior_predict(
#     prior_model,
#     newdata = new_experiment,
#     allow_new_levels = TRUE,
#     sample_new_levels = "gaussian"
#   )
  
#   # Use first posterior draw

#   new_experiment$log_rt <- yrep[sample(1:nrow(yrep), 1), ]
  
#   # Step 5: Fit model to generated data
#   if (use_brms_fit) {
#     fit <- brm(
#       log_rt ~ cond_MS + cond_M + cond_Pseudo + cond_Ph + 
#         (1 | subject) + (1 | item),
#       data = new_experiment,
#       family = gaussian(),
#       prior = c(
#         set_prior("normal(0, 1)", class = "b"),
#         set_prior("exponential(1)", class = "sd")
#       ),
#       chains = 2, iter = 1000, warmup = 500,
#       silent = 2, refresh = 0
#     )
    
#     post_summary <- fixef(fit)
    
#     results <- data.frame(
#       condition = c("MS", "M", "Pseudo", "Ph"),
#       estimate = post_summary[c("cond_MS", "cond_M", "cond_Pseudo", "cond_Ph"), "Estimate"],
#       lower = post_summary[c("cond_MS", "cond_M", "cond_Pseudo", "cond_Ph"), "Q2.5"],
#       upper = post_summary[c("cond_MS", "cond_M", "cond_Pseudo", "cond_Ph"), "Q97.5"],
#       significant = post_summary[c("cond_MS", "cond_M", "cond_Pseudo", "cond_Ph"), "Q97.5"] < 0
#     )
#   } else {
#     # FIX 3: Add bobyqa optimizer to address convergence warnings
#     fit <- lmer(
#       log_rt ~ cond_MS + cond_M + cond_Pseudo + cond_Ph + 
#         (1 | subject) + (1 | item),
#       data = new_experiment,
#       REML = FALSE,
#       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
#     )
    
#     coef_summary <- summary(fit)$coefficients
    
#     # Make sure we're creating the results dataframe correctly
#     results <- data.frame(
#       condition = c("MS", "M", "Pseudo", "Ph"),
#       estimate = coef_summary[c("cond_MS", "cond_M", "cond_Pseudo", "cond_Ph"), "Estimate"],
#       se = coef_summary[c("cond_MS", "cond_M", "cond_Pseudo", "cond_Ph"), "Std. Error"],
#       t_value = coef_summary[c("cond_MS", "cond_M", "cond_Pseudo", "cond_Ph"), "t value"],
#       stringsAsFactors = FALSE  # Important!
#     )
    
#     results$p_value <- 2 * pnorm(-abs(results$t_value))
#     results$significant <- results$p_value < 0.05
#   }
  
#   return(results)
# }



# ###########
# ## power funciton
# calculate_power <- function(n_subs, n_items_per_cond, mypriors, n_sim = 1000, use_brms = FALSE) {
  
#   sim_results <- replicate(n_sim, {
#     simulate_and_fit(n_subs, n_items_per_cond, mypriors, use_brms_fit = use_brms)
#   }, simplify = FALSE)
  
#   # Combine results
#   power_summary <- bind_rows(sim_results) %>%
#     group_by(condition) %>%
#     summarize(
#       power = mean(significant, na.rm = TRUE),
#       mean_estimate = mean(estimate, na.rm = TRUE),
#       sd_estimate = sd(estimate, na.rm = TRUE),
#       n_sims = n()
#     ) %>%
#     mutate(n_participants = n_subs)
  
#   return(power_summary)
# }



# # running across sizes

# # Sample sizes to test
# sample_sizes <- c(120)

# # Number of simulations (use 1000 for final analysis, 100 for testing)
# n_simulations <- 1000

# # Run power analysis
# power_results <- lapply(sample_sizes, function(n) {
#   calculate_power(
#     n_subs = n,
#     n_items_per_cond = 6,
#     mypriors = mypriors,
#     n_sim = n_simulations,
#     use_brms = FALSE  # Set TRUE for Bayesian approach (slower)
#   )
# })

# power_table <- bind_rows(power_results)



# power_curve <- ggplot(power_table, aes(x = n_participants, y = power, color = condition)) +
#  geom_line(size = 1.2) +
#  geom_point(size = 3) +
#  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red", size = 0.8) +
#  geom_hline(yintercept = 0.90, linetype = "dotted", color = "darkred", size = 0.8) +
#  scale_y_continuous(
#    labels = scales::percent,
#    limits = c(0, 1),
#    breaks = seq(0, 1, 0.1)
#  ) +
#  scale_x_continuous(breaks = sample_sizes) +
#  scale_color_brewer(palette = "Set1") +
#  labs(
#    title = "Power Analysis: Morphological Priming Study",
#    subtitle = sprintf("Based on Creemers et al. (2020) effect sizes | %d simulations per sample size", n_simulations),
#    x = "Number of Participants",
#    y = "Statistical Power",
#    color = "Condition",
#    caption = "Red dashed = 80% power | Red dotted = 90% power"
#  ) +
#  theme_minimal(base_size = 12) +
#  theme(
#    legend.position = "bottom",
#    plot.title = element_text(face = "bold")
#  )

# ggsave("power_curve.png", power_curve, width = 10, height = 6, dpi = 300)

# # Print summary results

# # Effect size info
# cat(sprintf("Target effect size (Pseudo condition): β = %.2f (log RT scale)\n", params$beta_pseudo))
# cat(sprintf("This corresponds to ~%.0f ms facilitation\n\n", abs(params$beta_pseudo) * 1000 * exp(params$intercept) / 1000))

# # Minimum N for 80% power
# min_n_80 <- power_table %>%
#   filter(power >= 0.80) %>%
#   group_by(condition) %>%
#   slice_min(n_participants) %>%
#   select(condition, n_participants, power)

# cat("Minimum N for 80% power:\n")
# if (nrow(min_n_80) > 0) {
#   print(min_n_80)
# } else {
#   cat("80% power not achieved at any tested sample size for any condition.\n")
#   cat("Consider testing larger sample sizes.\n")
# }

# # Check which conditions didn't reach 80% power
# conditions_tested <- unique(power_table$condition)
# conditions_reached <- unique(min_n_80$condition)
# conditions_not_reached <- setdiff(conditions_tested, conditions_reached)
# if (length(conditions_not_reached) > 0) {
#   cat(sprintf("\nConditions not reaching 80%% power at tested sample sizes: %s\n", 
#               paste(conditions_not_reached, collapse = ", ")))
# }

# # Power at each sample size
# for (n in sample_sizes) {
#   cat(sprintf("\n\nPower at n = %d:\n", n))
#   power_table %>% filter(n_participants == n) %>% print()
# }
