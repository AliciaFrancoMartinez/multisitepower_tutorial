# ============================================================================ #
# Script: A Tutorial on Multilevel Power Analysis 
#         for Multisite Studies with Repeated Measures 
# Author: Alicia Franco-Martínez
# Date:   2026-10-02
#
#
# Workflow:
# 1. Determine a preliminary study architecture
# 2. Identify and model your estimand
# 3. Consider the population you want to generalize to
# 4. Calculate power for this architecture with {simr}
# 5. Define a feasible architecture space and optimize it with {mlpwr}
# 6. Explore alternative effect-size and heterogeneity assumptions
# 7. Assess the selected architecture using MDES and precision
#
# Abbreviations used in object names:
#   partis / parti = participants / participant
#
# Tip: run this script in a fresh R session, 
# so that no objects from previous work interfere with it.
# ============================================================================ #


# ---- 0. Packages -------------------------------------------------------------

# install.packages(c("lme4", "lmerTest", "simr", "mlpwr"))

# *Note: while running this script you will probably see many
#        "boundary (singular) fit" messages, and {simr} may report some
#        "warnings" (these are convergence warnings, not the boundary
#        messages). Both usually come from the same source: with few sites,
#        the model cannot estimate well the variances and correlation of the
#        site-level random effects, so some of them are estimated as zero.
#        They do not affect the estimate of the effect of "condition" nor
#        the decision of its test, so they do not distort the power estimate.
#        See the tips at the end of Section 4 if they are very frequent.


# ---- 1. Determine a preliminary study architecture ---------------------------

# We begin with fixed level sizes before introducing an architecture space.
trials_size <- 32
partis_size <- 16
sites_size <- 12

# Generate all combinations of trials, participants, and sites (6,144 rows)
architecture <- expand.grid(trial = 1:trials_size,
                            participant = factor(1:partis_size),
                            site = factor(1:sites_size))

# *Note: Levels 2 and 3 should be 'factor' objects.                     
# You can chech it using str()
str(architecture)

# *Note: participants are numbered 1, 2, ..., partis_size within EACH site, so
#        participant "1" in site 1 is a different person from participant "1" 
#        in site 2. This is why the model below uses "site:participant" 
#        (participants nested within sites) instead of just "participant".
#        IMPORTANT! IF YOU DON'T SET site:participant, the model would think
#        that participant 1 in lab 1 is also participant 1 in lab 2!!!

# Predictor (independent variable or factor): "condition"
condition_pattern <- rep(c(-0.5, 0.5), length.out = trials_size)

# *Note: This coding (-0.5, 0.5) makes the intercept the grand mean 
#        and beta_condition the mean difference between conditions.

# Add the condition variable to the architecture design
architecture$condition <- rep(condition_pattern, 
                              times = partis_size * sites_size)


# ---- 2. Identify and model your estimand -------------------------------------

# For the tutorial example:
# Your estimand is the average trial-level effect of "condition": beta_condition


# ---- 3. Consider the population you want to generalize to --------------------

# For the tutorial example:
# You want to generalize beyond the recruited participants and sites, 
# so both levels are modeled as random:

# response ~ condition +
#            (1 + condition | site:participant) +
#            (1 + condition | site)


# ---- 4. Calculate power for this preliminary architecture with {simr} --------

# ----- a) Select justifiable values for the model parameters ----

# EQUATION OF THE MODEL:
# response = beta_intercept + beta_condition*condition +          # fixed effects
#            u_parti_intercept + u_parti_condition*condition +    # participant random effects
#            u_site_intercept + u_site_condition*condition +      # site random effects
#            epsilon                                              # residual error
#
# *Note: each parameter is named after the model term it belongs to
#        ("intercept" or the predictor "condition"), so the naming extends
#        directly to models with more predictors: For example, 
#        beta_gender, sd_parti_gender, for a predictor named "gender".


# PARAMETERS (the following values should be theoretically justified)

# Fixed effects
beta_intercept <- 0        # fixed intercept
beta_condition <- 0.15     # fixed effect of "condition" (your estimand)

# Random effects (their standard deviations: sd)
sd_parti_intercept <- 0.50 # parti-level random intercept
sd_parti_condition <- 0.55 # parti-level random slope of "condition"

sd_site_intercept <- 0.20  # site-level random intercept
sd_site_condition <- 0.15  # site-level random slope of "condition"

sd_residual <- 1.00        # residual error

# Variance-covariance matrices (covariances between effects fixed to 0)
V_parti <- matrix(c(sd_parti_intercept^2, 0,
                    0, sd_parti_condition^2),
                  nrow = 2, byrow = TRUE)

V_site <- matrix(c(sd_site_intercept^2, 0,
                   0, sd_site_condition^2),
                 nrow = 2, byrow = TRUE)

# MODEL:
model <- simr::makeLmer(
  formula = response ~ condition +
    (1 + condition | site:participant) +
    (1 + condition | site),
  fixef = c(beta_intercept, beta_condition),
  VarCorr = list(V_parti, V_site),
  sigma = sd_residual,
  data = architecture
)


summary(model) # <- Your model in the population


# ----- b) Simulate data from the model using the chosen architecture ----

set.seed(0)
simulated_response <- simr::doSim(model)


# ----- c) Fit the model and test whether the estimand is statistically significant ----

fitted_model <- simr::doFit(simulated_response, model)

p_value <- simr::doTest(fitted_model, simr::fixed("condition", "t"))

p_value < 0.05
#> [1] TRUE


# ----- d) Repeat the process and calculate the proportion of null-hypothesis rejections ----

# (!) This step can take several minutes. If you prefer not to run it,
#     skip it and load the pre-computed results instead:
# load("Results/power_preliminary.RData")

power_preliminary <- simr::powerSim(model,
                                    test = simr::fixed("condition", "t"),
                                    nsim = 1000, alpha = 0.05,
                                    seed = 0, progress = FALSE)
power_preliminary # this is the power estimate for the preliminary architecture
#> Power for predictor 'condition', (95% confidence interval):
#>       55.80% (52.66, 58.91)
#>
#> Test: t-test with Satterthwaite degrees of freedom (package lmerTest)
#>       Effect size for condition is 0.15
#>
#> Based on 1000 simulations, (104 warnings, 0 errors)
#> alpha = 0.05, nrow = 6144
#>
#> Time elapsed: 0 h 8 m 5 s

# Tip 1: for a robust power estimate, use many simulations, such as the
#        1000 used here. If you just want a quick look in your own console,
#        you can use fewer (e.g., nsim = 200), but the estimate will be
#        much less precise: its 95% CI will be about twice as wide.

# Tip 2: if most fits are singular or {simr} reports many (convergence)
#        warnings, which is common with few sites, you can:
#        - Fit the model without correlations between random effects,
#          writing "||" instead of "|" in the formula, e.g.:
#          (1 + condition || site:participant) + (1 + condition || site)
#          Here the data are generated with zero covariances, so this
#          model matches the model that generates the data.
#        - Or even fit a Bayesian model (e.g., with {brms}), whose priors 
#          handle better variances close to zero, although simulating power
#          with it is much slower.
#        To see the warnings recorded by {simr}: power_preliminary$warnings


# How to read it: the power of this preliminary architecture is low
# (far below the usual 80%-90%). This is because our model is quite
# demanding: it assumes a lot of variability between participants and
# between sites in the effect of "condition" (sd_parti_condition and
# sd_site_condition), which makes the estimand harder to detect.
# Therefore, we need a larger architecture: more trials, participants
# and/or sites. But which combination? That is what we explore next.


# ---- 5. Define a feasible architecture space and optimize it with {mlpwr} --

# These ranges can be replaced by the values considered feasible for your study.
# Here, the minimum values are those of the preliminary architecture
# (Section 1), since we already know that its power is too low.
architecture_space <- list(trials_size = c(32, 48),
                           partis_size = c(16, 36),
                           sites_size = c(12, 20))

# {mlpwr} searches this space without evaluating every possible architecture.


# ----- 5.Step I. Create a simulation function ("simfun"): ---------------------

# {mlpwr} requires a "simfun" that receives an architecture and
# returns a TRUE/FALSE result indicating whether the null was rejected,
# with this structure:
#
# simfun <- function(trials_size, partis_size, sites_size) {
#                 1st: Generate a dataset
#                 2nd: Test the estimand
#               }
#
# *Note: the argument names of "simfun" must match the names used in
#        "architecture_space".
#
# *Note: "simfun_multisite" does not receive the model parameters as arguments.
#        It reads the CURRENT values of beta_intercept, beta_condition,
#        sd_parti_intercept, etc. from your workspace every time it runs. 
#        Therefore, if you change one of those values (as we do in Section 6), 
#        the simulations change accordingly.

simfun_multisite <- function(trials_size, partis_size, sites_size) {
  
  # Create the architecture (same as before):
  architecture <- expand.grid(trial = 1:trials_size,
                              participant = factor(1:partis_size),
                              site = factor(1:sites_size))
  
  condition_pattern <- rep(c(-0.5, 0.5), length.out = trials_size)
  architecture$condition <- rep(condition_pattern, times = partis_size * sites_size)
  
  # Specify the model parameters (values already defined in the workspace):
  params <- list(beta = c(beta_intercept, beta_condition),
                 theta = c(sd_parti_intercept / sd_residual,
                           0,
                           sd_parti_condition / sd_residual,
                           sd_site_intercept / sd_residual,
                           0,
                           sd_site_condition / sd_residual),
                 sigma = sd_residual)
  
  # Simulate one dataset:
  response <- lme4::simulate.formula( ~ condition +
                                        (1 + condition | site:participant) +
                                        (1 + condition | site),
                                      newdata = architecture,
                                      newparams = params)
  architecture$response <- response[[1]] # Include the simulated response to data
  
  # Fit the model to the simulated dataset:
  fitted_model <- lmerTest::lmer(response ~ condition +
                                   (1 + condition | site:participant) +
                                   (1 + condition | site),
                                 data = architecture)
  
  # Test the estimand:
  p_value <- coef(summary(fitted_model))["condition", "Pr(>|t|)"]
  
  return(p_value < 0.05) # TRUE if detected
}

# Check that one simulation returns TRUE or FALSE.
# simfun_multisite(trials_size = trials_size,
#                  partis_size = partis_size,
#                  sites_size = sites_size)


# ----- 5.Step II. Create a cost function ("costfun"): -------------------------

# {mlpwr} also requires a "costfun" that specifies the relative costs
# associated with each architecture level and returns the total cost.

# Example of costs:
trial_cost <- 1    # One additional trial costs 1 unit.
parti_cost <- 25   # One additional participant costs 25 times the trial cost.
site_cost  <- 1000 # One additional site costs 1000 times the trial cost.

costfun_multisite <- function(trials_size, partis_size, sites_size) {

  total_trials <- sites_size * partis_size * trials_size
  total_partis <- sites_size * partis_size
  total_sites  <- sites_size

  # Cost of all trials, all partis, and all sites:
  total_trial_cost <- total_trials * trial_cost
  total_parti_cost <- total_partis * parti_cost
  total_site_cost  <- total_sites * site_cost

  # Total cost of the architecture:
  total_cost <- total_trial_cost + total_parti_cost + total_site_cost
  
  # Return: total cost of the architecture:
  return(total_trial_cost + total_parti_cost + total_site_cost)
}


# ----- 5.Step III. Determine a desired power and optimize ---------------------

desired_power <- 0.80 # let's say 80% power

# *Note: with our demanding assumptions, 90% power would require more sites
#        than those we consider feasible (at most 20).

# (!) This step can take several minutes (each evaluation fits a mixed model).
#     If you prefer not to run it, skip it and load the pre-computed
#     results instead:
# load("Results/optimal_architecture.RData")

set.seed(0)
optimal_architecture <- mlpwr::find.design(
  boundaries = architecture_space,
  simfun = simfun_multisite,
  costfun = costfun_multisite,
  power = desired_power,
  evaluations = 2000
)

# *Note: Here, we use the default settings for the surrogate-model search. 
#        {mlpwr} provides additional arguments to control its initialization 
#        and updating procedure, which may be useful when the default search 
#        does not perform satisfactorily; see the package documentation or tutorial 
#        (Zimmer et al., 2024; https:/doi.org/10.3758/s13428-023-02269-0)

# (!) If you just want to see how find.design() works in your own console,
#     you can try a faster run with fewer evaluations (e.g., evaluations = 200).
#     But with so few evaluations, each architecture is tried with only a
#     handful of simulations, so the selected architecture and its power
#     estimate are imprecise and can change with the seed.
#     For a good and robust analysis, use many evaluations: the 2000 used
#     here, or even the package default (evaluations = 4000).

summary(optimal_architecture)
#> Call:
#> mlpwr::find.design(simfun = simfun_multisite, boundaries = architecture_space,
#>     power = desired_power, evaluations = 2000, costfun = costfun_multisite)
#>
#> Design: trials_size = 32, partis_size = 19, sites_size = 17
#>
#> Power: 0.80071,  SE: 0.01183,  Cost: 35411
#> Evaluations: 2006,  Time: 1626.67,  Updates: 47
#> Surrogate: Gaussian process regression

# Inspect the architectures evaluated by the surrogate model:
mlpwr::simulations_data(optimal_architecture)
# (First rows only: 32 architectures in total; some columns omitted)
#>    trials_size partis_size sites_size  cost power_surrogate     power Evaluations
#> 1           32          16         12 22944       0.6440735 0.5588235          34
#> 2           42          32         12 37728       0.6951855 0.6176471          34
#> 3           38          20         14 31640       0.7059500 0.6176471          34
#> 4           45          21         16 39520       0.7750736 0.6470588          34
#> 5           32          17         15 29535       0.7251521 0.6764706          34
#> 6           32          21         16 35152       0.7855429 0.6862745         102
#> ...

# *Note: {mlpwr} may warn that "the predicted value consistently lies at the
#        edge of the parameter space". This warning can appear even when the 
#        final architecture is well inside the space, so check the selected 
#        values above: only if one of them equals a boundary (e.g., sites_size = 20) 
#        may the true optimum lie outside your architecture space. In that case, 
#        consider widening that boundary (if feasible) or increasing "evaluations".
#        Also check that the "Power" line is close to the desired power: 
#        if it is far below it, no architecture within your space reaches the 
#        desired power (as happens in the scenarios of Section 6).
#        This can happen whether or not the selected values are at a boundary.
#
# *Note: here, trials_size = 32 equals the lower limit of its range.
#        {mlpwr} prefers to invest in participants and sites rather than in trials:
#        under our assumptions, what limits power is the variability of the effect
#        of "condition" between participants and sites, not the trial-level noise,
#        so adding trials helps little. An architecture with fewer trials might be
#        slightly cheaper, but we decided not to go below the preliminary one.


# ---- 6. Explore alternative assumptions --------------------------------------

# The optimal architecture may change under different
# assumptions about the effect size or between-site heterogeneity.
#
# (!) Each scenario runs find.design() again, so each one takes several minutes.

# Save the baseline assumptions before exploring alternatives.
baseline_beta_condition <- beta_condition
baseline_sd_site_condition <- sd_site_condition


# Scenario 1: smaller effect size (smallerES)

beta_condition <- 0.10 # instead of the baseline 0.15
sd_site_condition <- baseline_sd_site_condition

# (!) If you prefer not to run this step, skip it and load the
#     pre-computed results instead:
# load("Results/smallerES_architecture.RData")

set.seed(0)
smallerES_architecture <- mlpwr::find.design(
  boundaries = architecture_space,
  simfun = simfun_multisite,
  costfun = costfun_multisite,
  power = desired_power,
  evaluations = 2000
)

summary(smallerES_architecture)
#> Call:
#> mlpwr::find.design(simfun = simfun_multisite, boundaries = architecture_space,
#>     power = desired_power, evaluations = 2000, costfun = costfun_multisite)
#>
#> Design: trials_size = 48, partis_size = 36, sites_size = 18
#>
#> Power: 0.58825,  SE: 0.0134,  Cost: 65304
#> Evaluations: 2006,  Time: 3386.07,  Updates: 47
#> Surrogate: Gaussian process regression


# Scenario 2: greater between-site heterogeneity (greaterBSH)

beta_condition <- baseline_beta_condition
sd_site_condition <- 0.25 # instead of the baseline 0.15

# (!) If you prefer not to run this step, skip it and load the
#     pre-computed results instead:
# load("Results/greaterBSH_architecture.RData")

set.seed(0)
greaterBSH_architecture <- mlpwr::find.design(
  boundaries = architecture_space,
  simfun = simfun_multisite,
  costfun = costfun_multisite,
  power = desired_power,
  evaluations = 2000
)

summary(greaterBSH_architecture)
#> Call:
#> mlpwr::find.design(simfun = simfun_multisite, boundaries = architecture_space,
#>     power = desired_power, evaluations = 2000, costfun = costfun_multisite)
#>
#> Design: trials_size = 48, partis_size = 36, sites_size = 18
#>
#> Power: 0.63044,  SE: 0.02911,  Cost: 65304
#> Evaluations: 2006,  Time: 2729.56,  Updates: 47
#> Surrogate: Gaussian process regression


# Compare the three scenarios side by side:
scenario_comparison <- transform(
  rbind(baseline   = optimal_architecture$final$design,
        smallerES  = smallerES_architecture$final$design,
        greaterBSH = greaterBSH_architecture$final$design),
  total_cost = costfun_multisite(trials_size, partis_size, sites_size),
  power = round(c(optimal_architecture$final$power,
                  smallerES_architecture$final$power,
                  greaterBSH_architecture$final$power), 2))


scenario_comparison
#>            trials_size partis_size sites_size total_cost power
#> baseline            32          19         17      35411  0.80
#> smallerES           48          36         18      65304  0.59
#> greaterBSH          48          36         18      65304  0.63

# How to read it:
#   - First, check the "power" column: in both alternative scenarios, the
#     estimated power (0.59 and 0.63) is far below 80%, even though {mlpwr}
#     pushes the architecture to (almost) the upper limits of the space
#     (48 trials, 36 participants and 18 of the 20 sites), and warns that
#     its prediction lies at the edge of the space. Under these more
#     demanding assumptions, no architecture within our feasible space
#     reaches the desired power, so the architectures reported for these
#     scenarios are not valid solutions.
#   - Why? A smaller effect requires more information overall. And greater
#     between-site heterogeneity in the effect adds variability that can
#     only be averaged out by sampling more SITES: adding trials or
#     participants within the same sites helps much less. Since our space
#     allows at most 20 sites, even its largest architecture falls short.
#   - To reach 80% power under these assumptions, we would need to widen
#     the architecture space (especially sites_size) or to accept a lower
#     power.


# Restore the baseline assumptions:
beta_condition <- baseline_beta_condition
sd_site_condition <- baseline_sd_site_condition


# ---- 7. Assess the selected architecture -------------------------------------

# Extract the architecture selected by {mlpwr}.
# *Note: we reuse the names of Section 1, so the preliminary architecture
#        (and its model) is overwritten by the selected one. If you want to
#        rerun Section 4 afterwards, run Section 1 again first.
trials_size <- optimal_architecture$final$design[["trials_size"]]
partis_size <- optimal_architecture$final$design[["partis_size"]]
sites_size  <- optimal_architecture$final$design[["sites_size"]]


# Create the selected architecture (as in Section 1):
architecture <- expand.grid(trial = 1:trials_size,
                            participant = factor(1:partis_size),
                            site = factor(1:sites_size))

condition_pattern <- rep(c(-0.5, 0.5), length.out = trials_size)

architecture$condition <- rep(condition_pattern,
                              times = partis_size * sites_size)


# ----- 7a. Robust estimation of power (and precision) -------------------------
# {mlpwr} has suggested an optimal architecture, but the power estimate for that
# specific architecture may still be based on relatively few simulations
# (especially with few evaluations).
#
# We therefore run a larger number of simulations for the selected architecture.
# From the same simulated studies, we can assess both:
#   - Power: how often the null hypothesis is rejected.
#   - Precision: how variable and uncertain the estimated effect is.

# Redefine the model using the selected architecture.
selected_model <- simr::makeLmer(
  formula = response ~ condition +
    (1 + condition | site:participant) +
    (1 + condition | site),
  fixef = c(beta_intercept, beta_condition),
  VarCorr = list(V_parti,
                 V_site),
  sigma = sd_residual,
  data = architecture
)

# Number of simulated studies used for the final assessment.
nsim_final <- 1000

# Create empty vectors to store the result of each simulated study.
# If a simulated study fails (e.g., the model cannot be fitted), its values
# stay NA and that study is excluded from the summaries below.
p_values  <- rep(NA, nsim_final)
estimates <- rep(NA, nsim_final)
ses       <- rep(NA, nsim_final)
dfs       <- rep(NA, nsim_final)

set.seed(0)

# (!) This loop fits 1000 mixed models, so it can take several minutes.
#     If you prefer not to run it, skip it and load the pre-computed
#     results instead
#     (p_values, estimates, ses and dfs):
# load("Results/final_simulations.RData")

for (i in 1:nsim_final) {
  #Progress bar
  cat("\rFitting simulation ", i, " of ", nsim_final, "\n")  #Progress bar
 
   tryCatch({

    # Simulate one dataset from the selected model and architecture:
    simulated_response <- simr::doSim(selected_model)

    # Fit the analysis model to the simulated dataset:
    fitted_model <- simr::doFit(simulated_response, selected_model)

    # Obtain the t-test with Satterthwaite df 
    # (the same test used by simr::fixed("condition", "t")):
    coefficient_table <- coef(summary(
      lmerTest::as_lmerModLmerTest(fitted_model)
    ))

    # IV. Store the p-value, the estimated effect, its standard error and df:
    p_values[i]  <- coefficient_table["condition", "Pr(>|t|)"]
    estimates[i] <- coefficient_table["condition", "Estimate"]
    ses[i]       <- coefficient_table["condition", "Std. Error"]
    dfs[i]       <- coefficient_table["condition", "df"]

  }, error = function(e) NULL) # a failed study simply stays NA
}


# POWER: proportion of simulated studies in which the null hypothesis was rejected.
power_final <- mean(p_values < 0.05, na.rm = TRUE)
power_final
#> [1] 0.763

# *Note: this is slightly lower than the 0.80 estimated by {mlpwr}: even
#        with 2000 evaluations, its estimate for the selected architecture
#        was based on only 170 simulations (Section 5).


# PRECISION: how much uncertainty there is around the estimated effect.
# We summarize it using:
#   - empirical_se: how much the estimated effect varies across simulated studies.
#   - mean_model_se: the average standard error reported by the fitted models.
#   - mean_95CI_width: the average width of the 95% confidence interval.
#   - empirical_coverage: the proportion of 95% confidence intervals
#     containing the true effect used to generate the data.
#
# *Note: the 95% CIs use the t distribution with Satterthwaite df, 
#        i.e., the same distribution used by the significance test. 
#        With few sites, these df are small, and using 1.96 instead 
#        would give CIs that are too narrow.

t_critical <- qt(0.975, df = dfs)
ci_lower   <- estimates - t_critical * ses
ci_upper   <- estimates + t_critical * ses

precision_results <- data.frame(
  mean_estimate = mean(estimates, na.rm = TRUE),
  empirical_se = sd(estimates, na.rm = TRUE),
  mean_model_se = mean(ses, na.rm = TRUE),
  mean_95CI_width = mean(ci_upper - ci_lower, na.rm = TRUE),
  empirical_coverage = mean(ci_lower <= beta_condition &
                              ci_upper >= beta_condition,
                            na.rm = TRUE),
  successful_simulations = sum(!is.na(p_values))
)

precision_results
#>   mean_estimate empirical_se mean_model_se mean_95CI_width empirical_coverage
#> 1     0.1484726   0.05359122    0.05116042       0.2162517              0.948
#>   successful_simulations
#> 1                   1000


# ----- 7b. Minimum detectable effect size (MDES) ------------------------------
# The MDES is the smallest effect size that the selected architecture can detect
# with the desired power.

effect_step <- 0.01
nsim_mdes <- 200 # (!) More simulations should be used for the final MDES.

current_effect <- beta_condition
current_power <- power_final

# *Note: the starting power comes from the 1000 simulations of Section 7a,
#        while the following steps use "nsim_mdes" simulations each, 
#        so the later power estimates are less precise.

mdes_results <- data.frame(effect_size = current_effect,
                          estimated_power = current_power)

# (!) The search below can take several minutes. If you prefer not to run it,
#     skip it and load the pre-computed results instead:
# load("Results/mdes_results.RData")
#
#     (If you load them, skip only the "repeat { ... }" loop below:
#     the other lines are fast and are needed to obtain the MDES.)

# Direction of the search:
#   - If the selected architecture REACHES the desired power for the assumed
#     effect, we DECREASE the effect until the power drops below it.
#   - If it does NOT reach it (as can happen when {mlpwr} used few
#     evaluations), we INCREASE the effect until the desired power is reached.
search_down <- current_power >= desired_power
direction   <- if (search_down) -1 else 1

max_iter <- 50  # Maximum number of effect sizes to test
iter <- 0       # Counter for the number of iterations

repeat {
  
  iter <- iter + 1
  
  # Move to the next candidate effect size:
  current_effect <- current_effect + direction * effect_step
  
  # Build the model using the current candidate effect size:
  mdes_model <- simr::makeLmer(
    response ~ condition +
      (1 + condition | site:participant) +
      (1 + condition | site),
    fixef = c(beta_intercept, current_effect),
    VarCorr = list(V_parti,
                   V_site),
    sigma = sd_residual,
    data = architecture
  )
  
  # Estimate power for the current effect size:
  mdes_result <- simr::powerSim(
    mdes_model,
    test = simr::fixed("condition", "t"),
    nsim = nsim_mdes,
    alpha = 0.05,
    seed = 0,
    progress = FALSE
  )
  
  # Compute the proportion of simulations in which the effect was detected:
  current_power <- mdes_result$x / mdes_result$n
  
  # Store the effect size tested and its estimated power:
  mdes_results <- rbind(
    mdes_results,
    data.frame(effect_size = current_effect,
               estimated_power = current_power)
  )
  
  # Stop once the search crosses the desired power.
  if (search_down && current_power < desired_power) break
  if (!search_down && current_power >= desired_power) break
  
  # Stop if the maximum number of iterations is reached:
  if (iter >= max_iter) {
    warning(
      "Maximum number of iterations reached. ",
      "Last effect tested: ", round(current_effect, 3),
      "; estimated power: ", round(current_power, 3), "."
    )
    break
  }
}
mdes_results
#>   effect_size estimated_power
#> 1        0.15           0.763
#> 2        0.16           0.835

# The MDES is the smallest effect size that still reaches the desired power:
#   - Searching down: the last effect fell below it, so the MDES is the
#     previous one.
#   - Searching up: the last effect is the first one that reaches it.
mdes <- if (search_down) {
  mdes_results$effect_size[nrow(mdes_results) - 1]
} else {
  mdes_results$effect_size[nrow(mdes_results)]
}

mdes
#> [1] 0.16

# How to read it: with the architecture selected by {mlpwr}, we would have 80%
# power only for effects of 0.16 or larger, not for the assumed 0.15. This is
# why the selected architecture should always be checked with more simulations
# (or found with more evaluations, as recommended in Section 5).


# ============================================================================ # 
# Congrats! You now know how to select your multisite study architecture :)
# ============================================================================ # 
