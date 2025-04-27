# =============================================
# 1. Load Required Libraries
# =============================================
if (!require("dplyr")) install.packages("dplyr")
if (!require("lubridate")) install.packages("lubridate")
if (!require("quadprog")) install.packages("quadprog")
if (!require("tidyr")) install.packages("tidyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("Matrix")) install.packages("Matrix")
if (!require("zoo")) install.packages("zoo")
library(dplyr)
library(lubridate)
library(quadprog)
library(tidyr)
library(ggplot2)
library(Matrix)
library(zoo)

# =============================================
# 2. Load Data
# =============================================

# Industry returns
industry_raw <- read.csv("./49_Industry_Portfolios.CSV")
industry <- industry_raw %>%
  mutate(Date = as.Date(paste0(Date, "01"), format = "%Y%m%d")) %>%
  mutate(across(-Date, ~ as.numeric(.) / 100))

# Fama-French factors
factors <- read.csv("./F-F_Research_Data_Factors.CSV") %>%
  mutate(Date = as.Date(paste0(X, "01"), format = "%Y%m%d")) %>%
  rename(RMRF = Mkt.RF) %>%
  select(Date, RMRF, RF) %>%
  mutate(across(c(RMRF, RF), ~ as.numeric(.) / 100))

# Macro uncertainty
macroU <- read.csv("./sydneymacrouncertainty.csv") %>%
  mutate(Date = parse_date_time(Date, orders = "m/Y")) %>%
  rename(MacroU = h.1) %>%
  select(Date, MacroU)

# 1-year treasury rate
treasury <- read.csv("./THREEFY1.csv") %>%
  rename(Date = observation_date, ShortRate = THREEFY1) %>%
  mutate(Date = as.Date(Date), ShortRate = as.numeric(ShortRate)) %>%
  select(Date, ShortRate)

# Trade policy uncertainty
tpu <- read.csv("./TPU INDEX MONTHLY.csv") %>%
  mutate(Date = as.Date(DATE, origin = "1899-12-30")) %>%
  select(Date, TPU)

# Economic policy uncertainty
epu <- read.csv("./EPU.csv") %>%
  mutate(Date = as.Date(paste(Year, Month, "01", sep = "-"))) %>%
  rename(EPU = News_Based_Policy_Uncert_Index) %>%
  select(Date, EPU)

# Merge all
df_factors <- factors %>%
  inner_join(macroU, by = "Date") %>%
  inner_join(tpu, by = "Date") %>%
  inner_join(treasury, by = "Date") %>%
  inner_join(epu, by = "Date")

full_data_monthly <- industry %>%
  inner_join(df_factors, by = "Date")

# =============================================
# 3. Preprocessing: Normalize MacroU, TPU, EPU
# =============================================
full_data_monthly <- full_data_monthly %>%
  mutate(
    MacroU_z = scale(MacroU)[,1],
    TPU_z = scale(TPU)[,1],
    EPU_z = scale(EPU)[,1]
  )

# =============================================
# 4. Define Industry Columns
# =============================================
industry_cols <- setdiff(names(full_data_monthly),
                         c("Date", "RMRF", "RF", "MacroU", "TPU", "ShortRate", "EPU", "MacroU_z", "TPU_z", "EPU_z"))

# =============================================
# 5. Optimization Function
# =============================================
optimize_multi_factor <- function(Sigma, beta_market, beta_macrou, beta_tpu, beta_epu,
                                  target_market_beta = 1,
                                  target_macrou = 0.0,
                                  target_tpu = 0.0,
                                  target_epu = 0.0) {
  n <- length(beta_market)
  Sigma_reg <- Sigma + diag(1e-4, n)
  
  Dmat <- Sigma_reg
  dvec <- rep(0, n)
  
  # Now allow short-selling
  Amat <- cbind(rep(1, n), beta_market, beta_macrou, beta_tpu, beta_epu)
  bvec <- c(1, target_market_beta, target_macrou, target_tpu, target_epu)
  meq <- 5  # all are equality
  
  result <- tryCatch({
    sol <- solve.QP(Dmat, dvec, Amat, bvec, meq)
    w <- sol$solution
    # ⚡ NO pmax, pmin anymore: allow negative weights (shorts)
    w/sum(w)  # normalize if needed
  }, error = function(e) {
    rep(1/n, n)
  })
  
  return(result)
}


# =============================================
# 6. Backtest Implementation
# =============================================
lookback_months <- 24
rebalance_every <- 3
portfolio_returns <- rep(NA, nrow(full_data_monthly))
portfolio_weights <- matrix(NA, nrow = nrow(full_data_monthly),
                            ncol = length(industry_cols))
colnames(portfolio_weights) <- industry_cols

for (t in (lookback_months + 1):nrow(full_data_monthly)) {
  if ((t - lookback_months - 1) %% rebalance_every != 0 && t > (lookback_months + 1)) {
    portfolio_weights[t, ] <- portfolio_weights[t - 1, ]
    portfolio_returns[t] <- sum(portfolio_weights[t, ] * full_data_monthly[t, industry_cols], na.rm = TRUE)
    next
  }
  
  df_window <- full_data_monthly[(t - lookback_months):(t - 1), ]
  
  valid_returns <- df_window[, industry_cols][complete.cases(df_window[, industry_cols]), ]
  if (nrow(valid_returns) < 12) {
    portfolio_weights[t, ] <- rep(1/length(industry_cols), length(industry_cols))
    portfolio_returns[t] <- mean(as.numeric(full_data_monthly[t, industry_cols]), na.rm = TRUE)
    next
  }
  
  Sigma_est <- cov(valid_returns) * 12
  Sigma_est <- Sigma_est + diag(1e-4, ncol(Sigma_est))
  
  beta_market <- sapply(industry_cols, function(col) {
    fit <- lm(df_window[[col]] ~ df_window$RMRF)
    0.5 * coef(fit)[2] + 0.5 * 1  # shrinkage to 1
  })
  
  beta_macrou <- sapply(industry_cols, function(col) {
    fit <- lm(df_window[[col]] ~ df_window$MacroU_z)
    coef(fit)[2]
  })
  
  beta_tpu <- sapply(industry_cols, function(col) {
    fit <- lm(df_window[[col]] ~ df_window$TPU_z)
    coef(fit)[2]
  })
  
  beta_epu <- sapply(industry_cols, function(col) {
    fit <- lm(df_window[[col]] ~ df_window$EPU_z)
    coef(fit)[2]
  })
  
  w_opt <- optimize_multi_factor(Sigma_est, beta_market, beta_macrou, beta_tpu, beta_epu)
  
  portfolio_weights[t, ] <- w_opt
  portfolio_returns[t] <- sum(w_opt * full_data_monthly[t, industry_cols], na.rm = TRUE)
}
names(full_data_monthly)
# =============================================
# 7. Performance Analysis
# =============================================
results <- data.frame(
  Date = full_data_monthly$Date,
  Portfolio_Return = portfolio_returns,
  Market_Return = full_data_monthly$RMRF + full_data_monthly$RF,
  RF = full_data_monthly$RF
) %>% na.omit() %>%
  mutate(
    Portfolio_Cumulative = cumprod(1 + Portfolio_Return),
    Market_Cumulative = cumprod(1 + Market_Return),
    Portfolio_Excess = Portfolio_Return - RF,
    Market_Excess = Market_Return - RF
  )

cat("Portfolio:\n",
    "Return =", round(mean(results$Portfolio_Excess) * 12, 4),
    "Vol =", round(sd(results$Portfolio_Return) * sqrt(12), 4),
    "Sharpe =", round(mean(results$Portfolio_Excess)/sd(results$Portfolio_Return) * sqrt(12), 2), "\n")

cat("Market:\n",
    "Return =", round(mean(results$Market_Excess) * 12, 4),
    "Vol =", round(sd(results$Market_Return) * sqrt(12), 4),
    "Sharpe =", round(mean(results$Market_Excess)/sd(results$Market_Return) * sqrt(12), 2), "\n")

# =============================================
# 8. Plot Cumulative Return
# =============================================
ggplot(results, aes(x = Date)) +
  geom_line(aes(y = Portfolio_Cumulative, color = "Optimized"), linewidth = 1) +
  geom_line(aes(y = Market_Cumulative, color = "Market"), linewidth = 1) +
  scale_color_manual(values = c("Optimized" = "blue", "Market" = "red")) +
  labs(title = "Monthly Cumulative Returns", y = "Growth of $1") +
  theme_minimal()

results

# =====================================
# 1. Pick 3 random industries
# =====================================
set.seed(123)  # for reproducibility
selected_industries <- sample(industry_cols, 3)
print(selected_industries)

# =====================================
# 2. Regressions for Portfolio, Market, and Selected Industries
# =====================================

# Factors you used
factor_vars <- c("RMRF", "MacroU_z", "TPU_z", "EPU_z")

# Prepare data subset
data_regression <- full_data_monthly %>%
  select(Date, all_of(factor_vars)) %>%
  inner_join(
    results %>% select(Date, Portfolio_Excess, Market_Excess),
    by = "Date"
  ) %>%
  inner_join(
    full_data_monthly %>% select(Date, all_of(selected_industries)),
    by = "Date"
  ) %>% na.omit()

# Function to regress and extract coefficients
get_betas <- function(dep_var, data) {
  formula <- as.formula(paste(dep_var, "~", paste(factor_vars, collapse = " + ")))
  fit <- lm(formula, data = data)
  coefs <- coef(fit)[-1]  # remove intercept
  return(coefs)
}

# Calculate betas
beta_portfolio <- get_betas("Portfolio_Excess", data_regression)
beta_market <- get_betas("Market_Excess", data_regression)

beta_industries <- sapply(selected_industries, function(ind) {
  get_betas(ind, data_regression)
})

# =====================================
# 3. Format output nicely
# =====================================
beta_compare <- data.frame(
  Factor = factor_vars,
  Portfolio_Beta = round(beta_portfolio, 5),
  Market_Beta = round(beta_market, 5)
)

for (i in seq_along(selected_industries)) {
  beta_compare[[selected_industries[i]]] <- round(beta_industries[, i], 5)
}

print(beta_compare)

# =============================================
# 1. Subset the data (2000–2024)
# =============================================
data_subset <- full_data_monthly %>%
  filter(Date >= as.Date("2000-01-01"), Date <= as.Date("2024-12-31"))

# =============================================
# 2. Identify Industry Columns
# =============================================
# Only keep the columns you know are industry returns
# Your industries are exactly the first 49 columns in the original 49-industry portfolio CSV
# So better select using the known pattern:
industry_cols <- setdiff(names(industry), "Date")

# =============================================
# 4. Loop: Regress and store betas
# =============================================
industry_beta_results <- data.frame(
  Industry = character(),
  RMRF = numeric(),
  MacroU_z = numeric(),
  TPU_z = numeric(),
  EPU_z = numeric(),
  stringsAsFactors = FALSE
)

for (ind in industry_cols) {
  model <- lm(as.formula(paste(ind, "~ RMRF + MacroU_z + TPU_z + EPU_z")), data = data_subset)
  
  betas <- coef(model)[-1]  # Exclude intercept
  
  industry_beta_results <- rbind(industry_beta_results, 
                                 data.frame(
                                   Industry = ind,
                                   RMRF = betas["RMRF"],
                                   MacroU_z = betas["MacroU_z"],
                                   TPU_z = betas["TPU_z"],
                                   EPU_z = betas["EPU_z"]
                                 ))
}

# =============================================
# 5. Calculate Mean Absolute Beta for each Industry
# =============================================
industry_beta_results <- industry_beta_results %>%
  mutate(
    Mean_Abs_Beta = rowMeans(abs(select(., RMRF, MacroU_z, TPU_z, EPU_z)))
  ) %>%
  arrange(desc(Mean_Abs_Beta))  # Rank descending

# =============================================
# 6. Show the Ranked Table
# =============================================
print(industry_beta_results)


# =============================================

# =============================================
# =============================================
# 1. Define a Simpler Optimization Function
# =============================================
optimize_variance_target_market <- function(Sigma, beta_market, 
                                            target_market_beta = 1.0) {
  n <- length(beta_market)
  Sigma_reg <- Sigma + diag(1e-4, n)
  
  Dmat <- Sigma_reg
  dvec <- rep(0, n)
  
  Amat <- cbind(rep(1, n), beta_market)
  bvec <- c(1, target_market_beta)
  meq <- 2
  
  result <- tryCatch({
    sol <- solve.QP(Dmat, dvec, Amat, bvec, meq)
    w <- sol$solution
    w <- pmax(w, 0)
    w <- pmin(w, 0.2)  # max 20% per sector
    w / sum(w)         # re-normalize
  }, error = function(e) {
    rep(1/n, n)
  })
  
  return(result)
}

# =============================================
# 2. Backtest with Only Variance Minimization
# =============================================

lookback_months <- 24
rebalance_every <- 3
portfolio_returns_varonly <- rep(NA, nrow(full_data_monthly))
portfolio_weights_varonly <- matrix(NA, nrow = nrow(full_data_monthly),
                                    ncol = length(industry_cols))
colnames(portfolio_weights_varonly) <- industry_cols

for (t in (lookback_months + 1):nrow(full_data_monthly)) {
  if ((t - lookback_months - 1) %% rebalance_every != 0 && t > (lookback_months + 1)) {
    portfolio_weights_varonly[t, ] <- portfolio_weights_varonly[t - 1, ]
    portfolio_returns_varonly[t] <- sum(portfolio_weights_varonly[t, ] * full_data_monthly[t, industry_cols], na.rm = TRUE)
    next
  }
  
  df_window <- full_data_monthly[(t - lookback_months):(t - 1), ]
  valid_returns <- df_window[, industry_cols][complete.cases(df_window[, industry_cols]), ]
  
  if (nrow(valid_returns) < 12) {
    portfolio_weights_varonly[t, ] <- rep(1/length(industry_cols), length(industry_cols))
    portfolio_returns_varonly[t] <- mean(as.numeric(full_data_monthly[t, industry_cols]), na.rm = TRUE)
    next
  }
  
  Sigma_est <- cov(valid_returns) * 12
  Sigma_est <- Sigma_est + diag(1e-4, ncol(Sigma_est))
  
  beta_market <- sapply(industry_cols, function(col) {
    fit <- lm(df_window[[col]] ~ df_window$RMRF)
    0.5 * coef(fit)[2] + 0.5 * 1  # shrinkage toward 1
  })
  
  w_opt <- optimize_variance_target_market(Sigma_est, beta_market)
  
  portfolio_weights_varonly[t, ] <- w_opt
  portfolio_returns_varonly[t] <- sum(w_opt * full_data_monthly[t, industry_cols], na.rm = TRUE)
}

# =============================================
# 3. Performance Analysis
# =============================================

results_varonly <- data.frame(
  Date = full_data_monthly$Date,
  Portfolio_Return = portfolio_returns_varonly,
  Market_Return = full_data_monthly$RMRF + full_data_monthly$RF,
  RF = full_data_monthly$RF
) %>% na.omit() %>%
  mutate(
    Portfolio_Cumulative = cumprod(1 + Portfolio_Return),
    Market_Cumulative = cumprod(1 + Market_Return),
    Portfolio_Excess = Portfolio_Return - RF,
    Market_Excess = Market_Return - RF
  )

cat("Variance-Minimized Portfolio:\n",
    "Return =", round(mean(results_varonly$Portfolio_Excess) * 12, 4),
    "Vol =", round(sd(results_varonly$Portfolio_Return) * sqrt(12), 4),
    "Sharpe =", round(mean(results_varonly$Portfolio_Excess)/sd(results_varonly$Portfolio_Return) * sqrt(12), 2), "\n")

cat("Market:\n",
    "Return =", round(mean(results_varonly$Market_Excess) * 12, 4),
    "Vol =", round(sd(results_varonly$Market_Return) * sqrt(12), 4),
    "Sharpe =", round(mean(results_varonly$Market_Excess)/sd(results_varonly$Market_Return) * sqrt(12), 2), "\n")

# =============================================
# 4. Plot Cumulative Returns
# =============================================
ggplot(results_varonly, aes(x = Date)) +
  geom_line(aes(y = Portfolio_Cumulative, color = "Variance-Minimized Portfolio"), linewidth = 1) +
  geom_line(aes(y = Market_Cumulative, color = "Market"), linewidth = 1) +
  scale_color_manual(values = c("Variance-Minimized Portfolio" = "blue", "Market" = "red")) +
  labs(title = "Monthly Cumulative Returns (Variance-Minimized Strategy)", y = "Growth of $1") +
  theme_minimal()

