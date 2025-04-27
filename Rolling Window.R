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

# Financial Uncertainty
finU<- read.csv('./FINU.csv')
finU <- finU %>%
  mutate(Date = parse_date_time(Date, orders = "m/Y")) %>%  # Parse month-year correctly
  mutate(Date = as.Date(paste0(year(Date), "-", sprintf("%02d", month(Date)), "-01")))  # Force first day of month

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
  inner_join(epu, by = "Date") %>%
  inner_join(finU, by = "Date")

full_data_monthly <- industry %>%
  inner_join(df_factors, by = "Date")

names(full_data_monthly)
names(df_factors)


# =============================================
# 1. Define rolling window settings
# =============================================
lookback_months <- 24  # 2 years window

# Industry columns
industry_cols <- setdiff(names(full_data_monthly),
                         c("Date", "RMRF", "RF", "MacroU", "TPU", "ShortRate", "EPU", "FinU"))

# Prepare storage for betas
rolling_betas <- list()

# =============================================
# 2. Rolling Regressions
# =============================================
for (i in (lookback_months + 1):nrow(full_data_monthly)) {
  
  window_data <- full_data_monthly[(i - lookback_months):(i - 1), ]
  
  betas_this_date <- data.frame(
    Date = full_data_monthly$Date[i],
    Industry = industry_cols,
    Beta_RMRF = NA,
    Beta_MacroU = NA,
    Beta_FinU = NA,
    Beta_TPU = NA,
    Beta_EPU = NA,
    Beta_ShortRate = NA
  )
  
  for (j in seq_along(industry_cols)) {
    ind <- industry_cols[j]
    
    # Excess return
    y <- window_data[[ind]] - window_data$RF
    x <- window_data %>%
      select(RMRF, MacroU, FinU, TPU, EPU, ShortRate)
    
    model <- lm(y ~ ., data = x)
    
    coefs <- coef(model)
    
    # Store betas
    betas_this_date$Beta_RMRF[j] <- coefs["RMRF"]
    betas_this_date$Beta_MacroU[j] <- coefs["MacroU"]
    betas_this_date$Beta_FinU[j] <- coefs["FinU"]
    betas_this_date$Beta_TPU[j] <- coefs["TPU"]
    betas_this_date$Beta_EPU[j] <- coefs["EPU"]
    betas_this_date$Beta_ShortRate[j] <- coefs["ShortRate"]
  }
  
  rolling_betas[[i - lookback_months]] <- betas_this_date
}

# =============================================
# 3. Combine into one big dataframe
# =============================================
beta_estimates <- do.call(rbind, rolling_betas)

# =============================================
# 4. Check result
# =============================================
head(beta_estimates)

# =============================================
# 1. Define optimization function
# =============================================
library(quadprog)

optimize_portfolio <- function(Sigma, betas, 
                               target_market_beta = 1, 
                               target_macroU = 0, 
                               target_finU = 0,
                               target_TPU = 0,
                               target_EPU = 0,
                               target_ShortRate = 0) {
  
  n <- nrow(betas)
  
  # Add a small regularization to covariance matrix
  Sigma_reg <- Sigma + diag(1e-6, n)
  
  Dmat <- Sigma_reg
  dvec <- rep(0, n)
  
  # Constraints: full investment and factor exposure targets
  Amat <- t(rbind(
    rep(1, n),
    betas$Beta_RMRF,
    betas$Beta_MacroU,
    betas$Beta_FinU,
    betas$Beta_TPU,
    betas$Beta_EPU,
    betas$Beta_ShortRate
  ))
  
  bvec <- c(1, target_market_beta, target_macroU, target_finU, target_TPU, target_EPU, target_ShortRate)
  
  meq <- length(bvec)  # all are equality constraints
  
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq)
  
  w_opt <- sol$solution
  w_opt <- w_opt / sum(w_opt)  # Normalize
  
  return(w_opt)
}

# =============================================
# 2. Setup dynamic backtest
# =============================================

# Parameters
lookback_months <- 24
rebalance_every <- 6  # rebalance every 6 months
n_assets <- length(industry_cols)

# Initialize storage
portfolio_returns <- rep(NA, nrow(full_data_monthly))
portfolio_weights <- matrix(NA, nrow = nrow(full_data_monthly), ncol = n_assets)
colnames(portfolio_weights) <- industry_cols

# Start loop
for (t in (lookback_months + 1):nrow(full_data_monthly)) {
  
  if ((t - lookback_months - 1) %% rebalance_every != 0 && t > (lookback_months + 1)) {
    # If not a rebalance month, use last weights
    portfolio_weights[t, ] <- portfolio_weights[t - 1, ]
    portfolio_returns[t] <- sum(portfolio_weights[t, ] * full_data_monthly[t, industry_cols], na.rm = TRUE)
    next
  }
  
  # Get the correct betas
  betas_this_date <- beta_estimates %>%
    filter(Date == full_data_monthly$Date[t]) %>%
    arrange(match(Industry, industry_cols))
  
  if (nrow(betas_this_date) != n_assets) {
    portfolio_weights[t, ] <- rep(1/n_assets, n_assets)
    portfolio_returns[t] <- mean(as.numeric(full_data_monthly[t, industry_cols]), na.rm = TRUE)
    next
  }
  
  # Estimate Sigma (Covariance) from past 24 months
  window_data <- full_data_monthly[(t - lookback_months):(t - 1), industry_cols]
  Sigma_est <- cov(window_data) * 12
  Sigma_est <- Sigma_est + diag(1e-6, ncol(Sigma_est))
  
  # Solve optimization
  w_opt <- optimize_portfolio(Sigma_est, betas_this_date,
                              target_market_beta = 0.7,   # 🛡️ Defensive tilt
                              target_macroU = 0,
                              target_finU = 0,
                              target_TPU = 0,
                              target_EPU = 0,
                              target_ShortRate = 0)
  
  portfolio_weights[t, ] <- w_opt
  portfolio_returns[t] <- sum(w_opt * full_data_monthly[t, industry_cols], na.rm = TRUE)
}

# =============================================
# 3. Build final results
# =================================
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

# =============================================
# 4. Quick plot
# =============================================
library(ggplot2)

ggplot(results, aes(x = Date)) +
  geom_line(aes(y = Portfolio_Cumulative, color = "Optimized Portfolio"), linewidth = 1) +
  geom_line(aes(y = Market_Cumulative, color = "Market"), linewidth = 1) +
  scale_color_manual(values = c("Optimized Portfolio" = "blue", "Market" = "red")) +
  labs(title = "Dynamic Portfolio vs Market Cumulative Return",
       y = "Growth of $1",
       x = "Date") +
  theme_minimal()

names(results)

# Portfolio performance
portfolio_annual_return <- mean(results$Portfolio_Return) * 12
portfolio_annual_vol <- sd(results$Portfolio_Return) * sqrt(12)
portfolio_sharpe <- mean(results$Portfolio_Excess) / sd(results$Portfolio_Return) * sqrt(12)

# Market performance
market_annual_return <- mean(results$Market_Return) * 12
market_annual_vol <- sd(results$Market_Return) * sqrt(12)
market_sharpe <- mean(results$Market_Excess) / sd(results$Market_Return) * sqrt(12)

# Display
cat("Optimized Portfolio:\n",
    "Annual Return =", round(portfolio_annual_return, 4),
    "Annual Vol =", round(portfolio_annual_vol, 4),
    "Sharpe =", round(portfolio_sharpe, 2), "\n")

cat("Market:\n",
    "Annual Return =", round(market_annual_return, 4),
    "Annual Vol =", round(market_annual_vol, 4),
    "Sharpe =", round(market_sharpe, 2), "\n")

analyze_crisis_metrics_full <- function(results, crisis_periods) {
  crisis_summary <- list()
  
  for (crisis_name in names(crisis_periods)) {
    period <- crisis_periods[[crisis_name]]
    
    crisis_data <- results %>%
      filter(Date >= period[1], Date <= period[2])
    
    n_months <- nrow(crisis_data)  # because returns are monthly
    
    # Cumulative Return
    portfolio_cum_return <- prod(1 + crisis_data$Portfolio_Return) - 1
    market_cum_return <- prod(1 + crisis_data$Market_Return) - 1
    
    # Annualized Return
    portfolio_annual_return <- (1 + portfolio_cum_return)^(12 / n_months) - 1
    market_annual_return <- (1 + market_cum_return)^(12 / n_months) - 1
    
    # Annualized Volatility
    portfolio_volatility <- sd(crisis_data$Portfolio_Return) * sqrt(12)
    market_volatility <- sd(crisis_data$Market_Return) * sqrt(12)
    
    # Sharpe Ratio (using excess returns over risk-free rate)
    portfolio_sharpe <- mean(crisis_data$Portfolio_Excess) / sd(crisis_data$Portfolio_Excess) * sqrt(12)
    market_sharpe <- mean(crisis_data$Market_Excess) / sd(crisis_data$Market_Excess) * sqrt(12)
    
    # Maximum Drawdown (based on cumulative returns)
    max_drawdown <- function(cum_return_series) {
      cummax_series <- cummax(cum_return_series)
      drawdowns <- (cummax_series - cum_return_series) / cummax_series
      max(drawdowns, na.rm = TRUE)
    }
    
    portfolio_mdd <- max_drawdown(crisis_data$Portfolio_Cumulative)
    market_mdd <- max_drawdown(crisis_data$Market_Cumulative)
    
    crisis_result <- data.frame(
      Crisis = crisis_name,
      Type = c("Portfolio", "Market"),
      CumulativeReturn = c(portfolio_cum_return, market_cum_return),
      AnnualizedReturn = c(portfolio_annual_return, market_annual_return),
      SharpeRatio = c(portfolio_sharpe, market_sharpe),
      AnnualizedVolatility = c(portfolio_volatility, market_volatility),
      MaxDrawdown = c(portfolio_mdd, market_mdd)
    )
    
    crisis_summary[[crisis_name]] <- crisis_result
  }
  
  do.call(rbind, crisis_summary)
}

crisis_periods_zoomed_extended <- list(
  Housing_Crisis_2008 = c(as.Date("2008-09-01"), as.Date("2009-03-31")),
  Trump_Trade_War = c(as.Date("2018-10-01"), as.Date("2018-12-31")),
  Covid_19_Crisis = c(as.Date("2020-02-01"), as.Date("2020-03-31")),
  Trade_Escalation_May2019 = c(as.Date("2019-05-01"), as.Date("2019-06-30")),
  
  Asian_Financial_Crisis = c(as.Date("1997-07-01"), as.Date("1998-06-30")),
  Dotcom_Bubble_Burst = c(as.Date("2000-03-01"), as.Date("2002-10-31")),
  Taper_Tantrum = c(as.Date("2013-05-01"), as.Date("2013-12-31"))
)


crisis_metrics_full_zoomed <- analyze_crisis_metrics_full(results, crisis_periods_zoomed_extended)

# =============================================
# View Results
# =============================================
print(crisis_metrics_full_zoomed)


# =============================================
# 1. Calculate Dynamic Portfolio Factor Exposures
# =============================================

# Initialize storage
portfolio_beta_TPU <- rep(NA, nrow(full_data_monthly))
portfolio_beta_MacroU <- rep(NA, nrow(full_data_monthly))
portfolio_beta_FinU <- rep(NA, nrow(full_data_monthly))

# Match industry columns
industry_cols <- setdiff(names(full_data_monthly),
                         c("Date", "RMRF", "RF", "MacroU", "TPU", "ShortRate", "EPU", "FinU"))

for (t in (lookback_months + 1):nrow(full_data_monthly)) {
  
  if (is.na(portfolio_weights[t, 1])) next  # Skip missing
  
  betas_this_date <- beta_estimates %>%
    filter(Date == full_data_monthly$Date[t]) %>%
    arrange(match(Industry, industry_cols))
  
  if (nrow(betas_this_date) != length(industry_cols)) next
  
  w_t <- portfolio_weights[t, ]
  
  portfolio_beta_TPU[t] <- sum(w_t * betas_this_date$Beta_TPU)
  portfolio_beta_MacroU[t] <- sum(w_t * betas_this_date$Beta_MacroU)
  portfolio_beta_FinU[t] <- sum(w_t * betas_this_date$Beta_FinU)
}

# Combine into dataframe
dynamic_betas <- data.frame(
  Date = full_data_monthly$Date,
  TPU_Beta = portfolio_beta_TPU,
  MacroU_Beta = portfolio_beta_MacroU,
  FinU_Beta = portfolio_beta_FinU
)
# =============================================
# 2. Plot Dynamic TPU, MacroU, FinU Betas
# =============================================

library(ggplot2)

dynamic_betas_long <- dynamic_betas %>%
  pivot_longer(cols = c(TPU_Beta, MacroU_Beta, FinU_Beta),
               names_to = "Factor",
               values_to = "Beta")

ggplot(dynamic_betas_long, aes(x = Date, y = Beta, color = Factor)) +
  geom_line(linewidth = 1) +
  labs(title = "Dynamic Portfolio Factor Betas Over Time",
       y = "Beta Exposure",
       x = "Date") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  scale_color_manual(values = c("TPU_Beta" = "red", "MacroU_Beta" = "blue", "FinU_Beta" = "green"))

