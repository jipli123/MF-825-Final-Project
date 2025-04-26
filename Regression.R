# =============================================
# 1. Load Required Libraries
# =============================================
if (!require("dplyr")) install.packages("dplyr")
if (!require("lubridate")) install.packages("lubridate")
if (!require("quadprog")) install.packages("quadprog")

library(dplyr)
library(lubridate)
library(quadprog)

# =============================================
# 2. Load and Preprocess Data
# =============================================

# --- Industry returns
industry_raw <- read.csv("./49_Industry_Portfolios.CSV")
industry <- industry_raw %>%
  mutate(Date = as.Date(paste0(Date, "01"), format = "%Y%m%d")) %>%
  mutate(across(-Date, ~ as.numeric(.) / 100)) # Make returns decimal

# --- Fama-French factors
factors <- read.csv("./F-F_Research_Data_Factors.CSV") %>%
  mutate(Date = as.Date(paste0(X, "01"), format = "%Y%m%d")) %>%
  rename(RMRF = Mkt.RF) %>%
  select(Date, RMRF, RF) %>%
  mutate(across(c(RMRF, RF), ~ as.numeric(.) / 100))

# --- Sydney macro uncertainty
macroU <- read.csv("./sydneymacrouncertainty.csv") %>%
  mutate(Date = parse_date_time(Date, orders = "m/Y")) %>%
  rename(MacroU = h.1) %>%
  select(Date, MacroU)

# --- 1-Year Treasury rate
treasury <- read.csv("./THREEFY1.csv") %>%
  rename(Date = observation_date, ShortRate = THREEFY1) %>%
  mutate(Date = as.Date(Date), ShortRate = as.numeric(ShortRate)) %>%
  select(Date, ShortRate)

# --- Trade Policy Uncertainty
tpu <- read.csv("./TPU INDEX MONTHLY.csv") %>%
  mutate(Date = as.Date(DATE, origin = "1899-12-30")) %>%
  select(Date, TPU)

# =============================================
# 3. Merge All Datasets
# =============================================

df_factors <- factors %>%
  inner_join(macroU, by = "Date") %>%
  inner_join(tpu, by = "Date") %>%
  inner_join(treasury, by = "Date")

df_all <- industry %>%
  inner_join(df_factors, by = "Date")

# =============================================
# 4. Compute Excess Returns
# =============================================

industry_cols <- setdiff(names(df_all), c("Date", "RMRF", "RF", "MacroU", "TPU", "ShortRate"))

df_all <- df_all %>%
  mutate(across(all_of(industry_cols), ~ . - RF))

# =============================================
# 5. Optimization Function: Mean-Variance using Market Beta as Proxy
# =============================================

optimize_meanvar_with_market_and_tpu <- function(Sigma, beta_market, beta_tpu, gamma = 50, target_market_beta = 1, target_tpu_beta = 0.1) {
  n <- length(beta_market)
  
  Dmat <- Sigma + diag(1e-6, n)  # stability
  dvec <- gamma * beta_market    # "reward" market beta
  
  Amat <- cbind(
    rep(1, n),           # sum weights = 1
    beta_market,         # market beta = 1
    beta_tpu             # TPU beta = 0.1
  )
  
  bvec <- c(1, target_market_beta, target_tpu_beta)
  
  meq <- 3  # all constraints are equality constraints
  
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq)
  w_opt <- sol$solution
  names(w_opt) <- names(beta_market)
  
  return(w_opt)
}

# =============================================
# 1. Load Daily Industry Return Data
# =============================================
daily_industry_raw <- read.csv("./49_Industry_Portfolios_Daily.csv")
names(daily_industry_raw)[1] <- "Date"   # <<< very important
daily_industry <- daily_industry_raw %>%
  mutate(Date = as.Date(as.character(Date), format = "%Y%m%d")) %>%
  mutate(across(-Date, ~ as.numeric(.) / 100))  # Make returns decimal

# =============================================
# 2. Backtest Setup
# =============================================

window_size <- 24  # months
start_idx <- which(df_all$Date >= df_all$Date[window_size])[1]
dates_backtest <- df_all$Date[(start_idx + 1):(nrow(df_all))]

returns_backtest <- c()

# New optimizer (with soft penalty)
optimize_soft_penalty <- function(Sigma, beta_market, beta_tpu,
                                  target_market_beta = 1,
                                  target_tpu_beta = 0.1,
                                  lambda = 5e-2) {
  n <- length(beta_market)
  
  Dmat <- Sigma + diag(lambda, n)
  dvec <- rep(0, n)
  
  Amat <- cbind(
    rep(1, n),
    beta_market,
    beta_tpu
  )
  bvec <- c(1, target_market_beta, target_tpu_beta)
  
  meq <- 3
  
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq)
  w_opt <- sol$solution
  names(w_opt) <- names(beta_market)
  return(w_opt)
}



# =============================================
# 2. Backtest Setup
# =============================================

window_size <- 24  # months
start_idx <- which(df_all$Date >= df_all$Date[window_size])[1]
dates_backtest <- df_all$Date[(start_idx + 1):(nrow(df_all))]

returns_backtest <- c()

# New optimizer (with soft penalty)
optimize_soft_penalty <- function(Sigma, beta_market, beta_tpu,
                                  target_market_beta = 1,
                                  target_tpu_beta = 0.1,
                                  lambda = 5e-2) {
  n <- length(beta_market)
  
  Dmat <- Sigma + diag(lambda, n)
  dvec <- rep(0, n)
  
  Amat <- cbind(
    rep(1, n),
    beta_market,
    beta_tpu
  )
  bvec <- c(1, target_market_beta, target_tpu_beta)
  
  meq <- 3
  
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq)
  w_opt <- sol$solution
  names(w_opt) <- names(beta_market)
  return(w_opt)
}

# =============================================
# 3. Backtesting Loop
# =============================================

for (t in seq_along(dates_backtest)) {
  date_now <- dates_backtest[t]
  
  idx_now <- which(df_all$Date == date_now)
  idx_start <- idx_now - window_size
  
  df_window_monthly <- df_all[idx_start:(idx_now-1), ]
  
  # Estimate monthly betas
  beta_market <- sapply(industry_cols, function(industry) {
    coef(lm(df_window_monthly[[industry]] ~ df_window_monthly$RMRF))[2]
  })
  
  beta_tpu <- sapply(industry_cols, function(industry) {
    coef(lm(df_window_monthly[[industry]] ~ df_window_monthly$TPU))[2]
  })
  
  # Estimate daily covariance
  window_start_date <- df_all$Date[idx_start]
  window_end_date <- df_all$Date[idx_now-1]
  
  df_window_daily <- daily_industry %>%
    filter(Date >= window_start_date & Date <= window_end_date)
  
  Sigma_daily <- cov(df_window_daily[, industry_cols]) * 252  # annualized
  
  # Optimization
  w_opt <- optimize_soft_penalty(Sigma_daily, beta_market, beta_tpu)
  
  # Next month return
  ret_next <- df_all[idx_now, industry_cols]
  returns_backtest[t] <- sum(w_opt * as.numeric(ret_next))
}

# =============================================
# 4. Performance Evaluation
# =============================================

cumulative_returns <- cumprod(1 + returns_backtest)
annualized_return <- mean(returns_backtest) * 12
annualized_vol <- sd(returns_backtest) * sqrt(12)
sharpe_ratio <- annualized_return / annualized_vol

cat("Annualized Return:", round(annualized_return, 4), "\n")
cat("Annualized Volatility:", round(annualized_vol, 4), "\n")
cat("Sharpe Ratio:", round(sharpe_ratio, 2), "\n")





# Plot cumulative return
plot(dates_backtest, cumulative_returns, type = "l", col = "blue",
     ylab = "Cumulative Return", xlab = "Date", main = "Cumulative Returns with Daily Sigma")

# Plot
plot(dates_backtest, cumulative_returns, type = "l", col = "blue",
     ylab = "Cumulative Return", xlab = "Date", main = "Cumulative Returns (Market Beta=1, TPU Beta=0.1)")


# =============================================
# Realized Beta Calculation Function
# =============================================

compute_realized_betas <- function(w_history, industry_cols, df_all, dates_backtest) {
  
  realized_market_beta <- c()
  realized_tpu_beta <- c()
  
  for (t in seq_along(dates_backtest)) {
    date_now <- dates_backtest[t]
    
    idx_now <- which(df_all$Date == date_now)
    idx_start <- idx_now - 24 + 1  # same rolling window size
    
    df_window <- df_all[idx_start:idx_now, ]
    
    # Regress portfolio return ~ Market (RMRF)
    port_ret <- rowSums(as.matrix(df_window[, industry_cols]) %*% w_history[[t]])
    
    model_market <- lm(port_ret ~ df_window$RMRF)
    realized_market_beta[t] <- coef(model_market)[2]
    
    # Regress portfolio return ~ TPU
    model_tpu <- lm(port_ret ~ df_window$TPU)
    realized_tpu_beta[t] <- coef(model_tpu)[2]
  }
  
  return(data.frame(Date = dates_backtest,
                    MarketBeta = realized_market_beta,
                    TPUBeta = realized_tpu_beta))
}

# =============================================
# Rebuild w_history
# =============================================

# Reconstruct w_history (weights at each t)
w_history <- list()

for (t in seq_along(dates_backtest)) {
  date_now <- dates_backtest[t]
  
  idx_now <- which(df_all$Date == date_now)
  idx_start <- idx_now - window_size
  
  df_window <- df_all[idx_start:(idx_now-1), ]
  
  Sigma_est <- cov(df_window[, industry_cols])
  
  beta_market <- sapply(industry_cols, function(industry) {
    coef(lm(df_window[[industry]] ~ df_window$RMRF))[2]
  })
  
  beta_tpu <- sapply(industry_cols, function(industry) {
    coef(lm(df_window[[industry]] ~ df_window$TPU))[2]
  })
  
  w_opt <- optimize_meanvar_with_market_and_tpu(Sigma_est, beta_market, beta_tpu,
                                                gamma = 50,
                                                target_market_beta = 1,
                                                target_tpu_beta = 0.1)
  
  w_history[[t]] <- w_opt
}

# =============================================
# Compute Realized Betas
# =============================================

realized_betas <- compute_realized_betas(w_history, industry_cols, df_all, dates_backtest)

# =============================================
# Plot Realized Market Beta
# =============================================

plot(realized_betas$Date, realized_betas$MarketBeta, type = "l", col = "blue",
     main = "Realized Portfolio Market Beta",
     xlab = "Date", ylab = "Market Beta")
abline(h = 1, col = "red", lty = 2)  # Target beta = 1

# =============================================
# Plot Realized TPU Beta
# =============================================

plot(realized_betas$Date, realized_betas$TPUBeta, type = "l", col = "purple",
     main = "Realized Portfolio TPU Beta",
     xlab = "Date", ylab = "TPU Beta")
abline(h = 0.1, col = "red", lty = 2)  # Target TPU beta = 0.1

sapply(1:10, function(k) {
  idx_now <- which(df_all$Date == dates_backtest[length(dates_backtest) - k])
  idx_start <- idx_now - window_size + 1
  var(df_all$RMRF[idx_start:idx_now])
})

summary(w_opt)
