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


full_data_monthly <- read.csv('Team3Data.csv')
full_data_monthly$Date <- as.Date(full_data_monthly$Date)

df_factors <- full_data_monthly %>%
  select(Date, RMRF, RF, MacroU, TPU, ShortRate, EPU, FinU)


# Industry columns
industry_cols <- setdiff(
  names(full_data_monthly),
  c("Date", "RMRF", "RF", "MacroU", "TPU", "ShortRate", "EPU", "FinU")
)


# =============================================
# 1. Define rolling window settings
# =============================================
lookback_months <- 24  # 2 years window

# Prepare storage for betas
rolling_betas <- list()
head(full_data_monthly)

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
    
    # Shrinkage parameter (how strong shrinkage is)
    lambda_shrink <- 0.2  #20% shrinkage
    
    # After getting the raw OLS coefficients:
    coefs_shrunk <- coefs  # initialize
    
    # Shrink each coefficient manually
    coefs_shrunk["RMRF"] <- (1 - lambda_shrink) * coefs["RMRF"] + lambda_shrink * 1
    coefs_shrunk["MacroU"] <- (1 - lambda_shrink) * coefs["MacroU"] + lambda_shrink * 0
    coefs_shrunk["FinU"] <- (1 - lambda_shrink) * coefs["FinU"] + lambda_shrink * 0
    coefs_shrunk["TPU"] <- (1 - lambda_shrink) * coefs["TPU"] + lambda_shrink * 0
    coefs_shrunk["EPU"] <- (1 - lambda_shrink) * coefs["EPU"] + lambda_shrink * 0
    coefs_shrunk["ShortRate"] <- (1 - lambda_shrink) * coefs["ShortRate"] + lambda_shrink * 0
    
    # Store shrunk betas
    betas_this_date$Beta_RMRF[j] <- coefs_shrunk["RMRF"]
    betas_this_date$Beta_MacroU[j] <- coefs_shrunk["MacroU"]
    betas_this_date$Beta_FinU[j] <- coefs_shrunk["FinU"]
    betas_this_date$Beta_TPU[j] <- coefs_shrunk["TPU"]
    betas_this_date$Beta_EPU[j] <- coefs_shrunk["EPU"]
    betas_this_date$Beta_ShortRate[j] <- coefs_shrunk["ShortRate"]
    
    
    # Store betas
    #betas_this_date$Beta_RMRF[j] <- coefs["RMRF"]
    #betas_this_date$Beta_MacroU[j] <- coefs["MacroU"]
    #betas_this_date$Beta_FinU[j] <- coefs["FinU"]
    #betas_this_date$Beta_TPU[j] <- coefs["TPU"]
    #betas_this_date$Beta_EPU[j] <- coefs["EPU"]
    #betas_this_date$Beta_ShortRate[j] <- coefs["ShortRate"]
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

optimize_portfolio_soft <- function(Sigma, betas, 
                                    target_market_beta = 0.5, 
                                    gamma_penalty = 1000) {
  
  n <- nrow(betas)
  
  # Add a small regularization to covariance matrix
  Sigma_reg <- Sigma + diag(1e-6, n)
  
  # Build penalty matrix for soft factor constraints
  P <- as.matrix(betas[, c("Beta_MacroU", "Beta_FinU", "Beta_TPU", "Beta_EPU", "Beta_ShortRate")])
  PenaltyMat <- gamma_penalty * (P %*% t(P))  # <-- CORRECTED
  
  Dmat <- Sigma_reg + PenaltyMat
  dvec <- rep(0, n)
  
  # Only hard constraints: full investment and market beta = 1
  Amat <- t(rbind(
    rep(1, n),
    betas$Beta_RMRF
  ))
  
  bvec <- c(1, target_market_beta)
  meq <- length(bvec)  # equality constraints
  
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq)
  
  w_opt <- sol$solution
  w_opt <- w_opt / sum(w_opt)  # Normalize weights again
  
  return(w_opt)
}


optimize_and_backtest <- function(target_market_beta_value) {
  
  # Initialize storage
  portfolio_returns <- rep(NA, nrow(full_data_monthly))
  portfolio_weights <- matrix(NA, nrow = nrow(full_data_monthly), ncol = n_assets)
  colnames(portfolio_weights) <- industry_cols
  
  # Loop
  for (t in (lookback_months + 1):nrow(full_data_monthly)) {
    
    if ((t - lookback_months - 1) %% rebalance_every != 0 && t > (lookback_months + 1)) {
      portfolio_weights[t, ] <- portfolio_weights[t - 1, ]
      portfolio_returns[t] <- sum(portfolio_weights[t, ] * full_data_monthly[t, industry_cols], na.rm = TRUE)
      next
    }
    
    betas_this_date <- beta_estimates %>%
      filter(Date == full_data_monthly$Date[t]) %>%
      arrange(match(Industry, industry_cols))
    
    if (nrow(betas_this_date) != n_assets) {
      portfolio_weights[t, ] <- rep(1/n_assets, n_assets)
      portfolio_returns[t] <- mean(as.numeric(full_data_monthly[t, industry_cols]), na.rm = TRUE)
      next
    }
    
    window_data <- full_data_monthly[(t - lookback_months):(t - 1), industry_cols]
    Sigma_est <- cov(window_data) * 12
    Sigma_est <- Sigma_est + diag(1e-6, ncol(Sigma_est))
    
    w_opt <- optimize_portfolio_soft(Sigma_est, betas_this_date,
                                     target_market_beta = target_market_beta_value,
                                     gamma_penalty = 10)
    
    portfolio_weights[t, ] <- w_opt
    portfolio_returns[t] <- sum(w_opt * full_data_monthly[t, industry_cols], na.rm = TRUE)
  }
  
  # Build final results
  results_tmp <- data.frame(
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
  
  # return BOTH results and weights
  return(list(results = results_tmp, weights = portfolio_weights))
}


# Run optimization for each target beta
n_assets <- length(industry_cols) #number of indusries
rebalance_every <- 6 #rebalance every 6 months

# Run optimization
opt_05 <- optimize_and_backtest(target_market_beta_value = 0.5)
opt_10 <- optimize_and_backtest(target_market_beta_value = 1.0)
opt_15 <- optimize_and_backtest(target_market_beta_value = 1.5)

# Now extract results and weights separately:
results_05 <- opt_05$results
portfolio_weights_05 <- opt_05$weights

results_10 <- opt_10$results
portfolio_weights_10 <- opt_10$weights

results_15 <- opt_15$results
portfolio_weights_15 <- opt_15$weights

names(full_data_monthly)
names(results_10)
names(results_15)
names(results_15)

# =============================================
# 4. Quick plot
# =============================================
library(ggplot2)
library(dplyr)
library(tidyr)

# Stack them together
results_combined <- bind_rows(
  results_05 %>% mutate(Portfolio = "Optimized (Beta = 0.5)"),
  results_10 %>% mutate(Portfolio = "Optimized (Beta = 1.0)"),
  results_15 %>% mutate(Portfolio = "Optimized (Beta = 1.5)")
)

# Reshape
results_long <- results_combined %>%
  select(Date, Portfolio, Portfolio_Cumulative, Market_Cumulative) %>%
  pivot_longer(cols = c(Portfolio_Cumulative, Market_Cumulative),
               names_to = "Series", values_to = "CumulativeReturn") %>%
  mutate(LineType = ifelse(Series == "Market_Cumulative", "Market", Portfolio))

library(ggplot2)

ggplot(results_long, aes(x = Date, y = CumulativeReturn, color = LineType)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(
    "Optimized (Beta = 0.5)" = "blue",
    "Optimized (Beta = 1.0)" = "green",
    "Optimized (Beta = 1.5)" = "purple",
    "Market" = "red"
  )) +
  labs(
    title = "Dynamic Portfolios (Different Market Betas) vs Market Cumulative Return",
    y = "Growth of $1",
    x = "Date",
    color = "Portfolio"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

compute_stats <- function(df, name) {
  ann_return <- round(mean(df$Portfolio_Return) * 12,2)
  ann_volatility <- round(sd(df$Portfolio_Return) * sqrt(12),2)
  sharpe_ratio <- round(mean(df$Portfolio_Excess)*12 /(sd(df$Portfolio_Excess)*sqrt(12)),2)
  
  return(data.frame(
    Portfolio = name,
    Annualized_Return = ann_return,
    Annualized_Volatility = ann_volatility,
    Sharpe_Ratio = sharpe_ratio
  ))
}

# For each optimized portfolio
stats_05 <- compute_stats(results_05, "Optimized (Beta = 0.5)")
stats_10 <- compute_stats(results_10, "Optimized (Beta = 1.0)")
stats_15 <- compute_stats(results_15, "Optimized (Beta = 1.5)")

# For Market (use Market_Return!)
compute_market_stats <- function(df) {
  ann_return <- round(mean(df$Market_Return) * 12,2)
  ann_volatility <- round(sd(df$Market_Return) * sqrt(12),2)
  sharpe_ratio <- round(mean(df$Market_Excess)*12 / (sd(df$Market_Excess)*sqrt(12)),2)
  
  return(data.frame(
    Portfolio = "Market",
    Annualized_Return = ann_return,
    Annualized_Volatility = ann_volatility,
    Sharpe_Ratio = sharpe_ratio
  ))
}

stats_market <- compute_market_stats(results_05)  # Market return is same across all

# Combine all
stats_all <- bind_rows(stats_05, stats_10, stats_15, stats_market)

# Print table
print(stats_all)

# Function to calculate portfolio betas dynamically
compute_portfolio_betas <- function(portfolio_weights, beta_estimates) {
  betas_list <- list()
  
  for (t in 1:nrow(portfolio_weights)) {
    w <- portfolio_weights[t, ]
    
    # Find matching betas
    beta_t <- beta_estimates %>%
      filter(Date == full_data_monthly$Date[t]) %>%
      arrange(match(Industry, industry_cols))
    
    if (nrow(beta_t) == length(industry_cols)) {
      beta_portfolio <- colSums(w * as.matrix(beta_t[, c(
        "Beta_RMRF", "Beta_MacroU", "Beta_FinU", "Beta_TPU", "Beta_EPU", "Beta_ShortRate"
      )]))
      
      betas_list[[length(betas_list) + 1]] <- data.frame(
        Date = full_data_monthly$Date[t],
        Beta_Market = beta_portfolio["Beta_RMRF"],
        Beta_MacroU = beta_portfolio["Beta_MacroU"],
        Beta_FinU = beta_portfolio["Beta_FinU"],
        Beta_TPU = beta_portfolio["Beta_TPU"],
        Beta_EPU = beta_portfolio["Beta_EPU"],
        Beta_ShortRate = beta_portfolio["Beta_ShortRate"]
      )
    }
  }
  
  betas_df <- do.call(rbind, betas_list)
  return(betas_df)
}


# Now apply to each set of portfolio weights
betas_05 <- compute_portfolio_betas(portfolio_weights_05, beta_estimates)
betas_10 <- compute_portfolio_betas(portfolio_weights_10, beta_estimates)
betas_15 <- compute_portfolio_betas(portfolio_weights_15, beta_estimates)

# From previous steps
betas_all <- bind_rows(
  betas_05 %>% mutate(Portfolio = "Beta = 0.5"),
  betas_10 %>% mutate(Portfolio = "Beta = 1.0"),
  betas_15 %>% mutate(Portfolio = "Beta = 1.5")
)

# Check if it worked
# 1. Plot for Macro Uncertainty Beta
p_macroU <- ggplot(betas_all, aes(x = Date, y = Beta_MacroU, color = Portfolio)) +
  geom_line() +
  labs(title = "Optimized Beta to Macro Uncertainty Over Time",
       x = "Date",
       y = "Beta to MacroU") +
  theme_minimal()

# 2. Plot for Financial Uncertainty Beta
p_finU <- ggplot(betas_all, aes(x = Date, y = Beta_FinU, color = Portfolio)) +
  geom_line() +
  labs(title = "Optimized Beta to Financial Uncertainty Over Time",
       x = "Date",
       y = "Beta to FinU") +
  theme_minimal()

# 3. Plot for Trade Policy Uncertainty Beta
p_tpu <- ggplot(betas_all, aes(x = Date, y = Beta_TPU, color = Portfolio)) +
  geom_line() +
  labs(title = "Optimized Beta to Trade Policy Uncertainty Over Time",
       x = "Date",
       y = "Beta to TPU") +
  theme_minimal()

# Display all together
library(patchwork)
(p_macroU / p_finU / p_tpu) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# =============================================
# Plotting Realized Betas for Each Industry Over Time
# =============================================
# Convert to long format for ggplot2
betas_long <- beta_estimates %>%
  pivot_longer(
    cols = starts_with("Beta_"),
    names_to = "Factor",
    values_to = "Beta"
  )

# Plot 1: Macro Uncertainty Beta Distribution Over Time
p_macroU <- betas_long %>%
  filter(Factor == "Beta_MacroU") %>%
  ggplot(aes(x = Date, y = Beta, group = Date)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "Macro Uncertainty Betas Across Industries Over Time",
       y = "Beta to MacroU", x = "Date") +
  theme_minimal()

# Plot 2: Trade Policy Uncertainty Beta Distribution Over Time
p_tpu <- betas_long %>%
  filter(Factor == "Beta_TPU") %>%
  ggplot(aes(x = Date, y = Beta, group = Date)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "Trade Policy Uncertainty Betas Across Industries Over Time",
       y = "Beta to TPU", x = "Date") +
  theme_minimal()

# Plot 3: Economic Policy Uncertainty Beta Distribution Over Time
p_epu <- betas_long %>%
  filter(Factor == "Beta_EPU") %>%
  ggplot(aes(x = Date, y = Beta, group = Date)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "Economic Policy Uncertainty Betas Across Industries Over Time",
       y = "Beta to EPU", x = "Date") +
  theme_minimal()

# Plot 4: Financial Uncertainty Beta Distribution Over Time
p_finu <- betas_long %>%
  filter(Factor == "Beta_FinU") %>%
  ggplot(aes(x = Date, y = Beta, group = Date)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "Financial Uncertainty Betas Across Industries Over Time",
       y = "Beta to FinU", x = "Date") +
  theme_minimal()

# Plot 5: ShortRate Beta Distribution Over Time
p_sr <- betas_long %>%
  filter(Factor == "Beta_ShortRate") %>%
  ggplot(aes(x = Date, y = Beta, group = Date)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "ShortRate Betas Across Industries Over Time",
       y = "Beta to ShortRate", x = "Date") +
  theme_minimal()

# Print the plots
p_macroU
p_tpu
p_epu
p_sr
p_finu

# Realized Beta comparison 
# =============================================
# 1. Prepare the returns data
# =============================================

# Add the portfolios together
all_returns <- bind_rows(
  results_05 %>% mutate(Portfolio = "Beta_0.5"),
  results_10 %>% mutate(Portfolio = "Beta_1.0"),
  results_15 %>% mutate(Portfolio = "Beta_1.5")
)

# Add Market as if it is a portfolio too
market_returns <- results_05 %>%
  mutate(
    Portfolio = "Market",
    Portfolio_Return = Market_Return
  )

# Combine
full_returns <- bind_rows(all_returns, market_returns)

# =============================================
# 2. Define Rolling Regression Function
# =============================================

run_rolling_betas <- function(portfolio_name, lookback_months = 24) {
  
  df_returns <- full_returns %>% filter(Portfolio == portfolio_name) %>%
    select(Date, Portfolio_Return)
  
  df <- df_returns %>%
    inner_join(df_factors, by = "Date") %>%
    select(Date, Portfolio_Return, RF, MacroU, FinU, TPU, EPU, ShortRate)
  
  betas_list <- list()
  
  for (i in (lookback_months + 1):nrow(df)) {
    
    window_data <- df[(i - lookback_months):(i - 1), ]
    
    y <- window_data$Portfolio_Return - window_data$RF  # Excess return
    x <- window_data %>%
      select(MacroU, FinU, TPU, EPU, ShortRate)
    
    model <- lm(y ~ ., data = x)
    
    coefs <- coef(model)
    
    betas_list[[length(betas_list) + 1]] <- data.frame(
      Date = df$Date[i],
      Beta_MacroU = coefs["MacroU"],
      Beta_FinU = coefs["FinU"],
      Beta_TPU = coefs["TPU"],
      Beta_EPU = coefs["EPU"],
      Beta_ShortRate = coefs["ShortRate"],
      Portfolio = portfolio_name
    )
  }
  
  betas_df <- do.call(rbind, betas_list)
  return(betas_df)
}

# =============================================
# 3. Run for each portfolio and market
# =============================================

betas_05_realized <- run_rolling_betas("Beta_0.5")
betas_10_realized <- run_rolling_betas("Beta_1.0")
betas_15_realized <- run_rolling_betas("Beta_1.5")
betas_market_realized <- run_rolling_betas("Market")

# Combine all
betas_realized_all <- bind_rows(
  betas_05_realized,
  betas_10_realized,
  betas_15_realized,
  betas_market_realized
)

# View
head(betas_realized_all)

# =============================================
# Corrected Boxplot + Line Plot
# =============================================
lookback_months <- 24  # Rolling 24-month window
rolling_betas_list <- list() # Initiate list
# ========================
# 2. Loop over industries
# ========================
for (ind in industry_cols) {
  
  betas_list <- list()
  
  for (i in (lookback_months + 1):nrow(full_data_monthly)) {
    
    window_data <- full_data_monthly[(i - lookback_months):(i - 1), ]
    
    y <- window_data[[ind]] - window_data$RF  # <-- excess returns if needed
    x <- window_data %>%
      select(MacroU, FinU, TPU, EPU, ShortRate)
    
    if (nrow(window_data) == lookback_months && !all(is.na(y))) {  # Only if full window and y not all NA
      model <- lm(y ~ ., data = x)
      coefs <- coef(model)
      
      betas_list[[length(betas_list) + 1]] <- data.frame(
        Date = full_data_monthly$Date[i],
        Industry = ind,
        Beta_MacroU = coefs["MacroU"],
        Beta_FinU = coefs["FinU"],
        Beta_TPU = coefs["TPU"],
        Beta_EPU = coefs["EPU"],
        Beta_ShortRate = coefs["ShortRate"]
      )
    }
  }
  
  rolling_betas_list[[ind]] <- do.call(rbind, betas_list)
}
# ========================
# 3. After loop: combine everything
# ========================
industry_betas_all <- do.call(rbind, rolling_betas_list)
industry_betas_long <- industry_betas_all %>%
  pivot_longer(
    cols = starts_with("Beta_"),
    names_to = "Factor",
    values_to = "Beta"
  )

# Combine optimized portfolios + market
portfolio_betas_all <- bind_rows(
  betas_all %>% mutate(Portfolio = case_when(
    Portfolio == "Beta = 0.5" ~ "Beta_0.5",
    Portfolio == "Beta = 1.0" ~ "Beta_1.0",
    Portfolio == "Beta = 1.5" ~ "Beta_1.5",
    TRUE ~ Portfolio
  )),
  betas_market_realized %>% mutate(Portfolio = "Market")
)

portfolio_betas_long <- portfolio_betas_all %>%
  pivot_longer(
    cols = starts_with("Beta_"),
    names_to = "Factor",
    values_to = "Beta"
  )

unique(portfolio_betas_long$Portfolio)

plot_factor_beta_lines <- function(factor_name, y_label) {
  
  # Filter industry and portfolio data
  industry_data <- industry_betas_long %>% filter(Factor == factor_name)
  portfolio_data <- portfolio_betas_long %>% filter(Factor == factor_name)
  
  # Combine beta values
  all_beta_values <- c(industry_data$Beta, portfolio_data$Beta)
  
  # Set limits based on 5th and 95th percentiles
  beta_min <- quantile(all_beta_values, 0.05, na.rm = TRUE)
  beta_max <- quantile(all_beta_values, 0.95, na.rm = TRUE)
  
  ggplot() +
    geom_boxplot(
      data = industry_data,
      aes(x = Date, y = Beta, group = Date),
      outlier.shape = NA,
      fill = "lightgray",
      color = "black",
      width = 0.5
    ) +
    geom_line(
      data = portfolio_data,
      aes(x = Date, y = Beta, color = Portfolio),
      size = 1
    ) +
    labs(title = paste("Industry vs Optimized Portfolios' Betas to", y_label),
         x = "Date", y = paste("Beta to", y_label)) +
    coord_cartesian(ylim = c(beta_min, beta_max)) +   # <<<<<< USE coord_cartesian better than ylim
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 14)) +
    scale_color_manual(values = c(
      "Beta_0.5" = "red",
      "Beta_1.0" = "green",
      "Beta_1.5" = "blue",
      "Market" = "purple"
    ))
}



# Now plot each factor
p_macroU_line <- plot_factor_beta_lines("Beta_MacroU", "Macro Uncertainty")
p_finU_line   <- plot_factor_beta_lines("Beta_FinU", "Financial Uncertainty")
p_tpu_line    <- plot_factor_beta_lines("Beta_TPU", "Trade Policy Uncertainty")
p_epu_line    <- plot_factor_beta_lines("Beta_EPU", "Economic Policy Uncertainty")
p_sr_line     <- plot_factor_beta_lines("Beta_ShortRate", "Short Rate")


# =============================================
# Display the plots
# =============================================
p_macroU_line
p_finU_line
p_tpu_line
p_epu_line
p_sr_line

library(dplyr)

# 1. Filter for the portfolios we care about
betas_comparison_all <- betas_realized_all %>%
  filter(Portfolio %in% c("Beta_0.5", "Beta_1.0", "Beta_1.5", "Market"))

# 2. Summarize the average beta for each factor
betas_summary_all <- betas_comparison_all %>%
  group_by(Portfolio) %>%
  summarise(
    Avg_Beta_MacroU = mean(Beta_MacroU, na.rm = TRUE),
    Avg_Beta_FinU = mean(Beta_FinU, na.rm = TRUE),
    Avg_Beta_TPU = mean(Beta_TPU, na.rm = TRUE),
    Avg_Beta_EPU = mean(Beta_EPU, na.rm = TRUE),
    Avg_Beta_ShortRate = mean(Beta_ShortRate, na.rm = TRUE)
  )

# 3. Extract Market average betas
market_betas <- betas_summary_all %>%
  filter(Portfolio == "Market") %>%
  select(-Portfolio) %>%
  as.list()

# 4. Compute Improvement Percentages
betas_summary_all <- betas_summary_all %>%
  mutate(
    MacroU_Improvement_Percent = (market_betas$Avg_Beta_MacroU - Avg_Beta_MacroU) / abs(market_betas$Avg_Beta_MacroU) * 100,
    FinU_Improvement_Percent = (market_betas$Avg_Beta_FinU - Avg_Beta_FinU) / abs(market_betas$Avg_Beta_FinU) * 100,
    TPU_Improvement_Percent = (market_betas$Avg_Beta_TPU - Avg_Beta_TPU) / abs(market_betas$Avg_Beta_TPU) * 100,
    EPU_Improvement_Percent = (market_betas$Avg_Beta_EPU - Avg_Beta_EPU) / abs(market_betas$Avg_Beta_EPU) * 100,
    ShortRate_Improvement_Percent = (market_betas$Avg_Beta_ShortRate - Avg_Beta_ShortRate) / abs(market_betas$Avg_Beta_ShortRate) * 100
  )

# 5. View result
print(betas_summary_all)

# Extract the market's betas
market_betas <- betas_summary_all %>%
  filter(Portfolio == "Market") %>%
  select(starts_with("Avg_Beta")) %>%
  as.list()

# Recompute the improvements correctly (use absolute values for comparison)
betas_summary_all <- betas_summary_all %>%
  mutate(
    MacroU_Improvement_Percent = (abs(market_betas$Avg_Beta_MacroU) - abs(Avg_Beta_MacroU)) / abs(market_betas$Avg_Beta_MacroU) * 100,
    FinU_Improvement_Percent = (abs(market_betas$Avg_Beta_FinU) - abs(Avg_Beta_FinU)) / abs(market_betas$Avg_Beta_FinU) * 100,
    TPU_Improvement_Percent = (abs(market_betas$Avg_Beta_TPU) - abs(Avg_Beta_TPU)) / abs(market_betas$Avg_Beta_TPU) * 100,
    EPU_Improvement_Percent = (abs(market_betas$Avg_Beta_EPU) - abs(Avg_Beta_EPU)) / abs(market_betas$Avg_Beta_EPU) * 100,
    ShortRate_Improvement_Percent = (abs(market_betas$Avg_Beta_ShortRate) - abs(Avg_Beta_ShortRate)) / abs(market_betas$Avg_Beta_ShortRate) * 100
  )

# View nicely only improvements
betas_improvement_table <- betas_summary_all %>%
  select(Portfolio,
         MacroU_Improvement_Percent,
         FinU_Improvement_Percent,
         TPU_Improvement_Percent,
         EPU_Improvement_Percent,
         ShortRate_Improvement_Percent)

print(betas_improvement_table)



# Filter only MacroU and FinU improvement
betas_improvement_focus <- betas_summary_all %>%
  select(Portfolio,
         MacroU_Improvement_Percent,
         FinU_Improvement_Percent)

print(betas_improvement_focus)

# Step 1: Reshape it to long format including MacroU, FinU, TPU, and EPU
betas_improvement_focus_long <- betas_summary_all %>%
  select(Portfolio, MacroU_Improvement_Percent, FinU_Improvement_Percent,
         TPU_Improvement_Percent, EPU_Improvement_Percent) %>%
  pivot_longer(
    cols = c(MacroU_Improvement_Percent, FinU_Improvement_Percent,
             TPU_Improvement_Percent, EPU_Improvement_Percent),
    names_to = "Factor",
    values_to = "Improvement_Percent"
  )

# Remove Market from the dataset
betas_improvement_focus_long_clean <- betas_improvement_focus_long %>%
  filter(Portfolio != "Market")

# Plot without Market
ggplot(betas_improvement_focus_long_clean, aes(x = Factor, y = Improvement_Percent, fill = Portfolio)) +
  geom_col(position = "dodge") +
  labs(title = "Improvement (%) of Each Optimized Portfolio Against Market",
       x = "Factor",
       y = "Improvement (%)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14)
  )

# ================================
# 1. Calculate Sensitivity Score
# ================================

# If you are using beta_estimates for industries:
industry_sensitivity <- beta_estimates %>%
  group_by(Industry) %>%
  summarize(
    Avg_Beta_MacroU = mean(Beta_MacroU, na.rm = TRUE),
    Avg_Beta_FinU = mean(Beta_FinU, na.rm = TRUE)
  ) %>%
  mutate(
    Uncertainty_Sensitivity = (abs(Avg_Beta_MacroU) + abs(Avg_Beta_FinU)) / 2
  )

# If you are using portfolio realized betas:
portfolio_sensitivity <- betas_realized_all %>%
  group_by(Portfolio) %>%
  summarize(
    Avg_Beta_MacroU = mean(Beta_MacroU, na.rm = TRUE),
    Avg_Beta_FinU = mean(Beta_FinU, na.rm = TRUE)
  ) %>%
  mutate(
    Uncertainty_Sensitivity = (abs(Avg_Beta_MacroU) + abs(Avg_Beta_FinU)) / 2
  )

# ================================
# 2. Combine and Prepare for Plot
# ================================

# Add a label to know which is portfolio or industry
industry_sensitivity <- industry_sensitivity %>%
  mutate(Type = "Industry", Name = Industry) %>%
  select(Name, Uncertainty_Sensitivity, Type)

portfolio_sensitivity <- portfolio_sensitivity %>%
  mutate(Type = "Portfolio", Name = Portfolio) %>%
  select(Name, Uncertainty_Sensitivity, Type)

# Combine
sensitivity_all <- bind_rows(industry_sensitivity, portfolio_sensitivity)

# ================================
# 3. Plot Ranking Bar Chart
# ================================

# Sort descending for plotting
sensitivity_all <- sensitivity_all %>%
  arrange(desc(Uncertainty_Sensitivity)) %>%
  mutate(Name = factor(Name, levels = rev(Name)))  # <<< reverse order!

# Plot
#ggplot(sensitivity_all, aes(x = Name, y = Uncertainty_Sensitivity, fill = Type)) +
#  geom_bar(stat = "identity") +
#  coord_flip() +
#  labs(title = "Uncertainty Sensitivity Score (Equal Weighted MacroU and FinU)",
#       x = "Portfolio / Industry",
#       y = "Sensitivity Score") +
#  theme_minimal() +
#  theme(plot.title = element_text(hjust = 0.5)) +
#  scale_fill_manual(values = c("Portfolio" = "steelblue", "Industry" = "lightgray"))

# 1. Average sensitivity for industries
industry_avg_sensitivity <- industry_betas_all %>%
  group_by(Industry) %>%
  summarise(Uncertainty_Sensitivity = mean((abs(Beta_MacroU) + abs(Beta_FinU) + abs(Beta_TPU) + abs(Beta_EPU)) / 4, na.rm = TRUE)) %>%
  mutate(Type = "Industry") %>%
  rename(Name = Industry)

# 2. Average sensitivity for portfolios
portfolio_avg_sensitivity <- betas_realized_all %>%
  filter(Portfolio %in% c("Beta_0.5", "Beta_1.0", "Beta_1.5", "Market")) %>%
  group_by(Portfolio) %>%
  summarise(Uncertainty_Sensitivity = mean((abs(Beta_MacroU) + abs(Beta_FinU) + abs(Beta_TPU) + abs(Beta_EPU)) / 4, na.rm = TRUE)) %>%
  mutate(Type = "Portfolio") %>%
  rename(Name = Portfolio)

# 3. Combine all
sensitivity_all <- bind_rows(industry_avg_sensitivity, portfolio_avg_sensitivity)

# 4. Sort descending
sensitivity_all <- sensitivity_all %>%
  arrange(desc(Uncertainty_Sensitivity)) %>%
  mutate(Name = factor(Name, levels = rev(Name)))

# 5. Plot
#ggplot(sensitivity_all, aes(x = Name, y = Uncertainty_Sensitivity, fill = Type)) +
#  geom_bar(stat = "identity") +
#  coord_flip() +
#  labs(title = "Uncertainty Sensitivity Score (Equal Weighted MacroU, FinU, TPU, EPU)",
#       x = "Portfolio / Industry",
#       y = "Sensitivity Score") +
#  theme_minimal() +
#  theme(plot.title = element_text(hjust = 0.5)) +
#  scale_fill_manual(values = c("Portfolio" = "steelblue", "Industry" = "lightgray"))

# 1. Average sensitivity for industries (EPU and TPU only)
industry_avg_sensitivity_epu_tpu <- industry_betas_all %>%
  group_by(Industry) %>%
  summarise(Uncertainty_Sensitivity = mean((abs(Beta_TPU) + abs(Beta_EPU)) / 2, na.rm = TRUE)) %>%
  mutate(Type = "Industry") %>%
  rename(Name = Industry)

# 2. Average sensitivity for portfolios (EPU and TPU only)
portfolio_avg_sensitivity_epu_tpu <- betas_realized_all %>%
  filter(Portfolio %in% c("Beta_0.5", "Beta_1.0", "Beta_1.5", "Market")) %>%
  group_by(Portfolio) %>%
  summarise(Uncertainty_Sensitivity = mean((abs(Beta_TPU) + abs(Beta_EPU)) / 2, na.rm = TRUE)) %>%
  mutate(Type = "Portfolio") %>%
  rename(Name = Portfolio)

# 3. Combine all
sensitivity_all_epu_tpu <- bind_rows(industry_avg_sensitivity_epu_tpu, portfolio_avg_sensitivity_epu_tpu)

# 4. Sort descending
sensitivity_all_epu_tpu <- sensitivity_all_epu_tpu %>%
  arrange(desc(Uncertainty_Sensitivity)) %>%
  mutate(Name = factor(Name, levels = rev(Name)))

# 5. Plot
#ggplot(sensitivity_all_epu_tpu, aes(x = Name, y = Uncertainty_Sensitivity, fill = Type)) +
#  geom_bar(stat = "identity") +
#  coord_flip() +
#  labs(title = "Uncertainty Sensitivity Score (Equal Weighted TPU and EPU)",
#       x = "Portfolio / Industry",
#       y = "Sensitivity Score") +
#  theme_minimal() +
#  theme(plot.title = element_text(hjust = 0.5)) +
#  scale_fill_manual(values = c("Portfolio" = "steelblue", "Industry" = "lightgray"))


# First, rename inside each dataset if needed
industry_betas_all_fixed <- industry_betas_all %>%
  rename(Portfolio = Industry)

betas_05_realized_fixed <- betas_05_realized %>%
  mutate(Portfolio = "Beta_0.5")

betas_10_realized_fixed <- betas_10_realized %>%
  mutate(Portfolio = "Beta_1.0")

betas_15_realized_fixed <- betas_15_realized %>%
  mutate(Portfolio = "Beta_1.5")

betas_market_realized_fixed <- betas_market_realized %>%
  mutate(Portfolio = "Market")

# Now safely bind everything together
betas_realized_all <- bind_rows(
  industry_betas_all_fixed,
  betas_05_realized_fixed,
  betas_10_realized_fixed,
  betas_15_realized_fixed,
  betas_market_realized_fixed
)


# Step 3: Standardize (scale) the factors
factors_scaled <- full_data_monthly %>%
  select(MacroU, FinU, TPU, EPU) %>%
  scale() %>% as.data.frame()

factor_sd <- apply(factors_scaled, 2, sd)  # Should be 1 already if scaled

# Step 4: Adjust betas
adjusted_betas_all <- betas_realized_all %>%
  mutate(
    Beta_MacroU = Beta_MacroU / sd(factors_scaled$MacroU),
    Beta_FinU   = Beta_FinU / sd(factors_scaled$FinU),
    Beta_TPU    = Beta_TPU / sd(factors_scaled$TPU),
    Beta_EPU    = Beta_EPU / sd(factors_scaled$EPU)
  )

# Step 5: Compute average absolute sensitivity per portfolio
sensitivity_summary <- adjusted_betas_all %>%
  group_by(Portfolio) %>%
  summarize(
    Sensitivity_Score = mean(abs(Beta_MacroU) + abs(Beta_FinU) + abs(Beta_TPU) + abs(Beta_EPU)) / 4,
    .groups = "drop"
  )

# Step 6: Add Type
sensitivity_summary <- sensitivity_summary %>%
  mutate(Type = ifelse(Portfolio %in% c("Beta_0.5", "Beta_1.0", "Beta_1.5", "Market"), "Portfolio", "Industry"))

# Step 7: Plot

ggplot(sensitivity_summary %>% arrange(Sensitivity_Score), 
       aes(x = reorder(Portfolio, Sensitivity_Score), y = Sensitivity_Score, fill = Type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Adjusted Uncertainty Sensitivity Score (Standardized Factors)",
       x = "Portfolio / Industry",
       y = "Sensitivity Score") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Portfolio" = "steelblue", "Industry" = "lightgray"))

# Step 5: Compute average absolute sensitivity for only TPU and EPU
sensitivity_summary_tpu_epu <- adjusted_betas_all %>%
  group_by(Portfolio) %>%
  summarize(
    Sensitivity_Score = mean(abs(Beta_TPU) + abs(Beta_EPU)) / 2,   # <--- ONLY TPU and EPU
    .groups = "drop"
  )

# Step 6: Add Type (Portfolio vs Industry)
sensitivity_summary_tpu_epu <- sensitivity_summary_tpu_epu %>%
  mutate(Type = ifelse(Portfolio %in% c("Beta_0.5", "Beta_1.0", "Beta_1.5", "Market"), "Portfolio", "Industry"))

# Step 7: Plot
ggplot(sensitivity_summary_tpu_epu %>% arrange(Sensitivity_Score), 
       aes(x = reorder(Portfolio, Sensitivity_Score), y = Sensitivity_Score, fill = Type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Adjusted Trade & Economic Policy Uncertainty Sensitivity (Standardized Factors)",
       x = "Portfolio / Industry",
       y = "Sensitivity Score") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Portfolio" = "steelblue", "Industry" = "lightgray"))



# ===========================
# 1. Compute realized beta over time
# ===========================

# For Beta = 0.5 portfolio
realized_beta_05 <- betas_05 %>%
  mutate(Target_Beta = 0.5) %>%
  select(Date, Realized_Beta = Beta_Market, Target_Beta) %>%
  mutate(Portfolio = "Beta_0.5")

# For Beta = 1.0 portfolio
realized_beta_10 <- betas_10 %>%
  mutate(Target_Beta = 1.0) %>%
  select(Date, Realized_Beta = Beta_Market, Target_Beta) %>%
  mutate(Portfolio = "Beta_1.0")

# For Beta = 1.5 portfolio
realized_beta_15 <- betas_15 %>%
  mutate(Target_Beta = 1.5) %>%
  select(Date, Realized_Beta = Beta_Market, Target_Beta) %>%
  mutate(Portfolio = "Beta_1.5")

# ===========================
# 2. Combine all
# ===========================

realized_betas_all <- bind_rows(realized_beta_05, realized_beta_10, realized_beta_15)

# ===========================
# 3. Calculate Average Realized Beta for each Portfolio
# ===========================

average_betas <- realized_betas_all %>%
  group_by(Portfolio) %>%
  summarize(Average_Realized_Beta = mean(Realized_Beta, na.rm = TRUE), .groups = "drop")

# Merge the average back for plotting
realized_betas_all <- realized_betas_all %>%
  left_join(average_betas, by = "Portfolio")

# ===========================
# 4. Plot realized vs target beta + average line
# ===========================

p <- ggplot(realized_betas_all, aes(x = Date, y = Realized_Beta, color = Portfolio)) +
  geom_line(size = 1) +  # <- This line plots the colored lines
  geom_hline(aes(yintercept = Target_Beta, linetype = "Target Beta"), color = "black", size = 1) +
  geom_hline(aes(yintercept = Average_Realized_Beta, linetype = "Average Realized Beta"), color = "gray30", size = 1) +
  facet_wrap(~ Portfolio, ncol = 1, scales = "free_y") +
  labs(
    title = "Realized Beta vs Target Beta",
    x = "Date",
    y = "Beta to Market (RMRF)",
    linetype = "Reference Line",
    color = "Portfolio"
  ) +
  scale_linetype_manual(values = c("Target Beta" = "dashed", "Average Realized Beta" = "dotted")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "bottom"
  )

# Print
print(p)

print(betas_summary_all)

# Assuming you already have factors_scaled from your earlier scaling
factor_sd <- apply(factors_scaled, 2, sd)

# Step 1: Normalize betas
betas_normalized <- betas_summary_all %>%
  mutate(
    Avg_Beta_MacroU = Avg_Beta_MacroU / sd(factors_scaled$MacroU),
    Avg_Beta_FinU = Avg_Beta_FinU / sd(factors_scaled$FinU),
    Avg_Beta_TPU = Avg_Beta_TPU / sd(factors_scaled$TPU),
    Avg_Beta_EPU = Avg_Beta_EPU / sd(factors_scaled$EPU)
  )

# Step 2: Reshape
betas_normalized_long <- betas_normalized %>%
  select(Portfolio, Avg_Beta_MacroU, Avg_Beta_FinU, Avg_Beta_TPU, Avg_Beta_EPU) %>%
  pivot_longer(
    cols = starts_with("Avg_Beta"),
    names_to = "Factor",
    values_to = "Normalized_Beta"
  ) %>%
  mutate(Factor = recode(Factor,
                         "Avg_Beta_MacroU" = "MacroU",
                         "Avg_Beta_FinU" = "FinU",
                         "Avg_Beta_TPU" = "TPU",
                         "Avg_Beta_EPU" = "EPU"))

# Step 3: Plot
ggplot(betas_normalized_long, aes(x = Portfolio, y = Normalized_Beta, fill = Factor)) +
  geom_col(position = "dodge") +
  labs(
    title = "Normalized Average Beta Exposure to Uncertainty Factors",
    x = "Portfolio",
    y = "Normalized Beta (per 1 SD of Factor)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "bottom"
  )




