library(lubridate)
library(zoo)
library(tidyverse)

# =============================
# 1. Load and preprocess data
# =============================

# --- TPU ---
df_tpu <- read.csv("TPU INDEX MONTHLY.csv") %>%
  mutate(Date = as.Date(DATE, origin = "1899-12-30")) %>%
  select(Date, TPU)

# --- EPU ---
df_epu <- read.csv("US POLICY UNCERTAINTY.csv") %>%
  mutate(Date = as.Date(paste(Year, Month, 1, sep = "-"))) %>%
  select(Date, EPU = News_Based_Policy_Uncert_Index)

# --- Macro Uncertainty ---
df_sydney <- read.csv("sydneymacrouncertainty.csv") %>%
  mutate(Date = parse_date_time(Date, orders = "m/Y")) %>%
  select(Date, MacroU = h.1)

# --- Treasury ---
df_treasury <- read.csv("THREEFY1.csv") %>%
  rename(Date = observation_date, Treasury1Y = THREEFY1) %>%
  mutate(Date = as.Date(Date))

# --- Combine macro datasets ---
df_combined <- df_tpu %>%
  inner_join(df_epu, by = "Date") %>%
  inner_join(df_sydney, by = "Date") %>%
  inner_join(df_treasury, by = "Date")

# --- Industry returns ---
df_indus <- read.csv("49_Industry_Portfolios.CSV", skip = 11, nrows = 1100) %>%
  rename(YearMonth = X) %>%
  mutate(Date = as.Date(paste0(YearMonth, "01"), format = "%Y%m%d")) %>%
  select(-starts_with("X"))

# --- Merge all together ---
df_final <- inner_join(df_combined, df_indus, by = "Date")

# =============================
# 2. Rolling regression method
# =============================

window_size <- 36
industry_cols <- colnames(df_final)[6:ncol(df_final)]
rolling_betas_list <- list()

for (ind in industry_cols) {
  df_reg <- df_final %>%
    select(Date, y = !!sym(ind), TPU, EPU, MacroU, Treasury1Y)
  
  roll_betas <- rollapply(
    data = df_reg[, -1],
    width = window_size,
    FUN = function(df) coef(lm(y ~ TPU + EPU + MacroU + Treasury1Y, data = as.data.frame(df))),
    by.column = FALSE,
    align = "right"
  )
  
  beta_df <- as.data.frame(roll_betas)
  beta_df$Date <- df_reg$Date[window_size:nrow(df_reg)]
  beta_df$Industry <- ind
  rolling_betas_list[[ind]] <- beta_df
}

df_rolling_betas <- bind_rows(rolling_betas_list)

# --- Summary stats for rolling ---
rolling_summary <- df_rolling_betas %>%
  group_by(Industry) %>%
  summarise(rolling_abs_beta = mean(abs(TPU), na.rm = TRUE))

# =============================
# 3. Exponentially Weighted Filtering (EWF)
# =============================

lambda <- 0.94
X_full <- df_final %>% select(TPU, EPU, MacroU, Treasury1Y)
X_full <- cbind(1, as.matrix(X_full))
colnames(X_full)[1] <- "Intercept"
n <- nrow(X_full)
p <- ncol(X_full)

all_betas_ewf <- list()

for (ind in industry_cols) {
  y_full <- df_final[[ind]]
  betas_ewf <- matrix(NA, nrow = n, ncol = p)
  colnames(betas_ewf) <- colnames(X_full)
  
  for (t in seq(p + 1, n)) {
    X_t <- X_full[1:t, ]
    y_t <- y_full[1:t]
    weights <- lambda^(rev(seq_len(t)) - 1)
    W <- diag(weights)
    
    XtWX <- t(X_t) %*% W %*% X_t
    XtWy <- t(X_t) %*% W %*% y_t
    
    beta_t <- tryCatch(solve(XtWX, XtWy), error = function(e) rep(NA, p))
    betas_ewf[t, ] <- beta_t
  }
  
  beta_df <- as.data.frame(betas_ewf)
  beta_df$Date <- df_final$Date
  beta_df$Industry <- ind
  all_betas_ewf[[ind]] <- beta_df
}

df_ewf_betas <- bind_rows(all_betas_ewf)

# --- Summary stats for EWF ---
ewf_summary <- df_ewf_betas %>%
  filter(!is.na(TPU)) %>%
  group_by(Industry) %>%
  summarise(ewf_abs_beta = mean(abs(TPU), na.rm = TRUE))

# =============================
# 4. Final Comparison Tables
# =============================

# --- Rolling Top/Bottom 5 ---
cat("\n Rolling Beta: Top 5 Most TPU-Sensitive Sectors\n")
print(rolling_summary %>% arrange(desc(rolling_abs_beta)) %>% head(5))

cat("\n Rolling Beta: Bottom 5 Least TPU-Sensitive Sectors\n")
print(rolling_summary %>% arrange(rolling_abs_beta) %>% head(5))

# --- EWF Top/Bottom 5 ---
cat("\n EWF Beta: Top 5 Most TPU-Sensitive Sectors\n")
print(ewf_summary %>% arrange(desc(ewf_abs_beta)) %>% head(5))

cat("\n EWF Beta: Bottom 5 Least TPU-Sensitive Sectors\n")
print(ewf_summary %>% arrange(ewf_abs_beta) %>% head(5))



# Ensure both summaries exist
# rolling_summary: from rolling betas
# ewf_summary: from exponentially weighted filtering

# --- Merge summaries ---
beta_comparison <- rolling_summary %>%
  inner_join(ewf_summary, by = "Industry")

# --- Optional: Sort by rolling or EWF beta ---
# beta_comparison <- beta_comparison %>% arrange(desc(rolling_abs_beta))
# beta_comparison <- beta_comparison %>% arrange(desc(ewf_abs_beta))

# --- View all betas ---
print(beta_comparison)

# --- Save or view nicely ---
View(beta_comparison)  # Open in RStudio viewer
write.csv(beta_comparison, "TPU_beta_comparison.csv", row.names = FALSE)  # Optional export


library(ggplot2)
library(dplyr)

# Choose the sector you want to inspect (e.g., "Steel")
sector_name <- "Gold"

# --- Rolling TPU betas ---
rolling_plot <- df_rolling_betas %>%
  filter(Industry == sector_name) %>%
  select(Date, TPU) %>%
  mutate(Method = "Rolling")

# --- EWF TPU betas ---
ewf_plot <- df_ewf_betas %>%
  filter(Industry == sector_name) %>%
  select(Date, TPU) %>%
  mutate(Method = "EWF")

# --- Combine for plotting ---
df_beta_plot <- bind_rows(rolling_plot, ewf_plot)

# --- Plot the TPU betas over time ---
ggplot(df_beta_plot, aes(x = Date, y = TPU, color = Method)) +
  geom_line(size = 1) +
  labs(
    title = paste("TPU Beta Over Time -", sector_name),
    subtitle = "Comparison of Rolling vs EWF Estimates",
    x = "Date", y = "TPU Beta"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Rolling" = "blue", "EWF" = "red")) +
  theme(legend.position = "bottom")


