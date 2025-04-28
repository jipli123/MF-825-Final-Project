library(shiny)
library(tidyverse)
library(quadprog)

##### Load data #####
industry_cols <- setdiff(names(full_data_monthly),
                         c("Date", "RMRF", "RF", "MacroU", "TPU", "ShortRate", "EPU", "FinU"))


##### Helper functions #####


estimate_betas_shrunk <- function(full_data_monthly, industry_cols, lookback_months = 24) {
  rolling_betas <- list()
  
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
      y <- window_data[[ind]] - window_data$RF
      x <- window_data %>% select(RMRF, MacroU, FinU, TPU, EPU, ShortRate)
      model <- lm(y ~ ., data = x)
      coefs <- coef(model)
      
      lambda_shrink <- 0.2
      
      coefs_shrunk <- coefs
      coefs_shrunk["RMRF"] <- (1 - lambda_shrink) * coefs["RMRF"] + lambda_shrink * 1
      coefs_shrunk["MacroU"] <- (1 - lambda_shrink) * coefs["MacroU"] + lambda_shrink * 0
      coefs_shrunk["FinU"] <- (1 - lambda_shrink) * coefs["FinU"] + lambda_shrink * 0
      coefs_shrunk["TPU"] <- (1 - lambda_shrink) * coefs["TPU"] + lambda_shrink * 0
      coefs_shrunk["EPU"] <- (1 - lambda_shrink) * coefs["EPU"] + lambda_shrink * 0
      coefs_shrunk["ShortRate"] <- (1 - lambda_shrink) * coefs["ShortRate"] + lambda_shrink * 0
      
      betas_this_date$Beta_RMRF[j] <- coefs_shrunk["RMRF"]
      betas_this_date$Beta_MacroU[j] <- coefs_shrunk["MacroU"]
      betas_this_date$Beta_FinU[j] <- coefs_shrunk["FinU"]
      betas_this_date$Beta_TPU[j] <- coefs_shrunk["TPU"]
      betas_this_date$Beta_EPU[j] <- coefs_shrunk["EPU"]
      betas_this_date$Beta_ShortRate[j] <- coefs_shrunk["ShortRate"]
    }
    rolling_betas[[i - lookback_months]] <- betas_this_date
  }
  
  beta_estimates <- do.call(rbind, rolling_betas)
  return(beta_estimates)
}

optimize_portfolio_soft <- function(Sigma, betas, target_market_beta = 1, gamma_penalty = 10) {
  n <- nrow(betas)
  Sigma_reg <- Sigma + diag(1e-6, n)
  P <- as.matrix(betas[, c("Beta_MacroU", "Beta_FinU", "Beta_TPU", "Beta_EPU", "Beta_ShortRate")])
  PenaltyMat <- gamma_penalty * (P %*% t(P))
  
  Dmat <- Sigma_reg + PenaltyMat
  dvec <- rep(0, n)
  
  Amat <- t(rbind(
    rep(1, n),
    betas$Beta_RMRF
  ))
  
  bvec <- c(1, target_market_beta)
  
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = length(bvec))
  w_opt <- sol$solution
  w_opt <- w_opt / sum(w_opt)
  
  return(w_opt)
}


##### Shiny App #####

ui <- fluidPage(
  titlePanel("Trade Policy Sensitive Portfolio"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("target_beta", "Target Market Beta:", min = 0.5, max = 1.5, value = 1.0, step = 0.05),
      selectInput("crisis_choice", "Select Crisis:", 
                  choices = c("Full Sample", "Dot-Com Crash", "Global Financial Crisis", "European Debt Crisis","China-U.S. Trade War","COVID-19 Crash"))
    ),
    mainPanel(
      plotOutput("cumulative_plot"),
      br(),
      verbatimTextOutput("performance_summary")
    )
  )
)

server <- function(input, output) {
  
  betas_ready <- reactive({
    estimate_betas_shrunk(full_data_monthly, industry_cols)
  })
  
  results_dynamic <- reactive({
    betas_est <- betas_ready()
    
    lookback_months <- 24
    rebalance_every <- 6
    n_assets <- length(industry_cols)
    
    portfolio_returns <- rep(NA, nrow(full_data_monthly))
    portfolio_weights <- matrix(NA, nrow = nrow(full_data_monthly), ncol = n_assets)
    colnames(portfolio_weights) <- industry_cols
    
    for (t in (lookback_months + 1):nrow(full_data_monthly)) {
      
      if ((t - lookback_months - 1) %% rebalance_every != 0 && t > (lookback_months + 1)) {
        portfolio_weights[t, ] <- portfolio_weights[t - 1, ]
        portfolio_returns[t] <- sum(portfolio_weights[t, ] * full_data_monthly[t, industry_cols], na.rm = TRUE)
        next
      }
      
      betas_this_date <- betas_est %>%
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
                                       target_market_beta = input$target_beta,
                                       gamma_penalty = 10)
      
      portfolio_weights[t, ] <- w_opt
      portfolio_returns[t] <- sum(w_opt * full_data_monthly[t, industry_cols], na.rm = TRUE)
    }
    
    results_df <- data.frame(
      Date = full_data_monthly$Date,
      Portfolio_Return = portfolio_returns,
      Market_Return = full_data_monthly$RMRF + full_data_monthly$RF,
      RF = full_data_monthly$RF
    ) %>%
      filter(!is.na(Portfolio_Return)) %>%
      mutate(
        Portfolio_Cumulative = cumprod(1 + Portfolio_Return),
        Market_Cumulative = cumprod(1 + Market_Return),
        Portfolio_Excess = Portfolio_Return - RF,
        Market_Excess = Market_Return - RF
      )
    
    return(results_df)
  })
  
  output$cumulative_plot <- renderPlot({
    results <- results_dynamic()
    
    results_filtered <- results
    title_suffix <- "Full Sample (1990-2024)"
    
    if (input$crisis_choice == "Dot-Com Crash") {
      results_filtered <- results %>% filter(Date >= as.Date("2000-03-01") & Date <= as.Date("2002-10-31"))
      title_suffix <- "Dot-Com Crash"
    } else if (input$crisis_choice == "Global Financial Crisis") {
      results_filtered <- results %>% filter(Date >= as.Date("2007-09-01") & Date <= as.Date("2009-06-30"))
      title_suffix <- "Global Financial Crisis"
    } else if (input$crisis_choice == "European Debt Crisis") {
      results_filtered <- results %>% filter(Date >= as.Date("2009-10-01") & Date <= as.Date("2012-06-30"))
      title_suffix <- "European Debt Crisis"
    } else if (input$crisis_choice == "COVID-19 Crash") {
      results_filtered <- results %>% filter(Date >= as.Date("2020-02-01") & Date <= as.Date("2020-09-30"))
      title_suffix <- "COVID-19 Crash"
    } else if (input$crisis_choice == "China-U.S. Trade War") {
      results_filtered <- results %>% filter(Date >= as.Date("2018-01-01") & Date <= as.Date("2024-12-31"))
      title_suffix <- "China-U.S. Trade War"
    }
    
    ggplot(results_filtered, aes(x = Date)) +
      geom_line(aes(y = Portfolio_Cumulative, color = "Optimized Portfolio"), linewidth = 1) +
      geom_line(aes(y = Market_Cumulative, color = "Market Portfolio"), linewidth = 1) +
      scale_color_manual(values = c("Optimized Portfolio" = "blue", "Market Portfolio" = "red")) +
      labs(
        title = paste0("Cumulative Returns During ", title_suffix, " (Target Beta = ", input$target_beta, ")"),
        y = "Growth of $1",
        x = "Date",
        color = "Portfolio"
      ) +
      theme_minimal()
  })
  
  output$performance_summary <- renderPrint({
    results <- results_dynamic()
    
    portfolio_annual_return <- mean(results$Portfolio_Return) * 12
    portfolio_annual_vol <- sd(results$Portfolio_Return) * sqrt(12)
    portfolio_sharpe <- mean(results$Portfolio_Excess) / sd(results$Portfolio_Return) * sqrt(12)
    
    market_annual_return <- mean(results$Market_Return) * 12
    market_annual_vol <- sd(results$Market_Return) * sqrt(12)
    market_sharpe <- mean(results$Market_Excess) / sd(results$Market_Return) * sqrt(12)
    
    cat("Optimized Portfolio:\n",
        "Annual Return =", round(portfolio_annual_return, 4),
        "| Annual Vol =", round(portfolio_annual_vol, 4),
        "| Sharpe =", round(portfolio_sharpe, 2), "\n\n",
        
        "Market Portfolio:\n",
        "Annual Return =", round(market_annual_return, 4),
        "| Annual Vol =", round(market_annual_vol, 4),
        "| Sharpe =", round(market_sharpe, 2))
  })
}

##### Run App #####

shinyApp(ui = ui, server = server)
