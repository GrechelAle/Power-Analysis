################################################################################
#Title: Application for simulation based power calculation 
#Author: Grechel Taucare
#Last: Update: 30/05/2025
################################################################################

library(shiny)
library(ggplot2)
library(dplyr)
library(lmerTest)


#UI


ui <- navbarPage("ITS Simulation-Based Power Analysis Tool",
                 
                 # About tab
                 tabPanel("Overview",
                          wellPanel(
                            h3("Overview"),
                            HTML("
    <p>Recognising the importance of long-term monitoring programs in supporting environmental and public health decisions, 
    it is essential to ensure that the results they produce are both reliable and actionable. 
    Policymakers and regulators depend on these data to guide decision-making, assess progress toward management goals, and respond to emerging challenges. To meet these needs, monitoring programs must be designed to detect meaningful trends with sufficient statistical power while using resources efficiently.</p>
    
    <p>This tool supports the design and evaluation of monitoring programs by estimating <strong>statistical power</strong>,the probability of detecting a trend if one truly exists, using a <strong>simulation-based approach</strong>. 
    It allows users to test how different combinations of monitoring duration,number of samples, number of sites, and trend magnitude influence the likelihood of detecting a trend.</p>
    
    <p>The tool uses a Monte Carlo simulation framework, where synthetic datasets are generated under realistic conditions and analysed using a statistical model. 
    Power is calculated as the proportion of simulations in which the trend is detected (e.g., 80% power means a significant result was found in 800 out of 1,000 simulations).</p>
    
    <p>Simulations are based on an <strong>Interrupted Time Series (ITS) model</strong> with a <strong>mixed-effects structure</strong>, 
    which is well suited for evaluating the impact of interventions or policies on outcomes over time. 
    The mixed-effects framework accounts for both overall trends (fixed effects) and site-level variability (random effects), 
    making it particularly appropriate for monitoring programs that involve repeated measurements across multiple sites or geographic locations.</p>
    
    <p><strong>The model used for the simulation is:</strong></p>
   <p>
  Y<sub>ij</sub> = β<sub>0</sub> + β<sub>1</sub>Time<sub>ij</sub> + β<sub>2</sub>Level<sub>ij</sub> + β<sub>3</sub>Trend<sub>ij</sub> + u<sub>j</sub> + ε<sub>ij</sub>
</p>

<p><strong>Where:</strong></p>
<ul>
  <li>Y<sub>ij</sub>: Outcome value for observation <em>i</em> at site <em>j</em></li>
  <li>β<sub>0</sub>: Intercept</li>
  <li>β<sub>1</sub>: Pre-intervention trend effect per unit time</li>
  <li>β<sub>2</sub>: Level change effect after intervention (fixed effect)</li>
  <li>β<sub>3</sub>: Trend change effect after the intervention</li>
  <li>u<sub>j</sub>: Random site effect that accounts for baseline differences across sites.</li>
  <li>ε<sub>ij</sub>: Residual error term.</li>
</ul>
  "),
                            
                            h4("Assumptions"),
                            HTML("
    <p>This tool assumes:</p>
    <ul>
  <li>Residual errors are normally distributed with mean zero and constant variance (homoscedasticity).</li>
  <li>Random effects (e.g., site-level intercepts) are normally distributed with mean zero and a specified standard deviation.</li>
  <li>Model coefficients, residual standard error, and random effect variance can be reasonably estimated by the user, typically based on historical data, published values, or regression diagnostics.</li>
  <li>Sampling intervals (e.g., time between observations) are consistent over the monitoring period.</li>
  <li>The intervention occurs at a clearly defined time point that applies uniformly across all sites.</li>
  <li>No other major changes happened at the same time as the intervention that could affect the results, unless those changes are included in the model.</li>
</ul>
    
   <p><strong>Note:</strong> If residual and random variability are uncertain, users are encouraged to test a range of plausible values to understand the sensitivity of results to this input.
    If your system involves non-linear patterns, irregular sampling, or other sources of complexity not captured by this model, a custom simulation approach may be more appropriate.</p>
  ")
            )
                          
                 ),
                 # Input values tab
                 tabPanel("Input values",
                          wellPanel(
                            h3("How to input values"),
                            HTML("
    <p>The input values are grouped into two categories: <strong>monitoring design parameters</strong> and <strong>model parameters</strong> which are used in the simulation-based power analysis. 
    The values entered here determine how the tool generates synthetic datasets and estimates the probability of detecting a specified trend over time.</p>
    
    <h4>Monitoring Design Parameters</h4>
    <ul>
      <li><strong>Number of Time Points:</strong> Total number of time points (e.g., years) in the monitoring program. Longer durations generally increase the statistical power to detect small trends.</li>
      <li><strong>Proportion of Time Before Intervention:</strong>  The proportion of the total monitoring duration that occurs before the intervention takes place. For example, if the intervention occurs after 6 years of a 10-year study, the value would be 0.6.</li>
      <li><strong>Number of Sites:</strong> Number of monitoring sites.</li>
      <li><strong>Number of Samples per Site:</strong> Number of samples collected at each site per time point.</li>
      <li><strong>β₂ (Targeted Level Change):</strong> The targeted immediate change in the outcome after the intervention, expressed as percentage of the baseline level (e.g., –15 for a 50% decrease).</li>
      <li><strong>β₃ (Targeted Trend Change):</strong> The targeted change after the intervention, expressed as the percentage of the pre-intervention trend (e.g., –0.50 for a 50% decrease in slope).</li>
    </ul>

    <h4>Model Parameters</h4>
    <ul>
      <li><strong>β₀ (Intercept):</strong> The baseline concentration when all predictors ar set to zero.</li>
      <li><strong>β₁ (Trend Before the interventio):</strong>Coefficient representing the linear trend in the outcome before the intervention.</li>
      <li><strong>Ramdom Error:</strong> Variation between monitoring sites that is not explained by fixed effects. This represents random intercepts in the mixed model.</li>
      <li><strong>Residual Standard Error:</strong> Unexplained variability in measurements within each site over time. This represents individual-level or observation-level error.</li>
     <li><strong>P-value Threshold:</strong> Significance level used to determine whether a trend is statistically detectable (commonly set at 0.05).</li>

    </ul>

<h4>Number of Simulations</h4>
<p>
  This defines how many synthetic datasets will be generated to estimate statistical power.
  Larger values provide more stable and reliable results but require longer computation time.
  A minimum of <strong>1,000 simulations</strong> is recommended.
</p>

<h4>Running the Analysis</h4>
<p>
  After specifying the input values and selecting the number of simulations, click the <strong>Power Analysis</strong> button. The app will:
</p>
<ol>
  <li>Simulate monitoring data based on the selected design,</li>
  <li>Fit the statistical model to each simulated dataset,</li>
  <li>Calculate statistical power as the proportion of simulations in which the targeted trend is statistically significant.</li>
</ol>

<p><strong>Interpretation Example</strong>:<br>
  If 1,000 simulations are run and 800 of them result in a statistically significant trend (p-value below the threshold), the estimated power is <strong>80%</strong>.
  This means the monitoring design has an 80% chance of detecting the specified trend if it truly exists.
</p>
")
      )
      ),
          
                 # Run tool tab
                 tabPanel("Run Tool",
                          
                          fluidPage(
                            
  sidebarLayout(
    sidebarPanel(
      h3(HTML("Monitoring Design Parameters")),
      #selectInput("method", HTML("Select Method:"),  choices = c("Trend", "Intervention"),selected = "Trend"),
      numericInput("N", "Number of Time Points:", 20, min = 1, step = 1),
      numericInput("fra_TP", "Proportion of Time Before Intervention", 0.5),
      numericInput("n_sites", "Number of Sites", 4, min = 2, step = 1),
      numericInput("n_per_site", "Number of Samples per Site:", 1, min = 1, step = 1),
      numericInput("B2", "β2 (Targeted Level Change %)", -15),
      numericInput("B3", "β3 (Trageted Trend Change %):", -55),
      h3(HTML("Model Parameters")),
      numericInput("B0", "β0 (Intercept):", 2115),
      numericInput("B1", "β1 (Trend Before Intervention):", 36),
      numericInput("SD_site", "Random Error:", 336.3),
      numericInput("SD_res", "Residual Error:", 386.2),
      numericInput("nsim", "Number of Simulations:", 1000, min = 100, step = 100),
      sliderInput("threshold", "P-value Threshold", min = 0, max = 1, value = 0.05, step = 0.001),
      
      actionButton("run", "Power Analysis")
    ),
    
    mainPanel(
      plotOutput("simulationPlot"),
      tableOutput("resultTable")
    )
  )
)
)
)
#Server
server <- function(input, output) {
  
  simulation_data <- eventReactive(input$run,{
    n_per_site <- input$n_per_site
    #2. Specifying the parameter values 
    B0 <- input$B0 # Intercept 
    B1= input$B1 # Trend before intervention 
    B2 = (input$B2 * input$B0) / 100 # Level change after the intervention 
    B3 = (input$B3 * input$B1) / 100 # Trend change after the intervention 
    SD_site=  input$SD_site # random effect site 
    SD_res=  input$SD_res # residual error 
    method = input$method
    
    # 3. Defining monitoring design parameters 
    
    N = input$N# Number of time points
    n_sites = input$n_sites #number of sites
    fra_TP <- input$fra_TP #  Fraction of time points before intervention
    Int_TP <- round(input$N* input$fra_TP)  # Intervention time point
    
    # Generate the selected simulation only when it's called
    Time <-  1:N
    Site <- letters[1:n_sites]
    
    x_time <- rep(Time,each = n_sites* n_per_site)
    x_site <- factor(rep(letters[1:n_sites],N * n_per_site)) #int handles scnerio up to 26 sites 
    x_level <- factor(ifelse(x_time <= Int_TP, "Before", "After"),levels = c("Before", "After"))
    x_trend <- ifelse(x_time <= Int_TP, 0, (x_time - Int_TP))
    
    df <- data.frame( x_site,x_time, x_level, x_trend)
    
    # Simulate site-level random effects for control group
    Ran_site <- rnorm(n_sites, 0, SD_site)
    eps <- rnorm(length(x_time), 0, SD_res)
    
    # Simulating methamphetamine excretion mg/dat/1000 people
    site_labels <- letters[1:n_sites]  # e.g., "a", "b", ...
    names(Ran_site) <- site_labels     # so you can look them up by site name
    expected <- B0 + Ran_site[as.character(x_site)] + 
      (B1 * x_time) + 
      (B2 * as.numeric(x_level)) + 
      (B3 * x_trend)
    
    df$simulated_values <- expected + eps
    
    df$Int_TP <- Int_TP
    
    return(df)
  })
  
  output$simulationPlot <- renderPlot({
    data <- simulation_data()
    ggplot(data, aes(x = x_time, y = simulated_values, color = x_site)) +
      geom_point(alpha = 0.7) +
      labs(title = "Example Simulated Data", 
           x = "Year", 
           y = "Outcome Value",
           color = "Sites") +
      theme_minimal()+
      geom_vline(xintercept = data$Int_TP, linetype = "dashed") +
      theme_bw()
  })
  

  
  
  # Event reactive for the result table, triggered by the button
  
  simulation_results <- eventReactive(input$run,{
    set.seed(123)
    nsim <- input$nsim
    threshold <- input$threshold
    n_per_site <- input$n_per_site
    
    # 2. Specifying the parameter values 
    B0 <- input$B0
    B1 <- input$B1
    B2 = (input$B2 * input$B0) / 100 # Level change after the intervention 
    B3 = (input$B3 * input$B1) / 100 # Trend change after the intervention
    SD_site <- input$SD_site
    SD_res <- input$SD_res
    method <- input$method
    
    # 3. Defining monitoring design parameters 
    N <- input$N
    n_sites <- input$n_sites
    fra_TP <- input$fra_TP
    Int_TP <- round(input$N * input$fra_TP)  
    
    # 4. Coefficients to extract and their display names
    coefficients <- c("x_levelAfter", "x_trend")
    display_names <- c("β2 (Level Effect) ", "β3 (Trend Effect)")
    true_values <- c(B2, B3)
    names(true_values) <- coefficients
    
    # 5. Preallocate storage
    results_matrix <- matrix(NA, nrow = nsim, ncol = 6 * length(coefficients),
                             dimnames = list(NULL, 
                                             as.vector(outer(coefficients, c("Estimate", "SE", "PValue", "Lower_CI", "Upper_CI", "Coverage"), paste, sep = "_"))))
    
    # 6. Simulation loop
    withProgress(message = "Running simulations...", value = 0, {
    for (i in 1:nsim) {
      # Generate the selected simulation only when it's called
      Time <- 1:N
      Site <- letters[1:n_sites]
      
      x_time <- rep(Time, each = n_sites * n_per_site)
      x_site <- factor(rep(letters[1:n_sites], N * n_per_site))
      x_level <- factor(ifelse(x_time <= Int_TP, "Before", "After"), levels = c("Before", "After"))
      x_trend <- ifelse(x_time <= Int_TP, 0, (x_time - Int_TP))
      
      df <- data.frame(x_site, x_time, x_level, x_trend)
      
      # Simulate site-level random effects for control group
      Ran_site <- rnorm(n_sites, 0, SD_site)
      eps <- rnorm(length(x_time), 0, SD_res)
      
      # Simulating methamphetamine excretion mg/dat/1000 people
      site_index <- rep(1:n_sites, times = n_per_site * N)
      expected <- B0 + Ran_site[site_index] + 
        (B1 * x_time) + 
        (B2 * as.numeric(x_level)) + 
        (B3 * x_trend)
      
      df$simulated_values <- expected + eps
      
      # 7. Fitting the regression model 
      model <- lmer(simulated_values ~ x_time + x_level + x_trend + (1 | x_site), data = df)
      model_summary <- coef(summary(model))
      
      # 8. Extract metrics efficiently for each coefficient
      for (coef in coefficients) {
        results_matrix[i, paste0(coef, "_Estimate")] <- model_summary[coef, "Estimate"]
        results_matrix[i, paste0(coef, "_SE")] <- model_summary[coef, "Std. Error"]
        results_matrix[i, paste0(coef, "_PValue")] <- model_summary[coef, "Pr(>|t|)"]
        
        # Confidence intervals
        results_matrix[i, paste0(coef, "_Lower_CI")] <- results_matrix[i, paste0(coef, "_Estimate")] - 1.96 * results_matrix[i, paste0(coef, "_SE")]
        results_matrix[i, paste0(coef, "_Upper_CI")] <- results_matrix[i, paste0(coef, "_Estimate")] + 1.96 * results_matrix[i, paste0(coef, "_SE")]
        
        # Coverage calculation
        true_val <- true_values[coef]
        results_matrix[i, paste0(coef, "_Coverage")] <- (true_val >= results_matrix[i, paste0(coef, "_Lower_CI")]) & 
          (true_val <= results_matrix[i, paste0(coef, "_Upper_CI")])
      }
      
      incProgress(1 / nsim)
    }
      
    })
    
    # 9. Summarize Results for Both Coefficients with Display Names
    final_results <- do.call(rbind, lapply(seq_along(coefficients), function(idx) {
      coef <- coefficients[idx]
      display_name <- display_names[idx]
      
      Mean_estimate <- mean(results_matrix[, paste0(coef, "_Estimate")], na.rm = TRUE)
      SD_estimate <- mean(results_matrix[, paste0(coef, "_SE")], na.rm = TRUE)
      power <- mean(results_matrix[, paste0(coef, "_PValue")] < threshold, na.rm = TRUE)
      coverage_prob <- mean(results_matrix[, paste0(coef, "_Coverage")], na.rm = TRUE)
      percent_bias <- ((Mean_estimate - true_values[coef]) / true_values[coef]) * 100
      
      data.frame(
        Coefficient = display_name,
        Mean.Estimate = Mean_estimate,
        SD.Estimate = SD_estimate,
        Power.Scenario = power,
        Percent.Bias = percent_bias
        #Coverage = coverage_prob
      )
    }))
    
    return(final_results)
  })
  
  
  
  # The table now renders only after button press
  output$resultTable <- renderTable({
    simulation_results()
  })
  
}

shinyApp(ui = ui, server = server)


