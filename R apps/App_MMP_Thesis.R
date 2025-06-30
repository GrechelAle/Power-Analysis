#Final application GBR 
################################################################################
#Title: Application for simulation based power calculation 
#Author: Grechel Taucare
#Last: Update: 22/06/2025
################################################################################


library(shiny)
library(ggplot2)
library(dplyr)

ui <- navbarPage("MMP Simulation-Based Power Analysis Tool",
                 
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
    
    <p>Simulations are based on a log-linear regression model designed to evaluate temporal variation in chemical concentrations while accounting for seasonal patterns and environmental factors. The model estimates the expected log-transformed concentration 
    log(Y<sub>t</sub>)  as a function of time (year), seasonality (month and month squared), and river flow conditions.</p>
    
    <p><strong>The model used for the simulation is:</strong></p>
    <p><em>log(Y<sub>t</sub>) = β₀ + β₁·Year + β₂·Month + β₃·Month² + β₄·RiverFlow + ε<sub>t</sub></em></p>

    <p><strong>Where:</strong></p>
    <ul>
      <li>log(Y<sub>t</sub>): Outcome concentration on the log scale</li>
      <li>β₀: Intercept</li>
      <li>β₁: Trend effect per unit time</li>
      <li>β₂: Month effect</li>
      <li>β₃: Quadratic month effect</li>
      <li>β₄: River flow effect</li>
      <li>ε<sub>t</sub>: Estimated error</li>
    </ul>
  "),
                          
                          h4("Assumptions"),
                          HTML("
    <p>This tool assumes:</p>
    <ul>
      <li>Estimated error is normally distributed with mean zero and constant variance</li>
      <li>Users can provide reasonable estimates of model coefficients and residual error, typically informed by historical data or regression diagnostics</li>
      <li>Sampling is consistent over time (e.g., same months each year, repeated at the same site)</li>
    </ul>
    
   <p><strong>Note:</strong> If residual variability is uncertain, users are encouraged to test a range of plausible values to understand the sensitivity of results to this input.
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
      <li><strong>Number of Years:</strong> Total duration of the monitoring program. Longer durations generally increase the statistical power to detect small trends.</li>
      <li><strong>Sampling Frequency (Months):</strong> Number of samples per year based on evenly spaced intervals (1 = monthly, 2 = bimonthly, 3 = quarterly, 4 = four-monthly).</li>
      <li><strong>Targeted Trend (β₁):</strong> The assumed or targeted trend on the log scale (e.g., –0.10 for a 10% annual decrease). This represents the magnitude of change the monitoring will aim to detect.</li>
      <li><strong>P-value Threshold:</strong> Significance level used to determine whether a trend is statistically detectable (commonly set at 0.05).</li>
    </ul>

    <h4>Model Parameters</h4>
     <p>
All model coefficients must be entered on the <strong>log scale</strong>. In a log-linear model, each 1-unit increase in a predictor (e.g., time) multiplies the expected value of the outcome (Y) by <code>exp(β)</code>. For small values of β, a quick approximation is:
</p>
<blockquote>
  <code>100 × β ≈ percentage change in the outcome</code>
</blockquote>
<p>
For example, <code>β = 0.05</code> implies a ~5% increase, and <code>β = –0.10</code> implies a ~10% decrease.
</p>
    <ul>
      <li><strong>β₀ (Intercept):</strong> The baseline concentration (on the log scale) when all predictors are zero.</li>
      <li><strong>β₂ (Month Effect):</strong> Coefficient capturing linear seasonal variatio(on the log scale)n.</li>
      <li><strong>β₃ (Month Squared Effect):</strong> Coefficient capturing non-linear (quadratic) seasonal effects (on the log scale).</li>
      <li><strong>β₄ (Flow Effect):</strong> Coefficient representing the influence of river flow on concentrations (on the log scale).</li>
      <li><strong>Mean Flow Value:</strong> The average river flow, reflecting typical site conditions (in the original units).</li>
      <li><strong>Flow Sd:</strong> Standard deviation of river flow, representing natural variability in flow (in original units).</li>
      <li><strong>Residual Standard Error:</strong> Represents the unexplained variability in the model. (on the log scale).</li>
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
</p>")
                  )
                 ),
                 
                 
                 # Run tool tab
                 tabPanel("Run Tool",

  fluidPage(
    sidebarLayout(
      sidebarPanel(
        h3(HTML("Monitoring Design Parameters")),
      numericInput("Nvec", "Number of Years:", 10, min = 1, step = 1),
      selectInput("freq", HTML("Sampling Frequency (Months):"), choices = c(1, 2, 3, 4), selected = 1),
      numericInput("B1_value", "β1 Targeted Trend :", -0.10),
      h3(HTML("Model Parameters")),
      numericInput("B0", "β0 (Intercept):", 2.24),
      numericInput("B2", "β2 (Month Effect):", -1.28),
      numericInput("B3", "β3 (Month Squared Effect):", 0.09),
      numericInput("B4", "β4 (Flow Effect):", 0.015),
      numericInput("flow_mean", "Mean Flow Value:", 111),
      numericInput("flow_sd", "Flow Sd:", 80),
      numericInput("sd", "Residual Error:", 1.36),
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

server <- function(input, output) {
  
  simulation_data <- eventReactive(input$run,{
    set.seed(123)
    N <- input$Nvec
    B1 <- input$B1_value
    f <- as.numeric(input$freq)
    freq <- seq(1, 12, by = f)
    B0 <- input$B0
    B2 <- input$B2
    B3 <- input$B3
    B4 <- input$B4
    sd <- input$sd
    flow_mean <- input$flow_mean
    flow_sd <- input$flow_sd
    
    # Generate the selected simulation only when it's called
    x_month <- rep(freq, times = N)
    x_year <- rep(1:N, each = length(freq))
    x_flow <- rnorm(length(x_year), mean = flow_mean, sd = flow_sd)
    y_det <- B0 + B1 * x_year + B2 * x_month + B3 * x_month^2 + B4 * x_flow
    log_y <- y_det + rnorm(length(x_year), sd = sd)
    y <- exp(log_y)
    
    data.frame(x_year, x_month, x_flow, y)
  })
  
  output$simulationPlot <- renderPlot({
    data <- simulation_data()
    ggplot(data, aes(x = x_year, y = log(y), color = as.factor(x_month))) +
      geom_point() +
      #geom_smooth(method = lm) +
      labs(title = "Example Simulated Dataset",
           x = "Number or years (N)",
           y = "Log Outcome (y)") +
      theme_bw()+
      theme(legend.position = "none")
  })
  
  # Event reactive for the result table, triggered by the button
  simulation_results <- eventReactive(input$run, {
    set.seed(123)
    nsim <- input$nsim
    N <- input$Nvec
    B1 <- input$B1_value
    f <- as.numeric(input$freq)
    freq <- seq(1, 12, by = f)
    B0 <- input$B0
    B2 <- input$B2
    B3 <- input$B3
    B4 <- input$B4
    sd <- input$sd
    threshold <- input$threshold
    flow_mean <- input$flow_mean
    flow_sd <- input$flow_sd
    
    # Simulation loop to generate data and calculate metrics
    pval <- numeric(nsim)
    b1val <- numeric(nsim)
    se <- numeric(nsim)
    coverage <- numeric(nsim)
    withProgress(message = "Running simulations...", value = 0, { 
    for (i in 1:nsim) {
      x_month <- rep(freq, times = N)
      x_year <- rep(1:N, each = length(freq))
      x_flow <- rnorm(length(x_year), mean = flow_mean, sd = flow_sd)
      y_det <- B0 + B1 * x_year + B2 * x_month + B3 * x_month^2 + B4 * x_flow
      log_y <- y_det + rnorm(length(x_year), sd = sd)
      y <- exp(log_y)
      
      model <- lm(log(y) ~ x_year + I(x_month^2) + x_flow)
      
      # Extract metrics
      b1val[i] <- coef(summary(model))["x_year", "Estimate"]
      se[i] <- coef(summary(model))["x_year", "Std. Error"]
      pval[i] <- coef(summary(model))["x_year", "Pr(>|t|)"]
      lower_CI <- b1val[i] - 1.96 * se[i]
      upper_CI <- b1val[i] + 1.96 * se[i]
      coverage[i] <- (B1 >= lower_CI) & (B1 <= upper_CI)
     
       incProgress(1 / nsim)  # Increment progress bar
    }
      
    })
    
    # Compute power, coverage, and bias
    Mean_estimate <- mean(b1val)
    SD_estimate <- mean(se)
    power <- sum(pval < threshold) / nsim
    coverage_prob <- mean(coverage)
    percent_bias <- ((mean(b1val) - B1) / B1) * 100
    
    data.frame(
      Coefficient = "(β1) Year",
      Mean.Estimate = Mean_estimate,
      #SD.Estimate = SD_estimate,
      Power.Scenario = power,
      Percent.Bias = percent_bias
      #Coverage = coverage_prob
    )
  })
  
  # The table now renders only after button press
  output$resultTable <- renderTable({
    simulation_results()
  })
}

shinyApp(ui = ui, server = server)
