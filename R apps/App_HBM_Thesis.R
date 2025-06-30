################################################################################
#Title: Application for simulation based power calculation 
#Author: Grechel Taucare
#Last: Update: 22/06/2025
################################################################################

library(shiny)
library(ggplot2)
library(dplyr)

# UI
ui <- navbarPage("HBM Simulation-Based Power Analysis Tool",
                 
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
    
    <p>Simulations are based on a log-linear regression model designed to evaluate temporal trends in chemical concentrations while accounting for demographic variability. The model estimates the expected log-transformed concentration
    log(Y<sub>t</sub>) as a function of time, age group, gender, and their interactions.</p>
    
    <p><strong>The model used for the simulation is:</strong></p>
    <p>
  log(Y<sub>t</sub>) = β<sub>0</sub> + β<sub>1</sub>t<sub>t</sub> + ∑<sup>5</sup><sub>j=1</sub>β<sub>2j</sub>(age_group<sub>tj</sub>) + 
  β<sub>3</sub>(gender<sub>t</sub>) + ∑<sup>5</sup><sub>j=1</sub>β<sub>4j</sub>(gender<sub>t</sub> × age_group<sub>tj</sub>) + ε<sub>t</sub>
</p>
<p><strong>Where:</strong></p>
<ul>
  <li>log(Y<sub>t</sub>): Outcome concentration on the log scale</li>
  <li>β<sub>0</sub>: Intercept</li>
  <li>β<sub>1</sub>: Trend effect per unit time</li>
  <li>β<sub>2j</sub>: Age group effect (for j = 1 to 5)</li>
  <li>β<sub>3</sub>: Gender effect</li>
  <li>β<sub>4j</sub>: Interaction effect between gender and each age group, allowing the effect of gender to vary by age group</li>
  <li>ε<sub>t</sub>: Error term</li>
</ul>
                          "),
                            h4("Assumptions"),
                            HTML("
    <p>This tool assumes:</p>
    <ul>
      <li>Residual errors are normally distributed with mean zero and constant variance (homoscedasticity).</li>
      <li>Model coefficients, residual standard error, and random effect variance can be reasonably estimated by the user, typically based on historical data, published values, or regression diagnostics.</li>
      <li>Sampling intervals (e.g., time between observations) are consistent over the monitoring period.</li>
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
                            
    <p>
The inputs are grouped into two categories: <strong>Monitoring Design Parameters</strong> and <strong>Model Parameters</strong>, which together define the simulation setup. These values determine how the tool generates synthetic datasets and estimates the probability of detecting a specified trend over time.
</p>

<h4>Monitoring Design Parameters</h4>
<ul>
  <li><strong>Number of Time Points:</strong> Total number of time points (e.g., years) in the monitoring program. Longer durations generally increase the statistical power to detect small trends.</li>
  
  <li><strong>Number of Age Groups:</strong> Number of distinct age groups (or other categorical subgroups). The system currently supports up to 6 categories.</li>
  
  <li><strong>Gender Categories:</strong> Fixed at 2 (e.g., male and female) in the current version.</li>
  
  <li><strong>Samples per Group:</strong> Number of samples collected <em>per age group, per gender</em>, at each time point.</li>
  
  <li><strong> β₁ (Targeted Trend):</strong> The assumed trend over time, entered on the <strong>log scale</strong>. For interpretation, a value of <code>–0.10</code> corresponds to an approximate 10% annual decrease, while <code>0.05</code> reflects a ~5% increase.</li>
  <li><strong>β₄ (Interaction Effects: Age Group × Gender):</strong> Coefficients representing how the effect of gender varies by age group. Enter one value per age group, starting with <code>0</code> for the reference group (on the log scale).<br>
  <em>Example:</em> <code>0,0.05,0.01,-0.03,-0.02,0.04</code></li>
  
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
  <li><strong>β₀ (Intercept):</strong> The baseline log-transformed outcome when all predictors (time, age group, gender) are at their reference levels.</li>
  
  <li><strong>β₂ (Age Group Effects):</strong> Coefficients for each age group (relative to the reference group). Enter one value per age group, starting with <code>0</code> for the reference group.<br>
  <em>Example:</em> <code>0,1.13,0.17,-0.14,-0.10,0.16</code></li>
  
  <li><strong>β₃ (Gender Effect):</strong> Coefficient for the effect of gender in the reference age group.</li>
  
  
  <li><strong>Residual Standard Error:</strong> Represents the unexplained variability in the model. This value is critical in determining statistical power.</li>
  <li><strong>P-value Threshold:</strong> The significance level for detecting a trend (commonly set at <code>0.05</code>).</li>

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
      numericInput("N_years", "Number of Time Points:", 5, min = 1),
      numericInput("n_age_groups", "Number of Age Groups:", 6, min = 2),
      tags$div(
        tags$label("Number of Genders:", `for` = "n_gender"),
        tags$input(id = "n_gender", type = "number", value = 2, disabled = NA, class = "form-control")
      )
      ,
      numericInput("n", "Samples per Group:", 3, min = 1),
      numericInput("B1", "β1 (Targeted Trend):", 0.05),
      textInput("B4", "β4j (Targeted Interaction Effects: Age Group × Gender):", "0, 0.05, 0.10,0.20, 0.30, 0.40"),
      h3(HTML("Model Parameters")),
      numericInput("B0", "β0 (Intercept):", 0.89),
      textInput("B2", "β2j (Age Group Effect):", "0,1.13,0.17,-0.14,-0.10,0.16"),
      numericInput("B3", "β3 (Gender Effect):", 0.13),
      numericInput("sigma", "Residual Error:", 0.38),
      numericInput("nsim", "Number of Simulations:", 1000, min = 1),
      sliderInput("threshold", "P-value Threshold", min = 0, max = 1, value = 0.05, step = 0.001),
      actionButton("run_simulation", "Power Analysis")
    ),
    mainPanel(
      plotOutput("plot_simulation"),
      tableOutput("resultTable")
    )
  )
)
)
)
# Server
server <- function(input, output) {
  simulation_data <- eventReactive(input$run_simulation, {
    set.seed(123)
    N_years = input$N_years
    n = input$n
    n_age_groups = input$n_age_groups
    n_gender = 2 #input$n_gender I leave two as fixed if have time i will change this to make it more greneral
    B0 = input$B0
    B1 = input$B1
    B3 = input$B3
    B2 = as.numeric(unlist(strsplit(input$B2, ",")))
    B4 = as.numeric(unlist(strsplit(input$B4, ",")))
    sigma = input$sigma
    
    x_year <- rep(seq(1, N_years), each = n_age_groups * n_gender * n)
    
    x_age_group <- rep(rep(1:n_age_groups, each = n_gender * n), times = length(x_year) / (n_age_groups * n_gender * n))
    
    x_gender <- rep(rep(1:n_gender, each = n), times = length(x_year) / (n_gender * n))
    
    x_gender <- as.numeric(x_gender) - 1
    
    Sample_ID <- unlist(lapply(1:(length(x_year) / n), function(x) seq_len(n)))
    
    log_y_det <- B0 + B1 * x_year + B2[x_age_group] + B3 * x_gender + B4[x_age_group] * x_gender
    
    log_y <- log_y_det + rnorm(length(x_year), sd = sigma)
    y <- exp(log_y)
    
    df <- data.frame(Sample_ID, x_year, x_age_group, x_gender, log_y_det, log_y, y) %>%
      mutate(x_age_group = factor(x_age_group, levels = 1:n_age_groups, labels = paste("Age Group", 1:n_age_groups)),
             x_gender = factor(x_gender, labels = c("Female", "Male")))
  })
  
  
  output$plot_simulation <- renderPlot({
    data <- simulation_data()
    ggplot(data, aes(x = x_year, y = y, colour = x_gender)) +
      geom_point() +
      facet_wrap(. ~ x_age_group) +
      theme_bw() +
      labs(title = "Example Simulated Data by Age Group and Gender", 
           x = "Number of Years(N)", 
           y = "Outcome Value (y)",
           color = "Sex")
  })
  
  simulation_results <- eventReactive(input$run_simulation, {
    set.seed(123)
    nsim <- input$nsim
    N_years <- input$N_years
    n <- input$n
    n_age_groups <- input$n_age_groups
    n_gender <- 2 #input$n_gender
    B0 <- input$B0
    B1 <- input$B1
    B3 <- input$B3
    B2 <- as.numeric(unlist(strsplit(input$B2, ",")))
    B4 <- as.numeric(unlist(strsplit(input$B4, ",")))
    sigma <- input$sigma
    threshold <- input$threshold
    
    # Define coefficients to extract
    coef_names <- c("x_year", paste0("x_genderMale:x_age_groupAge Group ", 2:n_age_groups))
    output_coef_names <- c("(β1) Year", paste("(β4) Male: Age Group", 2:n_age_groups))
    # Initialize results matrix
    results_matrix <- matrix(NA, nrow = nsim, ncol = length(coef_names)*3,
                             dimnames = list(NULL, paste0(rep(coef_names, each = 3),
                                                          c("_Estimate", "_SE", "_pvalue"))))
    withProgress(message = 'Running simulations...', value = 0, {
    for (i in 1:nsim) {
      x_year <- rep(seq(1, N_years), each = n_age_groups * n_gender * n)
      x_age_group <- rep(rep(1:n_age_groups, each = n_gender * n), 
                         times = length(x_year) / (n_age_groups * n_gender * n))
      x_gender <- rep(rep(1:n_gender, each = n), times = length(x_year) / (n_gender * n))
      x_gender <- as.numeric(x_gender) - 1
      Sample_ID <- unlist(lapply(1:(length(x_year) / n), function(x) seq_len(n)))
      
      log_y_det <- B0 + B1 * x_year + B2[x_age_group] + B3 * x_gender + B4[x_age_group] * x_gender
      log_y <- log_y_det + rnorm(length(x_year), sd = sigma)
      y <- exp(log_y)
      
      df <- data.frame(Sample_ID, x_year, x_age_group, x_gender, log_y_det, log_y, y) %>%
        mutate(x_age_group = factor(x_age_group, levels = 1:n_age_groups,
                                    labels = paste("Age Group", 1:n_age_groups)),
               x_gender = factor(x_gender, labels = c("Female", "Male")))
      
      model <- lm(log(y) ~ x_year + x_gender * x_age_group, data = df)
      
      model_summary <- coef(summary(model))
      
      # Extract all metrics efficiently
      for (coef in coef_names) {
        results_matrix[i, paste0(coef, "_Estimate")] <- model_summary[coef, "Estimate"]
        results_matrix[i, paste0(coef, "_SE")] <- model_summary[coef, "Std. Error"]
        results_matrix[i, paste0(coef, "_pvalue")] <- model_summary[coef, "Pr(>|t|)"]
      }
      # Increment progress
      incProgress(1 / nsim)
    }
    })
    
    # Compute summary metrics
    reference_coefs <- c(B1, B4[-1]) # B4[-1] excludes first age group (baseline)
    names(reference_coefs) <- coef_names
    
    final_results <- do.call(rbind, lapply(seq_along(coef_names), function(idx) {
      coef <- coef_names[idx]
      estimates <- results_matrix[, paste0(coef, "_Estimate")]
      ses <- results_matrix[, paste0(coef, "_SE")]
      pvalues <- results_matrix[, paste0(coef, "_pvalue")]
      
      coverage <- mean((estimates - 1.96 * ses <= reference_coefs[coef]) &
                         (estimates + 1.96 * ses >= reference_coefs[coef]))
      
      data.frame(
        Coefficient = output_coef_names[idx],
        Mean.Estimate = mean(estimates),
        #SD.Estimate = mean(ses),
        Power.Scenario = mean(pvalues < threshold),
        Percent.Bias = ((mean(estimates) - reference_coefs[coef]) / reference_coefs[coef]) * 100
        #Coverage = coverage
      )
    }))
    
    final_results
  })
  
  # Render the results in Shiny
  output$resultTable <- renderTable({
    req(simulation_results())
    simulation_results()
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
  