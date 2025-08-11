# Create the data
incidents <- c(0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	2,	3,	3,	3,	3,	3,	3,	3,	3,	4,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6)
incidents <- c(0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	1,	1,	1,	2,	3,	3,	4,	5,	5,	5,	6,	6,	6,	6,	7,	8,	8,	8,	10,	12,	12,	12,	12,	12,	13,	13,	15,	15,	15,	16,	18,	19,	20,	20,	20,	20,	21,	22,	24,	24,	26,	26,	26,	27,	27,	27,	28,	30,	30,	30,	30,	31,	31,	33,	33,	33,	33,	33,	33,	33,	34,	34,	34,	34,	34,	34,	34,	35,	35,	35,	35,	35,	35,	35,	36,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37, 37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37,	37)
#setwd("C:/Users/chenr/Desktop")
#data_hk <- read.csv("hk_incident_1.csv")
#incidents <- data_hk$cumu_cases_revise

# Create time variable
time <- 0:(length(incidents)-1)

# Combine into data frame
data <- data.frame(time = time, incidents = incidents)

# Load necessary packages
library(ggplot2)
library(minpack.lm)  # For robust nonlinear fitting

# Plot raw data
ggplot(data, aes(x = time, y = incidents)) + 
  geom_point() +
  labs(title = "Incidents Over Time", x = "Day", y = "Number of Incidents")

# Fit logistic growth model
# Model formula: y = Asym / (1 + exp((xmid - x)/scal))
# Where:
#   Asym = maximum value (carrying capacity)
#   xmid = time at inflection point
#   scal = growth rate

# Estimate starting parameters
Asym_start <- max(data$incidents) * 1.1 # Slightly above observed maximum
xmid_start <- data$time[which.max(diff(data$incidents))]  # Time of maximum growth
scal_start <- 1  # Default starting value

# Fit model
logistic_model <- nlsLM(
  incidents ~ Asym / (1 + exp((xmid - time)/scal)),
  data = data,
  start = list(Asym = Asym_start, xmid = xmid_start, scal = scal_start),
  control = nls.control(maxiter = 500, warnOnly = TRUE))

# View model summary
summary(logistic_model)

AIC(logistic_model)
BIC(logistic_model)
sigma(logistic_model)

library(minpack.lm)
library(ggplot2)
library(dplyr)


#### Period 1 Incidence ####
a<- 6.110298
c<- 0.4703313
k<- 0.07640862
m<- 68655968
a0<- 1.27391
data_incidence_period1 <- c()

for (t in 1:287) {
  Incident <- a*m*(1-c*exp(-k*t))^(m-1)*c*k*exp(-k*t)
  print(Incident)
  data_incidence_period1$Incident_mean[t] <- Incident
}
data_incidence_period1


#### Period 1 Incidence ####
a<- 37.82178
c<- -1.249739
k<- 0.05719832
m<- -10.547625
a0<- 8.134511

data_incidence_period2 <- c()
for (t in 1:133) {
  Incident <- a*m*(1-c*exp(-k*t))^(m-1)*c*k*exp(-k*t)
  print(Incident)
  data_incidence_period2$Incident_mean[t] <- Incident
}
data_incidence_period2


write.csv(data_incidence_period1,"data_incidence_period1.csv")
write.csv(data_incidence_period2,"data_incidence_period2.csv")

#HK_Fit_incidence.csv <- data_incidence_period1.csv+data_incidence_period2.csv

# Load required packages
library(surveillance)
library(tidyverse)
library(ggplot2)

# Load data
case_data <- read.csv("HK_Fit_incidence.csv") %>% 
  as_tibble() %>%
  mutate(dates = as.Date(dates))  

# Create sts object
sts_data <- sts(observed = case_data$reported_cases,
                start = c(as.numeric(format(min(case_data$dates), "%Y")),
                          as.numeric(format(min(case_data$dates), "%j"))),
                frequency = 365)

# Distribution parameters
# Incubation period (log-normal)
incu_meanlog <- 2.001    
incu_sdlog <- 0.234
max_incu <- 9+17  # Extended to cover combined distribution

# Reporting delay (log-normal)
report_meanlog <- 1.6057319
report_sdlog <- 0.5964992
max_report <- 9+17

# Create PMFs for each distribution
incu_pmf <- dlnorm(0:max_incu, meanlog = incu_meanlog, sdlog = incu_sdlog) %>% 
  {./sum(.)}  # Normalize

report_pmf <- dlnorm(0:max_report, meanlog = report_meanlog, sdlog = report_sdlog) %>% 
  {./sum(.)}  # Normalize

# Combine distributions through convolution
combined_pmf <- convolve(incu_pmf, rev(report_pmf), type = "o") %>% 
  {.[1:(max_incu + max_report + 1)]} %>%  # Trim to reasonable length
  {./sum(.)}  # Renormalize

# Run backprojection with combined delay distribution
results <- backprojNP(sts = sts_data,
                      incu.pmf = combined_pmf,  # Using combined distribution
                      iter = 1000,
                      alpha = 0.05)

# Extract and plot results
infection_incidence <- data.frame(
  date = case_data$dates,
  estimated = results@upperbound
  #lower_ci = results@ci[1,,1],
  #upper_ci = results@ci[2,,1]
)
infection_incidence
ggplot(infection_incidence, aes(x = date)) +
  geom_line(aes(y = estimated), color = "blue") +
  labs(title = "Infection Incidence with Combined Delay Distribution",
       x = "Date", y = "Estimated Infections") +
  theme_minimal()
# If ggplot is fail. 
# dev.off() 
# graphics.off() 
# write.csv(infection_incidence, "3_HK_estimated_infections.csv", row.names = FALSE)


# Set working directory and load data
setwd("C:/Users/chenr/Desktop")
dp <- read.csv("3_HK_estimated_infections.csv")
dp <- data.frame(dp)
incidence_data <- dp$estimated


# Parameters
incubation_meanlog <- 2.001
incubation_sdlog <- 0.234
max_days <- 9  # Maximum incubation period to consider

# Example incidence data (replace with your data)

# Calculate the probability that incubation lasts at least 'd' days
prob_still_incubating <- function(d) {
  1 - plnorm(d, meanlog = incubation_meanlog, sdlog = incubation_sdlog)
}

# Calculate daily prevalence
prevalence <- numeric(length(incidence_data))

for (t in 1:length(incidence_data)) {
  # For each day, sum up past cases still in incubation
  lookback_days <- min(t - 1, max_days)  # How far back to check
  days_ago <- 0:lookback_days  # Days since infection (0 = today)
  
  # Probability each past case is still incubating
  weights <- prob_still_incubating(days_ago)
  
  # Sum weighted incidence to get prevalence
  prevalence[t] <- sum(incidence_data[t - days_ago] * weights)
}

# Result
data.frame(
  Day = 1:length(incidence_data),
  Incidence = incidence_data,
  Prevalence = prevalence
)
prevalence
write.csv(prevalence, "prevalence.csv", row.names = FALSE)

# Set seed for reproducibility
set.seed(123)

library(fitdistrplus)

### Import data
setwd("C:/Users/chenr/Desktop")
data <- read.csv("HK_prevalence.csv")
imported_cases <- data$imported_case
len <- length(imported_cases)

### Initialize storage
imp_estimate_all <- matrix(NA, nrow = len, ncol = 100)
imp_estimate_all_high <- matrix(NA, nrow = len, ncol = 100)
imp_estimate_all_low <- matrix(NA, nrow = len, ncol = 100)

### Temporary storage for daily results
em_data <- data.frame(
  low = numeric(len),
  high = numeric(len),
  mean = numeric(len)
)

### Run simulations
for (seed in 1:100) {
  set.seed(seed)
  
  for (i in 1:len) {
    current_lambda <- imported_cases[i]
    
    # Generate Poisson samples
    # poiss_samples <- rpois(200, current_lambda)
    # Size bigger, error smaller
    poiss_samples <- rpois(100, current_lambda)
    sample_sum <- sum(poiss_samples)  # For poisson.test
    # if (sample_sum == 0) {
    #   poiss_ci <- qgamma(c(0.025, 0.975), shape = 0.5, rate = 200)
    # } else {
    #   poiss_ci <- poisson.test(sample_sum, T = 200)$conf.int
    # }
    # Get Poisson CI
    poiss_ci <- poisson.test(sample_sum, T = 200)$conf.int
    poiss_ci <- qgamma(
      c(0.025, 0.975),
      shape = 0.5 + sample_sum,  # Bayesian-ish adjustment for zeros
      rate = length(poiss_samples)
    )
    # Enforce: low ≥ 0, high ≥ mean (for THIS simulation)
    lower_ci <- pmax(poiss_ci[1], 0)
    # upper_ci <- pmax(poiss_ci[2], current_lambda)  # Ensures high ≥ mean
    upper_ci <- poiss_ci[2]
    
    
    
    
    em_data$low[i] <- lower_ci
    em_data$high[i] <- upper_ci
    em_data$mean[i] <- current_lambda
  }
  
  imp_estimate_all_high[, seed] <- em_data$high
  imp_estimate_all_low[, seed] <- em_data$low
  imp_estimate_all[, seed] <- em_data$mean
}

### Average across simulations (no further adjustment needed)
record <- data.frame(
  low = rowMeans(imp_estimate_all_low, na.rm = TRUE),
  high = rowMeans(imp_estimate_all_high, na.rm = TRUE),
  mean = rowMeans(imp_estimate_all, na.rm = TRUE)
)

### Verify: All high ≥ mean (prints FALSE if any violation)
print(any(record$high < record$mean))  # Should be FALSE

### View results
View(record)
sum(record$low)
sum(record$high)
sum(record$mean)
write.csv(record,"2_HK_risk.csv")
