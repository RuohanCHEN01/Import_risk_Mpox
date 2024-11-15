######################################################
### Title: Importation Risk of Mpox from Hong Kong to Shenzhen
### @Author: Chen Ruohan (HKU)
### R version 4.3.1 (2023-06-16 ucrt)
### Updated: 15/11/2023
######################################################

# Remove all variables
rm(list = ls())

# Methods
# install.packages("minpack.lm")
library(minpack.lm)
library(ggplot2)
library(dplyr)

### To calculate incidence
# Import the list of total cases
setwd("C:/OOOO/RA-HKU/Import_Risk")
data_hk <- read.csv("hk_incident.csv")
data_shenzhen <- read.csv("Shenzhen_incident.csv")
# View(data)
data_hk <- data.frame(data_hk)
data_shenzen <- data.frame(data_shenzhen)
data <- rbind(data_hk,data_shenzhen)
# Import the country list
country_list <- c("HK","Shenzhen")
# Fit sub-exponential distribution using nlsLM
cha_func <- function(t,a,c,k,m){
  a*(1-c*exp(-k*t))^m
}
# Create a empty dataset to store all incidence data
country_record <- c()
for (i in c(1)){
  # Get the special country data
  print(country_list[i])
  country_data <- data[which(data$location==country_list[i]),]
  # Create a new set to store the total cases and incident data
  cumu <- as.data.frame(country_data[,1:2])
  country_record <- as.data.frame(country_data[,1:2])
  # Preprocess the data format
  t <- 1:length(country_data$date)
  t <- as.numeric(unlist(t))
  y <- country_data$cumu_cases
  y <- as.numeric(unlist(y))
  a0 <- y[1]
  # Input the total cases to build model
  fit <- nlsLM(y ~ cha_func(t,a,c,k,m), start = list(a=80, c=0.4, k=0.02, m=20))
  summary(fit)
  # Parameter 95%CI
  conf_intervals <- coef(fit) + cbind(1.96 * sqrt(diag(vcov(fit))), 0 , -1.96 * sqrt(diag(vcov(fit))))
  matrix <- t(conf_intervals)
  matrix <- as.data.frame(matrix)
  # Matrix
  a_matric <- matrix[1]
  c_matric <- matrix[2]
  k_matric <- matrix[3]
  m_matric <- matrix[4]
  predicted <- predict(fit)
  plot(t,y,main=country_list[i])
  lines(t,y=predicted,col="darkgreen",lwd=3)
  # Calculate the incidence
  for (t in 1:length(country_data$date)) {
    Incident <- a_matric[2,1]*
      m_matric[2,1]*
      (1-c_matric[2,1]*
         exp(-k_matric[2,1]*t))^(m_matric[2,1]-1)*c_matric[2,1]*k_matric[2,1]*exp(-k_matric[2,1]*t)
    print(Incident)
    country_record$Incident_mean[t] <- Incident
  }
  View(country_record)
}
# write.csv(country_record,"HK_incidence.csv")

### To calculate prevalence
data <- read.csv("HK_incidence.csv")
data <- data.frame(data)
data_hk <- data[which(data$location=="HK"),]
N_hk <- 7413070
# HK Prevalence
t <- length(data_hk$Incident_mean)
t
data_hk$Incidence_7<- c(data_hk$Incident_mean[8:t],NA,NA,NA,NA,NA,NA,NA)
data_hk$Incidence_8<- c(data_hk$Incidence_7[2:t],NA)
data_hk$Incidence_9<- c(data_hk$Incidence_8[2:t],NA)
data_hk$Incidence_10<- c(data_hk$Incidence_9[2:t],NA)
data_hk$Incidence_11<- c(data_hk$Incidence_10[2:t],NA)
data_hk$Incidence_12<- c(data_hk$Incidence_11[2:t],NA)
data_hk$Incidence_13<- c(data_hk$Incidence_12[2:t],NA)
data_hk$pre <- data_hk$Incidence_7+data_hk$Incidence_8+data_hk$Incidence_9+data_hk$Incidence_10+data_hk$Incidence_11+data_hk$Incidence_12+data_hk$Incidence_13
data_hk$prevalence <- data_hk$pre/N_sz
prevalence_all <- rbind(data_hk[,c(1,2,12)],data_sz[,c(1,2,12)])
View(prevalence_all)
# write.csv(prevalence_all,"HK_prevalence.csv")

### To calculate risk
# Import the list of total cases
data <- read.csv("HK_prevalence.csv")
# View(data)
data <- as.data.frame(data)
import_case <- data$imported_case
len <- length(import_case)
em_data <- data$date
em_data <- data.frame(em_data)
record <- data$date
record <- data.frame(record)
# Calcuate imported cases 95%CI, within 100 iters
imp_estimate_all_low <- matrix(0, ncol = 100, nrow = len)
imp_estimate_all_high <- matrix(0, ncol = 100, nrow = len)
imp_estimate_all <- matrix(0, ncol = 100, nrow = len)
# 100 times
for (seeds in 1:100) {
  set.seed(seeds) 
  for (i in 1:len){
    import_case <- as.numeric(import_case)
    up_low <- t.test(rpois(200, import_case[i]))[4]$conf.int
    em_data$low[i] <- up_low[1]
    em_data$high[i] <- up_low[2]
    em_data$mean[i] <-import_case[i]
  }
  imp_estimate_all_high[, seeds] <- em_data$high
  imp_estimate_all_low[, seeds] <- em_data$low
  imp_estimate_all[, seeds] <- em_data$mean
}
imp_estimate_all <- as.data.frame(imp_estimate_all)
imp_estimate_all_high <- as.data.frame(imp_estimate_all_high)
imp_estimate_all_low <- as.data.frame(imp_estimate_all_low)
imp_estimate_high <- rowMeans(imp_estimate_all_high, na.rm = TRUE)
imp_estimate_low <- rowMeans(imp_estimate_all_low, na.rm = TRUE)
imp_estimate <- rowMeans(imp_estimate_all, na.rm = TRUE)
# Output
record$low <- imp_estimate_low
record$high <- imp_estimate_high
record$mean <- imp_estimate
# write.csv(record,"HK_risk.csv")

