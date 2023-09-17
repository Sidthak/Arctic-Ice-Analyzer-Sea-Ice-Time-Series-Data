setwd("/Users/cdmstudent/Desktop/SeaICE")
library(readr)
library(dplyr)
library(ggplot2)   # For qplot
library(fBasics)
library(lubridate) # for mdy date conversion
library(ggfortify) # for ts autoplot
library(zoo)
library(tseries)
library(forecast)
source("backtest.R")
library(lmtest)
seaice <- read_csv("/Users/cdmstudent/Desktop/SeaICE/seaice.csv")
head(seaice)
seaice <- seaice[c(-5, -6)]
seaice$Date <- as.Date(with(seaice, paste(Year, Month, Day, sep = '-')), "%Y-%m-%d")
seaice <- seaice %>% select(-Year, -Month, -Day)
head(seaice)
south_seaice <- filter(seaice, hemisphere == "south")
head(south_seaice)
# Histogram south_seaice
ggplot(south_seaice, aes(x = Extent)) +
  geom_histogram( fill = "lightblue", color = "black") +
  labs(title = "Histogram of south_seaice Sea Ice Extent", x = "Extent", y = "Frequency")

# Q-Q plot south_seaice
ggplot(south_seaice, aes(sample = Extent)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "Q-Q Plot of south_seaice Sea Ice Extent", x = "Theoretical Quantiles", y = "Sample Quantiles")

# Perform Jarque-Bera test
jb_test <- jarque.bera.test(south_seaice$Extent)

# Print the test results
print(jb_test)


south_seaice_ts = ts(south_seaice $Extent, start=c(1978, 10),end  = c(2019, 5), frequency=12)
autoplot(south_seaice_ts)
#######################
# Divide the time series into three segments

segment1 <- window(south_seaice_ts, start = c(1978, 10), end = c(1986, 12))
segment2 <- window(south_seaice_ts, start = c(1987, 1), end = c(1994, 12))
segment3 <- window(south_seaice_ts, start = c(1995, 1), end = c(2002, 12))
segment4 <- window(south_seaice_ts, start = c(2003, 1), end = c(2010, 12))
segment5 <- window(south_seaice_ts, start = c(2011, 1), end = c(2018, 12))
segment6 <- window(south_seaice_ts, start = c(2019, 1), end = c(2019, 5))

# Calculate the mean and variance for each segment
mean_segment1 <- mean(segment1)
var_segment1 <- var(segment1)

mean_segment2 <- mean(segment2)
var_segment2 <- var(segment2)

mean_segment3 <- mean(segment3)
var_segment3 <- var(segment3)

mean_segment4 <- mean(segment4)
var_segment4 <- var(segment4)

mean_segment5 <- mean(segment5)
var_segment5 <- var(segment5)

mean_segment6 <- mean(segment6)
var_segment6 <- var(segment6)
#################

autoplot(decompose(south_seaice_ts))
acf_result<-acf(south_seaice_ts)

install.packages("forecast")
library(forecast)
############
# Perform seasonal decomposition
decomposed <- stl(south_seaice_ts, s.window = "periodic")

# Plot the seasonal component
plot(decomposed$time.series[, "seasonal"], main = "Seasonal Component")

# Perform Loess (S-PLUS) test
loess_test <- decomposed$time.series[, "seasonal"]
loess_test <- loess_test[!is.na(loess_test)]

# Check for seasonality using Loess (S-PLUS) test
if (length(unique(loess_test)) == 1) {
  cat("No significant seasonality detected using Loess (S-PLUS) test.\n")
} else {
  cat("Significant seasonality detected using Loess (S-PLUS) test.\n")
}



####################################################
# Apply Fourier analysis
fourier_result <- fft(south_seaice_ts)

# Calculate the power spectrum
power_spectrum <- Mod(fourier_result)^2

# Plot the power spectrum
plot(power_spectrum, type = "h", main = "Power Spectrum")

# Identify the dominant frequencies
dominant_frequencies <- sort.list(power_spectrum, decreasing = TRUE)

# Check for seasonality using dominant frequencies
if (length(dominant_frequencies) > 1) {
  cat("Significant seasonality detected using Fourier analysis.\n")
  cat("Dominant frequency/ies:", names(dominant_frequencies)[1], "\n")
} else {
  cat("No significant seasonality detected using Fourier analysis.\n")
}
#######################################
# Access the statistics output

acf_result$acf  # Autocorrelation values
acf_result$lag  # Lag values
acf_result$ci   # Confidence intervals

# Print the complete output
print(acf_result)

pacf_result<-pacf(south_seaice_ts)

print(pacf_result)


#############################################



Box.test(south_seaice_ts,type='Ljung')



# Dickey-Fuller unit root test
adf.test(south_seaice_ts)
# Perform Dickey-Fuller unit root test
adf_result <- adf.test(south_seaice_ts)

# KPSS unit root test
kpss.test(south_seaice_ts)


south_seaice_diff <- diff(south_seaice_ts)

autoplot(south_seaice_diff)
hist(south_seaice_diff, main = "Histogram of South_seaice_diff")

# Plot QQ plot
qqnorm(south_seaice_diff, main = "QQ Plot of South_seaice_diff")
qqline(south_seaice_diff)
acf(south_seaice_diff)
pacf(south_seaice_diff)
eacf(south_seaice_diff)

adf.test(south_seaice_diff)
kpss.test(south_seaice_diff)

#########################double diff 
double_south_seaice_diff <- diff(south_seaice_diff)

third_south_diff <- diff(double_south_seaice_diff)

adf.test(third_south_diff)
kpss.test(third_south_diff)

autoplot(south_seaice_diff)
hist(south_seaice_diff, main = "Histogram of South_seaice_diff")

# Plot QQ plot
qqnorm(south_seaice_diff, main = "QQ Plot of South_seaice_diff")
qqline(south_seaice_diff)
acf(south_seaice_diff)
pacf(south_seaice_diff)
eacf(south_seaice_diff)

adf.test(south_seaice_diff)
kpss.test(south_seaice_diff)

###################################
sqrt_series <- sqrt(south_seaice_diff)

# ADF test after square root transformation
adf_result_sqrt <- adf.test(sqrt_series)
adf_p_value_sqrt <- adf_result_sqrt$p.value

# KPSS test after square root transformation
kpss_result_sqrt <- kpss.test(sqrt_series)
kpss_p_value_sqrt <- kpss_result_sqrt$p.value



model1 = Arima(south_seaice_diff, order=c(1, 0, 3), seasonal = list(order=c(2,0,0), seasonal=12))
summary(model1)
coeftest(model1)
Box.test(model1$residuals, lag=10, type="Ljung")

model2 = Arima(south_seaice_diff, order=c(2, 0, 1), seasonal = list(order=c(1,0,0), seasonal=12))
summary(model2)
coeftest(model2)
plot(model2$residuals)
hist(model2$residuals)
Box.test(model2$residuals, lag=10, type="Ljung")


model_aic <- auto.arima(south_seaice_diff, ic = "aic")
# View the AIC-selected model
print(model_aic)
coeftest(model_aic)
hist(model_aic$residuals)
plot(model_aic$residuals)
Box.test(model_aic$residuals, lag=10, type="Ljung")
# Assuming you have a time series model fitted and stored in 'model'
# Assuming the residuals are stored in a variable 'residuals'
residuals<-model_aic$residuals
# Calculate squared residuals
#squared_residuals <- residuals^2


# Create a time index corresponding to the time points
time_index <- 1:length(residuals)

# Create the scatter plot
plot(time_index, residuals, type = "p", pch = 16, xlab = "Time", ylab = "Residuals", main = "Scatter plot of Residuals")

# Create a lagged version of squared residuals
lagged_squared_residuals <- c(NA, residuals[-length(residuals)])

# Perform the Breusch-Pagan test with the lagged squared residuals as a regressor
bp_test <- bptest(residuals ~ lagged_squared_residuals)

# Print the test results
print(bp_test)
hist(model_aic$residuals)


model_bic <- auto.arima(south_seaice_diff, ic = "bic")
# View the BIC-selected model
print(model_bic)
coeftest(model_bic)
hist(model_bic$residuals)
Box.test(model_bic$residuals, lag=10, type="Ljung")

n <- length(south_seaice_diff)
b1 = backtest(model1,south_seaice_diff, h=1, orig=.8*n)
b2 = backtest(model2,south_seaice_diff, h=1, orig=.8*n)
b3 = backtest(model_aic,south_seaice_diff, h=1, orig=.8*n)
b4=backtest(model_bic,south_seaice_diff, h=1, orig=.8*n)
#################################################
##################################################
forecast_model <- forecast(model2, h = 10)
summary(forecast_model)
plot(forecast_model)


#########################ARIMA without seasonal term trial  


# Fit ARIMA model using AIC without seasonality 
model_aic_noseason<- auto.arima(south_seaice_diff, ic = "aic", seasonal = FALSE)

# Print AIC-selected model
print(model_aic_noseason)
coeftest(model_aic_noseason)

Box.test(model_aic_noseason$residuals, lag = 10, type = "Ljung")

# Fit SARIMA model using BIC without seasonality 
model_bic_noseason <- auto.arima(south_seaice_diff, ic = "bic", seasonal = FALSE)

# Print BIC-selected model 
print(model_bic_noseason)
coeftest(model_bic_noseason)
hist(model_bic_noseason$residuals)
Box.test(model_bic_noseason$residuals, lag = 10, type = "Ljung")


forecast_model <- forecast(model_aic, h = 20)
summary(forecast_model)
plot(forecast_model)
