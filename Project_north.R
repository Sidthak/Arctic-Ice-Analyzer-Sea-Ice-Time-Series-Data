 library(readr)
 library(dplyr)
 library(ggplot2)   # For qplot
 library(fBasics)
 library(lubridate) # for mdy date conversion
 library(ggfortify) # for ts autoplot
 library(zoo)
 library(tseries)
 seaice <- read_csv("~/Desktop/Week 4/seaice.csv")
head(seaice)
seaice <- seaice[c(-5, -6)]
seaice$Date <- as.Date(with(seaice, paste(Year, Month, Day, sep = '-')), "%Y-%m-%d")
seaice <- seaice %>% select(-Year, -Month, -Day)
head(seaice)
tail(seaice)
north_seaice <- filter(seaice, hemisphere == "north")
head(north_seaice)
# Histogram north_seaice
ggplot(north_seaice, aes(x = Extent)) +
  geom_histogram( fill = "lightblue", color = "black") +
  labs(title = "Histogram of north_seaice Sea Ice Extent", x = "Extent", y = "Frequency")

# Q-Q plot north_seaice
ggplot(north_seaice, aes(sample = Extent)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "Q-Q Plot of north_seaice Sea Ice Extent", x = "Theoretical Quantiles", y = "Sample Quantiles")


# Perform Jarque-Bera test
jb_test <- jarque.bera.test(north_seaice$Extent)

# Print the test results
print(jb_test)


north_seaice_ts = ts(north_seaice $Extent, start=c(1978, 10),end  = c(2019, 5), frequency=12)
autoplot(north_seaice_ts)
autoplot(decompose(north_seaice_ts))
acf(north_seaice_ts)
pacf(north_seaice_ts)
Box.test(north_seaice_ts, type = "Ljung-Box")

# Dickey-Fuller unit root test
adf.test(north_seaice_ts)

# KPSS unit root test
kpss.test(north_seaice_ts)

#########################
north_seaice_diff <- diff(north_seaice_ts)

autoplot(north_seaice_diff)
hist(north_seaice_diff, main = "Histogram of north_seaice_diff")

# Plot QQ plot
qqnorm(north_seaice_diff, main = "QQ Plot of north_seaice_diff")
qqline(north_seaice_diff)
acf(north_seaice_diff)
pacf(north_seaice_diff)
eacf(north_seaice_diff)

adf.test(north_seaice_diff)
kpss.test(north_seaice_diff)

model1 = Arima(north_seaice_diff, order=c(2, 0, 2))
summary(model1)
coeftest(model1)
Box.test(model1$residuals, lag=10, type="Ljung")

model2 = Arima(north_seaice_diff, order=c(2, 0, 1))
summary(model2)
coeftest(model2)
Box.test(model2$residuals, lag=10, type="Ljung")


model_aic <- auto.arima(north_seaice_diff, ic = "aic")
# View the AIC-selected model
print(model_aic)
coeftest(model_aic)
Box.test(model_aic$residuals, lag=10, type="Ljung")

model_bic <- auto.arima(north_seaice_diff, ic = "bic")
# View the BIC-selected model
print(model_bic)
coeftest(model_bic)
Box.test(model_bic$residuals, lag=10, type="Ljung")

n <- length(north_seaice_diff)
b1 = backtest(model1,north_seaice_diff, h=1, orig=.8*n)
b2 = backtest(model2,north_seaice_diff, h=1, orig=.8*n)
b3 = backtest(model_aic,north_seaice_diff, h=1, orig=.8*n)
b4=backtest(model_bic,north_seaice_diff, h=1, orig=.8*n)

###############
forecast_model <- forecast(model2, h = 10)
summary(forecast_model)
plot(forecast_model)


#Tentative SARIMA Model (North pole)

mn = Arima(north_seaice_diff, order = c(2, 0, 0), seasonal = list(order = c(0, 1, 2), period = 4))
mn
coeftest(mn)
acf(mn$residuals)

Box.test(mn$residuals, 4, "Ljung-Box")

fn = forecast(mn, h = 10)
plot(fn)
summary(fn)