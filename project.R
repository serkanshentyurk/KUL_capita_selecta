library(ggplot2)
library(ggfortify)
library(forecast)
library(stats)
source('/Users/Serkan/Desktop/Academia/KUL/2023 - 2024/2024 Spring/Capita Selecta in Statistics/Project/arima.sim.freq2.r')
source('/Users/Serkan/Desktop/Academia/KUL/2023 - 2024/2024 Spring/Capita Selecta in Statistics/Project/PredictTS.R')

N = 128
set.seed(0833419)
fs = 4
### Q1
data_output = arima.sim.freq2(n = N, 
                              model = list(Poles = c(0.9*exp(2*pi*1i*0.125),0.8*exp(2*pi*1i*0.25)),  # 45 - 90
                                           Zeros = 0.7*exp(2*pi*1i*0.1875)), # 67.5
                              smodel = list(Poles = 0.95,
                                            Zeros = 0.7*exp(2*pi*1i*0.1), # 36
                                            S = 6), 
                              D = 1, std = 2, fs = fs, unit = "day")
data_output



### Q2
data = data_output$data$TimeSeries
ts_data = ts(data, frequency = 4)

auto_ARIMA = auto.arima(ts_data)
auto_ARIMA  # ARIMA(2,1,2)(0,0,2)[4] with drift 
A = c(1, -auto_ARIMA$coef[1:2])
B = c(1,auto_ARIMA$coef[3:4])
C = c(1,auto_ARIMA$coef[5:6])

poles = data.frame(Poles_fp = 1/polyroot(A), RadiusPoles = abs(1/polyroot(A)), AnglePoles = Arg(1/polyroot(A))/2/pi*fs)
zeros = data.frame(Zeros_fp = 1/polyroot(B), RadiusZeros = abs(1/polyroot(B)), AngleZeros = Arg(1/polyroot(B))/2/pi*fs)
seasonal_zeros = data.frame(Zeros_fp = 1/polyroot(C), RadiusZeros = abs(1/polyroot(C)), AngleZeros = Arg(1/polyroot(C))/2/pi*fs)
poles
zeros
seasonal_zeros

s_period = spectrum(ts_data, method = 'pgram', fast = FALSE, log ='no')
s_period_fitted = spectrum(auto_ARIMA$fitted, method ='pgram', fast = FALSE, log = 'no')
Spectral = data.frame(freq = s_period$freq, Density = s_period$spec, AR = s_period_fitted$spec)
plot(ggplot(Spectral, aes(x=freq)) + geom_line(aes(y=log(Density),colour='Original Data')) + geom_line(aes(y=log(AR), colour = 'Fitted Data')))

plot(ts_data, col = "blue")
lines(auto_ARIMA$fitted, col = "red")
legend("topright", legend = c("Data", "Fitted"), col = c("blue", "red"), lty = 1)


### Q3
m_ARIMA = arima(ts_data, order = c(4,1,2), seasonal = list(order = c(1,1,2), period =6))
m_ARIMA
A = c(1,-m_ARIMA$coef[1:4])
B = c(1, m_ARIMA$coef[5:6])
C = c(1,-m_ARIMA$coef[7:7])
D = c(1, m_ARIMA$coef[8:9])

poles = data.frame(Poles_fp = 1/polyroot(A), RadiusPoles = abs(1/polyroot(A)), AnglePoles = Arg(1/polyroot(A))/2/pi*fs)
zeros = data.frame(Zeros_fp = 1/polyroot(B), RadiusZeros = abs(1/polyroot(B)), AngleZeros = Arg(1/polyroot(B))/2/pi*fs)
seasonal_poles = data.frame(Poles_fp = 1/polyroot(C), RadiusPoles = abs(1/polyroot(C)), AnglePoles = Arg(1/polyroot(C))/2/pi*fs)
seasonal_zeros = data.frame(Zeros_fp = 1/polyroot(D), RadiusZeros = abs(1/polyroot(D)), AngleZeros = Arg(1/polyroot(D))/2/pi*fs)
poles
zeros
seasonal_poles
seasonal_zeros



### Q4
data_cropped = data[1:90]
data_rest = data[91:128]
ts_data_cropped = ts(data_cropped, frequency = fs)

auto_ARIMA_cropped = auto.arima(ts_data_cropped)
auto_ARIMA_cropped
PredictTS(data_cropped, auto_ARIMA_cropped, 100, 38, fs, data)
autoplot(forecast(auto_ARIMA_cropped,h=38)) + 
  autolayer(ts_data, series = "Real Data", color = 'red')



### Q5
n_cropped = 90
plot(autoplot(ts_data_cropped))
ts_data_df1 = diff(ts_data_cropped,1)
plot(autoplot(ts_data_df1))

s_period = spectrum(ts_data_df1, method = 'pgram')
s_ar = spectrum(ts_data_df1, method = 'ar', n.freq=(n_cropped/2+1))
Spectral = data.frame(freq = s_period$freq, Period = s_period$spec, AR = s_ar$spec[2:(n_cropped/2+1)])
plot(ggplot(Spectral, aes(x=freq)) + geom_line(aes(y=Period,colour='periodogram')) + geom_line(aes(y=AR, colour = 'AR = 8')))

acf(ts_data_df1, lag.max = 100)
acf(ts_data_df1, lag.max = 100, type = 'partial')


ts_data_df1_df4 = diff(ts_data_df1, 4)
plot(autoplot(ts_data_df1_df4))

s_period = spectrum(ts_data_df1_df4, method = 'pgram')
s_ar = spectrum(ts_data_df1_df4, method = 'ar', n.freq=(n_cropped/2+1))
Spectral = data.frame(freq = s_period$freq, Period = s_period$spec, AR = s_ar$spec[2:(n_cropped/2+1)])
plot(ggplot(Spectral, aes(x=freq)) + geom_line(aes(y=Period,colour='periodogram')) + geom_line(aes(y=AR, colour = 'AR = 7')))

acf(ts_data_df1_df4, lag.max = 100)
acf(ts_data_df1_df4, lag.max = 100, type = 'partial')

arima_result = arima(ts_data_df1_df4, order = c(4,0,4))
arima_result$x <- ts_data_df1_df4 
arima_result

A = c(1,-arima_result$coef[1:4])
B = c(1, arima_result$coef[5:8])

poles = data.frame(Poles_fp = 1/polyroot(A), RadiusPoles = abs(1/polyroot(A)), AnglePoles = Arg(1/polyroot(A))/2/pi*fs)
zeros = data.frame(Zeros_fp = 1/polyroot(B), RadiusZeros = abs(1/polyroot(B)), AngleZeros = Arg(1/polyroot(B))/2/pi*fs)
poles
zeros

PredictTS(ts_data_df1_df4, arima_result, 100, 38, fs,diff(diff(data,1),4))
autoplot(forecast(arima_result,h=38)) +
  autolayer(diff(diff(ts(data,frequency=4),1),4), series = "Real Data", color = 'red')
