PredictTS <- function(data,model,Sim,h,Fs, data_all){
  library(ggplot2)
  library(ggfortify)
  library(forecast)
  N <- length(data)
  ARIMA_sim <- matrix(NA,2*N,Sim); c <- matrix(NA,Sim,1)
  for (k in 1:Sim){
    ARIMA_sim[,k] <- c(data,simulate(model, nsim= N, future = TRUE)); c[k] <- cor(data,ARIMA_sim[(N+1):(2*N),k])}
    SE_max_cor <- sqrt((1-max(c)^2)/(N-2))
    ARIMA_sim <- matrix(NA,N+h,Sim);
  while (sum(is.na(ARIMA_sim[1,]))>0) {
    temp <- simulate(model, nsim= N, future = TRUE);
    if (cor(data,temp)>max(c) + qnorm(0.01, 0,1)*SE_max_cor){
      ARIMA_sim[,min(which(is.na(ARIMA_sim[1,])))] <- c(data,temp[1:h])}
  }
  # Forecast_data <- data.frame(time = seq(from = 0, to = ((N+h)-1)/Fs, by = 1/Fs), mean = rowMeans(ARIMA_sim), lower = apply(ARIMA_sim,1,quantile, probs = alpha/2), upper = apply(ARIMA_sim,1,quantile,probs = 1-alpha/2)) 
  Forecast_data <- data.frame(time = seq(from = 0, to = ((N+h)-1)/Fs, by = 1/Fs), mean = rowMeans(ARIMA_sim), lower = rowMeans(ARIMA_sim)-1.96*apply(ARIMA_sim,1,sd), upper = rowMeans(ARIMA_sim)+1.96*apply(ARIMA_sim,1,sd))
  Forecast_data$real = data_all
  plot(ggplot(Forecast_data, aes(x=time)) + geom_line(aes(y=mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.3) + geom_line(aes(y=real), color='red'))
  return(as.list(Forecast_data))
  }