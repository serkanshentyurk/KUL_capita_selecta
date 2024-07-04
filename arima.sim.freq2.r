arima.sim.freq2 <- function(n, model, smodel, D, SD, std, fs, unit)
  {#This function simulates an SARIMA process - Seasonal Autoregressive 
  #Integrated Moving Average of length n. This simulations is conducted in the
  #frequency domain. Please provide in model a list of Poles and Zeros w.r.t. z^-1.
  #For instance model = list (Poles = c(0.93*exp(2*pi*1i*0.125)),0.85*exp(2*pi*1i*0.25), Zeros = 0.10*exp(2*pi*1i*0.1875))
  #Note that it is not necessary to provide the complex conjugates.
  #In the smodel one can provide a list of Poles and Zeros w.r.t. z^-S with $S$ the seasonal order.
  #For instance smodel = list(Poles = 0.8, Zeros = 0.8*exp(2*pi*1i*0.5) , S = 6)
  #One can provide the order of the differentiation in D and the order of the seasonal Differentiator in SD.
  #Finally the standard deviation of the innovation process is provided in std and the sampling frequency in fs (i.e. the number of measurements per unit).
  #In unit one can specify the time base e.g. "days" which implies that the time series is sampled equidistantly
  #at a rate of fs-measurements per day.
  #Try as an example: arima.sim.freq(n = 128, model = list (Poles = c(0.9*exp(2*pi*1i*0.125),0.8*exp(2*pi*1i*0.25)), Zeros = 0.7*exp(2*pi*1i*0.1875)), std = 2, fs = 4, unit = "day")  #The output is ...
  library(polynom)
  library(ggforce)
  library(ggplot2)
  t <- (0:(n-1))/fs; freqUnit <- paste("Frequency [1/",unit,"]"); unit <- paste("Time [",unit,"]");
  
  if (missing(model)==0){
  if (length(model$Poles)>0){DT_Poles <- unique(c(model$Poles,Conj(model$Poles)))} else {DT_Poles <- NA} 
  if (length(model$Zeros)>0){DT_Zeros <- unique(c(model$Zeros,Conj(model$Zeros)))} else {DT_Zeros <- NA}}
  else {DT_Poles <- NA; DT_Zeros <- NA}
  DT_PolesOriginal <- DT_Poles; DT_ZerosOriginal <- DT_Zeros

  if (missing(smodel)==0){
  DT_SPoles <- NA; DT_SZeros <- NA
  if (length(smodel$Poles)>0){
    for (k in 1:length(smodel$Poles)){
  DT_SPoles <- c(DT_SPoles,abs(smodel$Poles[k])^(1/smodel$S)*exp(1i*(2*pi * (seq(smodel$S)-1)+ Arg(smodel$Poles[k]))/smodel$S))}
  DT_SPoles <- DT_SPoles[-1]}
  if (length(smodel$Zeros)>0){
  for (k in 1:length(smodel$Zeros)){
    DT_SZeros <- c(DT_SZeros,abs(smodel$Zeros[k])^(1/smodel$S)*exp(1i*(2*pi * (seq(smodel$S)-1)+ Arg(smodel$Zeros[k]))/smodel$S))}
  DT_SZeros <- DT_SZeros[-1]} 
  if (sum(is.na(DT_Poles))==1){DT_Poles <- DT_SPoles} else if (sum(is.na(DT_SPoles))==1) {DT_Poles <- DT_Poles} else {DT_Poles <- c(DT_Poles,DT_SPoles)}
  if (sum(is.na(DT_Zeros))==1){DT_Zeros <- DT_SZeros} else if (sum(is.na(DT_SZeros))==1) {DT_Zeros <- DT_Zeros} else {DT_Zeros <- c(DT_Zeros,DT_SZeros)}
  } else {DT_SPoles <- NA; DT_SZeros <- NA}
  DT_Poles <- unique(DT_Poles); DT_Zeros <- unique(DT_Zeros);
  if (sum(is.na(DT_Poles))==0){A <- as.polynomial(rev(poly.calc(DT_Poles))); tau = abs(log(1/10)/log(max(abs(DT_Poles)))/fs);} else {A <- as.polynomial(1); tau <- 1}
  if (sum(is.na(DT_Zeros))==0){B <- as.polynomial(rev(poly.calc(DT_Zeros)))} else {B <- as.polynomial(1)}
  if (tau < 1){tau <- 1}
  z <- exp(-2*pi*1i*(0:(n-1))/n);
  H <- as.function(B)(z)/as.function(A)(z);
  h <- Re(1/n*fft(H, inverse = TRUE));
  data <- data.frame(time = t, freq = (0:(n-1))/n*fs, PSD = 20*log10(abs(H)), FRF = H, ImpRes = h);
  plot_PSD <- ggplot(data,aes(x=freq)) + geom_line(aes(y=PSD), color = "red") + labs(x=freqUnit,y="Magnitude [dB]") + coord_cartesian(xlim = c(0,fs/2))
  plot_ImpRes <- ggplot(data, aes(x=time)) + geom_line(aes(y=ImpRes), color = "red") + labs(x=unit, y="Amplitude") + coord_cartesian(xlim = c(0,max(t)/2))
  u <- rnorm(round(n*tau),0,std); data$innovations <- u[(length(u)-n+1):length(u)]; U <- fft(data$innovations, inverse = FALSE);
  y <- Re(1/n*fft(H*U, inverse = TRUE));
  if (missing(D)==0){
   DT_PolesOriginal = c(DT_Poles,rep(1,D));
   data$TimeSeries <- diffinv(y,differences = D)[-seq(D)]
   }
  else {data$TimeSeries <- y}
  plot_TS <- ggplot(data, aes(x=time)) + geom_line(aes(y=TimeSeries, color = "Time Series")) + labs(x=unit, y="Data")
  ComplexPlane <- data.frame(Reals = Re(c(DT_PolesOriginal,DT_ZerosOriginal, DT_SPoles, DT_SZeros)), Imaginaries = Im(c(DT_PolesOriginal,DT_ZerosOriginal, DT_SPoles, DT_SZeros)), Type = c(rep("Poles",length(DT_PolesOriginal)), rep("Zeros", length(DT_ZerosOriginal)), rep("Seasonal Poles",length(DT_SPoles)),rep("Seasonal Zeros", length(DT_SZeros))));
  if (missing(D)==0){ComplexPlane$Type[ComplexPlane$Reals==1] <- paste("Differentiator (order", D, ")")}
  plot_Configuration <- ggplot(ComplexPlane, aes(x=Reals)) + geom_point(aes(y = Imaginaries, color = Type, shape = Type)) + scale_shape_manual(values=c(4, 4, 1, 1, 4)) + geom_circle(aes(x0=0,y0=0,r=1)) + labs(x = "Real Part", y = "Imaginary Part")
  return(list(PSD = plot_PSD,ImpRes = plot_ImpRes,TS = plot_TS, Configuration = plot_Configuration, data = data, TimeConstant = tau, AR = A, MA = B))
  }