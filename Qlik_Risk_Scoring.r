# -----------------clearing cache & memory----------------------

gc() 

# -------------inculding Package librarys--------------

library(data.table); # data handle
library(TTR); # Data Transformation
library(moments); # Data formats
library(smooth); # Moving Average
library(greybox);
library(Mcomp);
library(reshape2);
library(forecast); # Forecast capabilities
library(dplyr);
library(gmodels); # statistical methods & models
library(jtrans); # Johnson Transformation

# ------------define functions-------------------

riskscore<-function (x, s = 1, type = 1){
  if(type == 1){
    s <- 1 + s
    bias = 1
    score <- 1 / (1 + exp(-s*(x*bias) + s))
    
  } else if (type == 2){
    s <- 1 + s
    bias = -1
    score <- 1 / (1 + exp(-s*(x*bias) + s))
    
  } else if (type == 3){
    score <- 2 / (1 + exp(-s*(x^2))) - 1
  };
  
  return(score)
};

# -------------loading data from Qlik Sense input--------------

Input.Data<-as.data.table(read.csv(q$inputFile,encoding=q$fileEncoding))

# -------------getting additional arguments from Qlik input-variables for R evaluation-------------
Argument.1<-q$Arg1 # KPIName
Argument.2<-q$Arg2 # MovingPeriods
Argument.3<-q$Arg3 # Risk Score Type
Argument.4<-q$Arg4 # Sensitivity Rsik Score

# ------------performing evaluations---------------
# Vorbereitung
Argument.2<-as.numeric(Argument.2)
Argument.3<-as.numeric(Argument.3)
Argument.4<-as.numeric(Argument.4)

# Datensatz für KPI extrahieren
Data.KPI<-as.numeric(gsub(",", ".",unlist(subset(Input.Data,select = Argument.1))))
Data.KPI[is.na(Data.KPI)] <- 0 # NA's durch 0 ersetzen

# Forecast mit Trend-Komponente
fcP<-30 # forecasting Periods
Data.KPI.ts <- ts(Data.KPI, end = length(Data.KPI), frequency=7) # Daten zu Zeitreihe formatieren
Data.KPI.Decomposed<-decompose(Data.KPI.ts, type = "additive", filter = NULL) # Dekomposition der Zeitreihe

# fourier Transformation und Forecast
Data.KPI.z <- fourier(ts(Data.KPI, frequency=30), K=7) # Zeitreihe
Data.KPI.zf <- fourier(ts(Data.KPI, frequency=30), K=7, h=fcP) # Zeitreihe Forecast
Data.KPI.fit <- auto.arima(Data.KPI.ts, xreg=Data.KPI.z, seasonal=FALSE) # Fitting Modell für Froecast
Data.KPI.fc <- forecast(Data.KPI.fit, xreg=Data.KPI.zf, h=fcP) # Arima Forecast mit Fitted Modell
Data.KPI.fc<-c(unname(unlist(Data.KPI),force=TRUE),c(Data.KPI.fc[["mean"]])) # Transformation zu Vector

Data.KPI.fc.js<-round(as.numeric(jtrans(Data.KPI.fc)$transformed),4) # Forecast zu Normalverteilung transformieren
Data.KPI.fc.js[is.na(Data.KPI.fc.js)]<-0 

# Gleitende Durchschnitte
SMA<-round(SMA(ts(Data.KPI.fc), n = Argument.2), 4) # Simple Moving Average
EMA<-round(EMA(ts(Data.KPI.fc), n = Argument.2), 4) # Exponential Moving Average
SMA[is.na(SMA)]<-Data.KPI.fc
EMA[is.na(EMA)]<-Data.KPI.fc

# Dekomposition der Forecast Zeitreihe
Data.KPI.fc.ts<-ts(Data.KPI.fc, end = length(Data.KPI), frequency=7)
Data.KPI.fc.Decomposed<-decompose(Data.KPI.fc.ts,type = "additive", filter = NULL) # Dekomposition

# Abweichung vom gleitenden Durchschnitt
Data.KPI.ABW<-Data.KPI.fc-EMA # Abweichung des aktuellen Wertes (auch Forecast) vom EMA
Data.KPI.ABW.Z<-round(as.numeric(jtrans(Data.KPI.ABW)$transformed),4) # Normalisierung der ABweichung vom EMA
Data.KPI.ABW.Z[is.na(Data.KPI.ABW.Z)]<-0 # NA's mit nullen auffüllen

# Risk Scoring
Data.KPI.RiskScore <- riskscore(x = Data.KPI.ABW.Z, s = Argument.4, type = Argument.3)
Data.KPI.RiskScore.ts <- riskscore(x = Data.KPI.fc.js, s = Argument.4, type = Argument.3)

# --------adding all data components to one data table for saving and exporting------------
Output.Data<-data.table(TrendComponent = Data.KPI.fc.Decomposed$trend,
                        RandomComponent = Data.KPI.Decomposed$random,
                        SeasonalComponent = Data.KPI.Decomposed$season,
                        SimpleMovingAverage = SMA,
                        ExponentialMovingAverage = EMA,
                        TransformedEMA = Data.KPI.fc.js,
                        ForecastedMean = Data.KPI.fc,
                        RiskScoreR = Data.KPI.RiskScore,
                        RiskScoreTS = Data.KPI.RiskScore.ts
)

# -------------writing output CSV for Loading into Qlik Sense---------------
write.csv(Output.Data,file=q$outputFile)

# END