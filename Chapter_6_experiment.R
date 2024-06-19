
# Libraries

library(SLBDD)
library(forecast)
# Install odpc if not installed yet
# install.packages("odpc_2.0.5.tar.gz", repos = NULL, type = "source")
library(odpc)
library(MLmetrics)
library(tictoc)

# Load dfmpcN routine 
source("dfmpcN.R")

# Seed

set.seed(100510766)

# Data

data=FREDMDApril19

## Train-test split

train=data[1:705,] 
test=data[706:710,]
Housing=c('HOUST','HOUSTNE','HOUSTMW','HOUSTS','HOUSTW','PERMIT','PERMITNE',
          'PERMITMW','PERMITS','PERMITW')
Housing_index=which(colnames(data)%in%Housing)


y_true=as.matrix(test[,Housing_index])

## Scaling
  
x=as.matrix(train)
xs=scale(x)
means=apply(x[,Housing_index],2,mean)
sds=apply(x[,Housing_index],2,sd)

# DFM

tic()
dfm_model=dfmpc(xs)

t=dim(train)[1]
h=5

r=dfm_model$r
A=dfm_model$L
f=dfm_model$F
a=dfm_model$E
univ_models_f=dfm_model$MarmaF
univ_models_a=dfm_model$MarmaE

Pf_t_plus_h=matrix(NA,nrow=r,ncol=h)
Pa_t_plus_h=matrix(NA,nrow=dim(train)[2],ncol=h)


for(i in 1:r){
  model_f=arima(f[,i],order=c(univ_models_f[i,1],univ_models_f[i,6],univ_models_f[i,2]),
              seasonal=list(order=c(univ_models_f[i,3],univ_models_f[i,7],univ_models_f[i,4]),
                            period=univ_models_f[i,5]))
  Pf_t_plus_h[i,]=forecast(model_f,h=h)$mean
}
for(i in 1:dim(train)[2]){
  model_a=arima(a[,i],order=c(univ_models_a[i,1],univ_models_a[i,6],univ_models_a[i,2]),
                seasonal=list(order=c(univ_models_a[i,3],univ_models_a[i,7],univ_models_a[i,4]),
                              period=univ_models_a[i,5]))
  Pa_t_plus_h[i,]=forecast(model_a,h=h)$mean
}

Py_t_plus_h=t(A%*%(Pf_t_plus_h)+Pa_t_plus_h)

y_pred_dfm=Py_t_plus_h[,Housing_index]
for(i in 1:ncol(y_pred_dfm)){
  y_pred_dfm[,i]=y_pred_dfm[,i]*sds[i]+means[i]
}
time_dfm=toc()$callback_msg

## MAE and MSE
MAE(y_pred_dfm,y_true)
MSE(y_pred_dfm,y_true)

write.csv(y_pred_dfm,file="y_DFM.csv", row.names = FALSE)

# DFMN

tic()
dfmN_model=dfmpcN(xs)

t=dim(train)[1]
h=5

r=dfmN_model$r
A=dfmN_model$L
f=dfmN_model$F
a=dfmN_model$E
univ_models_f=dfmN_model$MarmaF
univ_models_a=dfmN_model$MarmaE

Pf_t_plus_h_N=matrix(NA,nrow=r,ncol=h)
Pa_t_plus_h_N=matrix(NA,nrow=dim(train)[2],ncol=h)


for(i in 1:r){
  model_f=arima(f[,i],order=c(univ_models_f[i,1],univ_models_f[i,6],univ_models_f[i,2]),
                seasonal=list(order=c(univ_models_f[i,3],univ_models_f[i,7],univ_models_f[i,4]),
                              period=univ_models_f[i,5]))
  Pf_t_plus_h_N[i,]=forecast(model_f,h=h)$mean
}
for(i in 1:dim(train)[2]){
  model_a=arima(a[,i],order=c(univ_models_a[i,1],univ_models_a[i,6],univ_models_a[i,2]),
                seasonal=list(order=c(univ_models_a[i,3],univ_models_a[i,7],univ_models_a[i,4]),
                              period=univ_models_a[i,5]))
  Pa_t_plus_h_N[i,]=forecast(model_a,h=h)$mean
}

Py_t_plus_h_N=t(A%*%(Pf_t_plus_h_N)+Pa_t_plus_h_N)

y_pred_dfmN=Py_t_plus_h_N[,Housing_index]
for(i in 1:ncol(y_pred_dfmN)){
  y_pred_dfmN[,i]=y_pred_dfmN[,i]*sds[i]+means[i]
}
time_dfmN=toc()$callback_msg


## MAE and MSE
MAE(y_pred_dfmN,y_true)
MSE(y_pred_dfmN,y_true)

write.csv(y_pred_dfmN,file="y_DFMN.csv", row.names = FALSE)

# GDFM

tic()
gdfm_model=crit.odpc(xs,ncores = 10)

# Valores predichos Housing con ruido

y=forecast.odpcs(gdfm_model,h=5)

y_pred_gdfm_factor=y[,Housing_index]
for(i in 1:ncol(y_pred_gdfm_factor)){
  y_pred_gdfm_factor[,i]=y_pred_gdfm_factor[,i]*sds[i]+means[i]
}

y_pred_gdfm=y[,Housing_index]+t(Pa_t_plus_h)[,Housing_index]
for(i in 1:ncol(y_pred_gdfm)){
  y_pred_gdfm[,i]=y_pred_gdfm[,i]*sds[i]+means[i]
}

y_pred_gdfmN=y[,Housing_index]+t(Pa_t_plus_h_N)[,Housing_index]
for(i in 1:ncol(y_pred_gdfmN)){
  y_pred_gdfmN[,i]=y_pred_gdfmN[,i]*sds[i]+means[i]
}
time_gdfm=toc()$callback_msg

## MAE and MSE
MAE(y_pred_gdfm_factor,y_true)
MSE(y_pred_gdfm_factor,y_true)

## MAE and MSE
MAE(y_pred_gdfm,y_true)
MSE(y_pred_gdfm,y_true)

## MAE and MSE
MAE(y_pred_gdfmN,y_true)
MSE(y_pred_gdfmN,y_true)


write.csv(y_pred_gdfm_factor,file="y_GDFM_factor.csv", row.names = FALSE)
write.csv(y_pred_gdfm,file="y_GDFM.csv", row.names = FALSE)
write.csv(y_pred_gdfmN,file="y_GDFMN.csv", row.names = FALSE)



