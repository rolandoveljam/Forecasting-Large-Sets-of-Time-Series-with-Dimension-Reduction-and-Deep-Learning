####################

library(SLBDD)
library(forecast)
library(caret)
library(ggplot2)
library(tidyr)
library(progress)
library(gridExtra)
library(polynom)
library(ModelMetrics)
library(tictoc)
library(astsa)

####################

ARIMA_lab <- function(vec){
  vec[is.na(vec)]=0
  if(vec[3]==0){
    # AR model
    lab=paste("AR(",vec[1],")","_(",vec[2],",",vec[4],")",sep="")
  }
  else if(vec[1]==0){
    # MA model
    lab=paste("MA(",vec[3],")","_(",vec[2],",",vec[4],")",sep="")
  }
  else{
    # ARIMA model
    lab=paste("ARMA(",vec[1],",",vec[3],")","_(",vec[2],",",vec[4],")",sep="")
  }
  return(lab)
}

coefs_mse <- function(ref_order,pred_order,ref_coefs,pred_coefs){
  # Coefficients do not include intercept (include.mean=FALSE)
  ar_ref=ref_order[1]
  ar_pred=pred_order[1]
  ma_ref=ref_order[3]
  ma_pred=pred_order[3]
  final_coefs=numeric(ar_ref+ma_ref)
  if((ar_ref==ar_pred)&(ma_ref==ma_pred)){
    # Correct prediction
    final_coefs=pred_coefs
  }else{
    if(ma_ref==0){
      # The real model is AR
      if(ar_ref<=ar_pred){
        # Real model has lower or equal AR order => ignore last coefficients if necessary
        final_coefs=pred_coefs[1:ar_ref]
      } else{
        if(ma_pred==0){
          # Real model has higher AR order + no pred MA => include zeros as last coefficients
          final_coefs=c(pred_coefs,rep(0,ar_ref-ar_pred))
        }else{
          # Real model has higher AR order + pred MA => invert MA 
          roots=abs(1/polyroot(c(1,tail(pred_coefs,ma_pred))))
          if(ar_pred==0){
            p=1
          }else{
            p=polynomial(c(1,pred_coefs[1:ar_pred]))
          }
          for(i in 1:ma_pred){
            r=roots[i]
            p=p*polynomial(c(1,r,r^2,r^3,r^4,r^5))
          }
          final_coefs=as.vector(p[2:(ar_ref+1)])
        }
      }
    }
    if(ar_ref==0){
      # The real model is MA
      if(ma_ref<=ma_pred){
        # Real model has lower or equal MA order => ignore last coefficients if necessary
        final_coefs=pred_coefs[(1+ar_pred):(ar_pred+ma_ref)]
      } else{
        if(ar_pred==0){
          # Real model has higher MA order + no pred AR => include zeros as last coefficients
          final_coefs=c(pred_coefs,rep(0,ma_ref-ma_pred))
        }else{
          # Real model has higher MA order + pred AR => invert AR 
          roots=abs(1/polyroot(c(1,pred_coefs[1:ar_pred])))
          if(ma_pred==0){
            p=1
          }else{
            p=polynomial(c(1,tail(pred_coefs,ma_pred)))
          }
          for(i in 1:ar_pred){
            r=roots[i]
            p=p*polynomial(c(1,r,r^2,r^3,r^4,r^5))
          }
          final_coefs=as.vector(p[2:(ma_ref+1)])
        }
      }
    }
    if((ar_ref!=0)&(ma_ref!=0)){
      # The real model is ARMA
      if(ma_pred>0){
        if(ar_pred>0){
          final_coefs=c(pred_coefs[1:min(ar_pred,ar_ref)],rep(0,max(ar_ref-ar_pred,0)),
                      pred_coefs[(ar_pred+1):(ar_pred+min(ma_ref,ma_pred))],rep(0,max(ma_ref-ma_pred,0)))
        }else{
          final_coefs=c(rep(0,ar_ref),pred_coefs[(ar_pred+1):(ar_pred+min(ma_ref,ma_pred))],rep(0,max(ma_ref-ma_pred,0)))
        }
      }else{
        if(ar_pred>0){
          final_coefs=c(pred_coefs[1:min(ar_pred,ar_ref)],rep(0,max(ar_ref-ar_pred,0)),rep(0,ma_ref))
        }else{
          final_coefs=rep(0,ar_ref+ma_ref)
        }
      }
    }
  }
  return(mse(ref_coefs,final_coefs))
}
  

CMfact <- function(obs, ref){
  obs_fact=as.factor(obs)
  ref_fact=as.factor(ref)
  # Levels sorted alphabetically to follow the correct order
  lev=sort(union(obs_fact,ref_fact))
  CM=confusionMatrix(obs_fact,factor(ref_fact,levels = lev))
  return(CM)
}

tableprop <- function(vec){
  t=sort(summary(as.factor(vec))/length(vec), decreasing=T)
  a=min(which(cumsum(t)>=0.8))
  names=names(t[1:a])
  table=data.frame(Model=names,Percentage=as.numeric(t[1:a]))
  rownames(table)=NULL
  return(table)
}

simulation <- function(n, l, ar_coef=numeric(0), ma_coef=numeric(0), d=0, D=0){
  pb <- progress_bar$new(
    format = "[:bar] :percent ETA: :eta",
    total = n
  )
  ref_order=c(length(ar_coef),d,length(ma_coef))
  ref_coefs=c(ar_coef,ma_coef)
  SLBDD <- character(n)
  SLBDD_mse <- numeric(n)
  SLBDD_s2 <- numeric(n)
  forecast <- character(n)
  forecast_mse <- numeric(n)
  forecast_s2 <- numeric(n)
  SLBDD_mae <- numeric(n)
  forecast_mae <- numeric(n)
  if(D==0){
    # No seasonality
    for(i in 1:n){
      # Simulation of the series
      y <- arima.sim(n=l+1,list(order=ref_order,ar=ar_coef,ma=ma_coef))
      series=y[1:l]
      nextstep=y[l+1]
      # Tackle error produced by arimaSpec() 
      # Error in m1$residuals : $ operator is invalid for atomic vectors
      while(inherits(try({SLBDD_order=arimaSpec(series)$order
        SLBDD_model=arima(series,order=SLBDD_order,include.mean = FALSE)},
        silent = T), "try-error")) {
        # It produces an error, we change the series
        y <- arima.sim(n=l+1,list(order=ref_order,ar=ar_coef,ma=ma_coef))
        series=y[1:l]
        nextstep=y[l+1]
      } 
      # Label predicted by the SLBDD function arimaSpec()
      SLBDD[i] <- ARIMA_lab(c(SLBDD_order,0))
      # MSE from the coefficients predicted by SLBDD
      SLBDD_mse[i] <- coefs_mse(ref_order,SLBDD_order,ref_coefs,as.vector(SLBDD_model$coef))
      # s2 from SLBDD
      SLBDD_s2[i] <- SLBDD_model$sigma2
      # MAE for prediction h=1 with SLBDD model 
      SLBDD_mae[i] <- mae(nextstep,forecast(SLBDD_model,h=1)$mean[1])
      # Label predicted by the forecast function auto.arima()
      forecast_model <- auto.arima(series,allowmean=FALSE,allowdrift = FALSE)
      forecast[i] <- ARIMA_lab(c(arimaorder(forecast_model),0))
      # MSE from the coefficients predicted by forecast
      forecast_mse[i] <- coefs_mse(ref_order,arimaorder(forecast_model),ref_coefs,as.vector(forecast_model$coef))
      # s2 from forecast
      forecast_s2[i] <- forecast_model$sigma2
      # MAE for prediction h=1 with forecast model
      forecast_mae[i] <- mae(nextstep,forecast(forecast_model,h=1)$mean[1])
      # Update the progress bar
      pb$tick() 
    }
  }else{
    # Seasonality
    # Change sign to match notation in astsa package
    ma_coef=-ma_coef
    for(i in 1:n){
      # Simulation of the series
      y <- sarima.sim(n=l+1,ar=ar_coef,ma=ma_coef,d=d,D=D,S=12)
      series=y[1:l]
      nextstep=y[l+1]
      # Tackle error produced by sarimaSpec() 
      # Error in m1$residuals : $ operator is invalid for atomic vectors
      while(inherits(try({SLBDD_order=sarimaSpec(series)$order[c(1:3,5)]
      SLBDD_model=arima(series,order=SLBDD_order[1:3],seasonal=list(order=c(0,SLBDD_order[4],0),period=12),include.mean = FALSE)},
      silent = T), "try-error")) {
        # It produces an error, we change the series
        y <- sarima.sim(n=l+1,ar=ar_coef,ma=ma_coef,d=d,D=D,S=12)
        series=y[1:l]
        nextstep=y[l+1]
      } 
      # Label predicted by the SLBDD function arimaSpec()
      SLBDD[i] <- ARIMA_lab(SLBDD_order)
      # MSE from the coefficients predicted by SLBDD
      SLBDD_mse[i] <- coefs_mse(ref_order,SLBDD_order[1:3],ref_coefs,as.vector(SLBDD_model$coef))
      # s2 from SLBDD
      SLBDD_s2[i] <- SLBDD_model$sigma2
      # MAE for prediction h=1 with SLBDD model 
      SLBDD_mae[i] <- mae(nextstep,forecast(SLBDD_model,h=1)$mean[1])
      # Label predicted by the forecast function auto.arima()
      forecast_model <- auto.arima(series,allowmean=FALSE,allowdrift = FALSE)
      forecast[i] <- ARIMA_lab(arimaorder(forecast_model)[c(1:3,5)])
      # MSE from the coefficients predicted by forecast
      forecast_mse[i] <- coefs_mse(ref_order,arimaorder(forecast_model)[1:3],ref_coefs,as.vector(forecast_model$coef))
      # s2 from forecast
      forecast_s2[i] <- forecast_model$sigma2
      # MAE for prediction h=1 with forecast model
      forecast_mae[i] <- mae(nextstep,forecast(forecast_model,h=1)$mean[1])
      # Update the progress bar
      pb$tick() 
    }
  }
  # Return 
  sol=rbind(SLBDD,SLBDD_mse,SLBDD_s2,SLBDD_mae,forecast,forecast_mse,forecast_s2,forecast_mae)
  return(sol)
}

####################

set.seed(100510766)

# Parameters

length_series_vec <- c(50,100,200)
pairs_dD_differences_vec <- pairlist(c(0,0),c(1,0),c(1,1))
n_series <- 200

# Final results array

# SLBDD or forecast
# 5 AR + 5 MA + 3 ARMA
# length series
# lentgh differences
# metric

complete_results=array(NA,dim=c(13,length(length_series_vec),length(pairs_dD_differences_vec),4,2))

dimnames(complete_results)[[1]] <- c("AR(1)", "AR(2)", "AR(3)","AR(4)","AR(5)","MA(1)", 
                                     "MA(2)", "MA(3)","MA(4)","MA(5)",
                                     "ARMA(1,1)", "ARMA(2,1)", "ARMA(1,2)")
dimnames(complete_results)[[2]] <- length_series_vec
dimnames(complete_results)[[3]] <- as.character(pairs_dD_differences_vec)
dimnames(complete_results)[[4]] <- c("Acc","MSE","s2","MAE")
dimnames(complete_results)[[5]] <- c("SLBDD","forecast")

# Start of the simulation

tic()
n_dif=0
for(differences in pairs_dD_differences_vec){
  n_dif=n_dif+1
  d=differences[1]
  D=differences[2]
  n_length=0
  for(length_series in length_series_vec){
    n_length=n_length+1
    
    # Model simulations  
    
    # AR models
    
    # AR(1) with G=0.3
    
    a=simulation(n_series, length_series, ar_coef=c(0.3),d=d,D=D)
    AR1_0.3_SLBDD <- mean(a[1,]==paste("AR(1)","_(",d,",",D,")",sep=""))
    AR1_0.3_SLBDD_mse <- mean(as.numeric(a[2,]))
    AR1_0.3_SLBDD_s2 <- mean(as.numeric(a[3,]))
    AR1_0.3_SLBDD_mae <- mean(as.numeric(a[4,]))
    AR1_0.3_forecast <- mean(a[5,]==paste("AR(1)","_(",d,",",D,")",sep=""))
    AR1_0.3_forecast_mse <- mean(as.numeric(a[6,]))
    AR1_0.3_forecast_s2 <- mean(as.numeric(a[7,]))
    AR1_0.3_forecast_mae <- mean(as.numeric(a[8,]))
    
    # AR(1) with G=0.7
    
    a=simulation(n_series, length_series, ar_coef=c(0.7),d=d,D=D)
    AR1_0.7_SLBDD <- mean(a[1,]==paste("AR(1)","_(",d,",",D,")",sep=""))
    AR1_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    AR1_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    AR1_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    AR1_0.7_forecast <- mean(a[5,]==paste("AR(1)","_(",d,",",D,")",sep=""))
    AR1_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    AR1_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    AR1_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # Joint AR(1) simulation 
    
    complete_results[1,n_length,n_dif,"Acc","SLBDD"]=mean(c(AR1_0.3_SLBDD,AR1_0.7_SLBDD))
    complete_results[1,n_length,n_dif,"MSE","SLBDD"]=mean(c(AR1_0.3_SLBDD_mse,AR1_0.7_SLBDD_mse))
    complete_results[1,n_length,n_dif,"s2","SLBDD"]=mean(c(AR1_0.3_SLBDD_s2,AR1_0.7_SLBDD_s2))
    complete_results[1,n_length,n_dif,"MAE","SLBDD"]=mean(c(AR1_0.3_SLBDD_mae,AR1_0.7_SLBDD_mae))
    
    complete_results[1,n_length,n_dif,"Acc","forecast"]=mean(c(AR1_0.3_forecast,AR1_0.7_forecast))
    complete_results[1,n_length,n_dif,"MSE","forecast"]=mean(c(AR1_0.3_forecast_mse,AR1_0.7_forecast_mse))
    complete_results[1,n_length,n_dif,"s2","forecast"]=mean(c(AR1_0.3_forecast_s2,AR1_0.7_forecast_s2))
    complete_results[1,n_length,n_dif,"MAE","forecast"]=mean(c(AR1_0.3_forecast_mae,AR1_0.7_forecast_mae))
    
    # AR(2) with G={0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)=(1-0.7*B)^2=(1-1.4*B+0.49*B^2)
    # > 1/polyroot(c(1,-1.4,0.49))
    # [1] 0.7-0i 0.7+0i
    
    a=simulation(n_series, length_series, ar_coef=c(1.4, -0.49),d=d,D=D)
    AR2_0.7_0.7_SLBDD <- mean(a[1,]==paste("AR(2)","_(",d,",",D,")",sep=""))
    AR2_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    AR2_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    AR2_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    AR2_0.7_0.7_forecast <- mean(a[5,]==paste("AR(2)","_(",d,",",D,")",sep=""))
    AR2_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    AR2_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    AR2_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # AR(2) with G={0.3, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)=(1-0.3*B)*(1-0.7*B)=(1-1*B+0.21*B^2)
    # > 1/polyroot(c(1,-1,0.21))
    # [1] 0.7-0i 0.3+0i
    
    a=simulation(n_series, length_series, ar_coef=c(1, -0.21),d=d,D=D)
    AR2_0.3_0.7_SLBDD <- mean(a[1,]==paste("AR(2)","_(",d,",",D,")",sep=""))
    AR2_0.3_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    AR2_0.3_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    AR2_0.3_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    AR2_0.3_0.7_forecast <- mean(a[5,]==paste("AR(2)","_(",d,",",D,")",sep=""))
    AR2_0.3_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    AR2_0.3_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    AR2_0.3_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # Joint AR(2) simulation
    
    complete_results[2,n_length,n_dif,"Acc","SLBDD"]=mean(c(AR2_0.7_0.7_SLBDD,AR2_0.3_0.7_SLBDD))
    complete_results[2,n_length,n_dif,"MSE","SLBDD"]=mean(c(AR2_0.7_0.7_SLBDD_mse,AR2_0.3_0.7_SLBDD_mse))
    complete_results[2,n_length,n_dif,"s2","SLBDD"]=mean(c(AR2_0.7_0.7_SLBDD_s2,AR2_0.3_0.7_SLBDD_s2))
    complete_results[2,n_length,n_dif,"MAE","SLBDD"]=mean(c(AR2_0.7_0.7_SLBDD_mae,AR2_0.3_0.7_SLBDD_mae))
    
    complete_results[2,n_length,n_dif,"Acc","forecast"]=mean(c(AR2_0.7_0.7_forecast,AR2_0.3_0.7_forecast))
    complete_results[2,n_length,n_dif,"MSE","forecast"]=mean(c(AR2_0.7_0.7_forecast_mse,AR2_0.3_0.7_forecast_mse))
    complete_results[2,n_length,n_dif,"s2","forecast"]=mean(c(AR2_0.7_0.7_forecast_s2,AR2_0.3_0.7_forecast_s2))
    complete_results[2,n_length,n_dif,"MAE","forecast"]=mean(c(AR2_0.7_0.7_forecast_mae,AR2_0.3_0.7_forecast_mae))
    
    # AR(3) with G={0.7, 0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)*(1-G3*B)=(1-0.7*B)^3=(1-2.1*B+1.47*B^2-0.343*B^3)
    # > 1/polyroot(c(1,-2.1,1.47,-0.343))
    # [1] 0.7+0i 0.7-0i 0.7-0i
    
    a=simulation(n_series, length_series, ar_coef=c(2.1, -1.47, 0.343),d=d,D=D)
    AR3_0.7_0.7_0.7_SLBDD <- mean(a[1,]==paste("AR(3)","_(",d,",",D,")",sep=""))
    AR3_0.7_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    AR3_0.7_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    AR3_0.7_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    AR3_0.7_0.7_0.7_forecast <- mean(a[5,]==paste("AR(3)","_(",d,",",D,")",sep=""))
    AR3_0.7_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    AR3_0.7_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    AR3_0.7_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # AR(3) with G={0.3, 0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)*(1-G3*B)=(1-0.3*B)*(1-0.7*B)^2=(1-1.7*B+0.91*B^2-0.147*B^3)
    # > 1/polyroot(c(1,-1.7,0.91,-0.147))
    # [1] 0.7+0i 0.7-0i 0.3-0i
    
    a=simulation(n_series, length_series, ar_coef=c(1.7,-0.91,0.147),d=d,D=D)
    AR3_0.3_0.7_0.7_SLBDD <- mean(a[1,]==paste("AR(3)","_(",d,",",D,")",sep=""))
    AR3_0.3_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    AR3_0.3_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    AR3_0.3_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    AR3_0.3_0.7_0.7_forecast <- mean(a[5,]==paste("AR(3)","_(",d,",",D,")",sep=""))
    AR3_0.3_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    AR3_0.3_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    AR3_0.3_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # Joint AR(3) simulation 
    
    complete_results[3,n_length,n_dif,"Acc","SLBDD"]=mean(c(AR3_0.7_0.7_0.7_SLBDD,AR3_0.3_0.7_0.7_SLBDD))
    complete_results[3,n_length,n_dif,"MSE","SLBDD"]=mean(c(AR3_0.7_0.7_0.7_SLBDD_mse,AR3_0.3_0.7_0.7_SLBDD_mse))
    complete_results[3,n_length,n_dif,"s2","SLBDD"]=mean(c(AR3_0.7_0.7_0.7_SLBDD_s2,AR3_0.3_0.7_0.7_SLBDD_s2))
    complete_results[3,n_length,n_dif,"MAE","SLBDD"]=mean(c(AR3_0.7_0.7_0.7_SLBDD_mae,AR3_0.3_0.7_0.7_SLBDD_mae))
    
    complete_results[3,n_length,n_dif,"Acc","forecast"]=mean(c(AR3_0.7_0.7_0.7_forecast,AR3_0.3_0.7_0.7_forecast))
    complete_results[3,n_length,n_dif,"MSE","forecast"]=mean(c(AR3_0.7_0.7_0.7_forecast_mse,AR3_0.3_0.7_0.7_forecast_mse))
    complete_results[3,n_length,n_dif,"s2","forecast"]=mean(c(AR3_0.7_0.7_0.7_forecast_s2,AR3_0.3_0.7_0.7_forecast_s2))
    complete_results[3,n_length,n_dif,"MAE","forecast"]=mean(c(AR3_0.7_0.7_0.7_forecast_mae,AR3_0.3_0.7_0.7_forecast_mae))
    
    # AR(4) with G={0.7, 0.7, 0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)*(1-G3*B)*(1-G4*B)=(1-0.7*B)^4=(1-2.8*B+2.94*B^2-1.372*B^3+0.2401*B^4)
    # > 1/polyroot(c(1,-2.8,2.94,-1.372,0.2401))
    # [1] 0.7-0i 0.7-0i 0.7+0i 0.7+0i
    
    a=simulation(n_series, length_series, ar_coef=c(2.8,-2.94,1.372,-0.2401),d=d,D=D)
    AR4_0.7_0.7_0.7_0.7_SLBDD <- mean(a[1,]==paste("AR(4)","_(",d,",",D,")",sep=""))
    AR4_0.7_0.7_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    AR4_0.7_0.7_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    AR4_0.7_0.7_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    AR4_0.7_0.7_0.7_0.7_forecast <- mean(a[5,]==paste("AR(4)","_(",d,",",D,")",sep=""))
    AR4_0.7_0.7_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    AR4_0.7_0.7_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    AR4_0.7_0.7_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # AR(4) with G={0.3, 0.7, 0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)*(1-G3*B)*(1-G4*B)=(1-0.3*B)*(1-0.7*B)^3=(1-2.4*B+2.1*B^2-0.784*B^3+0.1029*B^4)
    # > 1/polyroot(c(1,-2.4,2.1,-0.784,0.1029))
    # [1] 0.7+0i 0.7-0i 0.7-0i 0.3+0i
    
    a=simulation(n_series, length_series, ar_coef=c(2.4,-2.1,0.784,-0.1029),d=d,D=D)
    AR4_0.3_0.7_0.7_0.7_SLBDD <- mean(a[1,]==paste("AR(4)","_(",d,",",D,")",sep=""))
    AR4_0.3_0.7_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    AR4_0.3_0.7_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    AR4_0.3_0.7_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    AR4_0.3_0.7_0.7_0.7_forecast <- mean(a[5,]==paste("AR(4)","_(",d,",",D,")",sep=""))
    AR4_0.3_0.7_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    AR4_0.3_0.7_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    AR4_0.3_0.7_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # Joint AR(4) simulation
    
    complete_results[4,n_length,n_dif,"Acc","SLBDD"]=mean(c(AR4_0.7_0.7_0.7_0.7_SLBDD,AR4_0.3_0.7_0.7_0.7_SLBDD))
    complete_results[4,n_length,n_dif,"MSE","SLBDD"]=mean(c(AR4_0.7_0.7_0.7_0.7_SLBDD_mse,AR4_0.3_0.7_0.7_0.7_SLBDD_mse))
    complete_results[4,n_length,n_dif,"s2","SLBDD"]=mean(c(AR4_0.7_0.7_0.7_0.7_SLBDD_s2,AR4_0.3_0.7_0.7_0.7_SLBDD_s2))
    complete_results[4,n_length,n_dif,"MAE","SLBDD"]=mean(c(AR4_0.7_0.7_0.7_0.7_SLBDD_mae,AR4_0.3_0.7_0.7_0.7_SLBDD_mae))
    
    complete_results[4,n_length,n_dif,"Acc","forecast"]=mean(c(AR4_0.7_0.7_0.7_0.7_forecast,AR4_0.3_0.7_0.7_0.7_forecast))
    complete_results[4,n_length,n_dif,"MSE","forecast"]=mean(c(AR4_0.7_0.7_0.7_0.7_forecast_mse,AR4_0.3_0.7_0.7_0.7_forecast_mse))
    complete_results[4,n_length,n_dif,"s2","forecast"]=mean(c(AR4_0.7_0.7_0.7_0.7_forecast_s2,AR4_0.3_0.7_0.7_0.7_forecast_s2))
    complete_results[4,n_length,n_dif,"MAE","forecast"]=mean(c(AR4_0.7_0.7_0.7_0.7_forecast_mae,AR4_0.3_0.7_0.7_0.7_forecast_mae))
    
    # AR(5) with G={0.7, 0.7, 0.7, 0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)*(1-G3*B)*(1-G4*B)*(1-G5*B)=(1-0.7*B)^5=(1-3.5*B+4.9*B^2-3.43*B^3+1.2005*B^4-0.16807*B^5)
    # > 1/polyroot(c(1,-3.5,4.9,-3.43,1.2005,-0.16807))
    # [1] 0.7+0i 0.7-0i 0.7-0i 0.7+0i 0.7+0i
    
    a=simulation(n_series, length_series, ar_coef=c(3.5,-4.9,3.43,-1.2005,0.16807),d=d,D=D)
    AR5_0.7_0.7_0.7_0.7_0.7_SLBDD <- mean(a[1,]==paste("AR(5)","_(",d,",",D,")",sep=""))
    AR5_0.7_0.7_0.7_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    AR5_0.7_0.7_0.7_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    AR5_0.7_0.7_0.7_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    AR5_0.7_0.7_0.7_0.7_0.7_forecast <- mean(a[5,]==paste("AR(5)","_(",d,",",D,")",sep=""))
    AR5_0.7_0.7_0.7_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    AR5_0.7_0.7_0.7_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    AR5_0.7_0.7_0.7_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # AR(5) with G={0.3, 0.7, 0.7, 0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)*(1-G3*B)*(1-G4*B)*(1-G5*B)=(1-0.3*B)*(1-0.7*B)^4=(1-3.1*B+3.78*B^2-2.254*B^3+0.6517*B^4-0.07203*B^5)
    # > 1/polyroot(c(1,-3.1,3.78,-2.254,0.6517,-0.07203))
    # [1] 0.3+0i 0.7-0i 0.7-0i 0.7+0i 0.7+0i
    
    a=simulation(n_series, length_series, ar_coef=c(3.1,-3.78,2.254,-0.6517,0.07203),d=d,D=D)
    AR5_0.3_0.7_0.7_0.7_0.7_SLBDD <- mean(a[1,]==paste("AR(5)","_(",d,",",D,")",sep=""))
    AR5_0.3_0.7_0.7_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    AR5_0.3_0.7_0.7_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    AR5_0.3_0.7_0.7_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    AR5_0.3_0.7_0.7_0.7_0.7_forecast <- mean(a[5,]==paste("AR(5)","_(",d,",",D,")",sep=""))
    AR5_0.3_0.7_0.7_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    AR5_0.3_0.7_0.7_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    AR5_0.3_0.7_0.7_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # Joint AR(5) simulation and graphics
    
    complete_results[5,n_length,n_dif,"Acc","SLBDD"]=mean(c(AR5_0.7_0.7_0.7_0.7_0.7_SLBDD,AR5_0.3_0.7_0.7_0.7_0.7_SLBDD))
    complete_results[5,n_length,n_dif,"MSE","SLBDD"]=mean(c(AR5_0.7_0.7_0.7_0.7_0.7_SLBDD_mse,AR5_0.3_0.7_0.7_0.7_0.7_SLBDD_mse))
    complete_results[5,n_length,n_dif,"s2","SLBDD"]=mean(c(AR5_0.7_0.7_0.7_0.7_0.7_SLBDD_s2,AR5_0.3_0.7_0.7_0.7_0.7_SLBDD_s2))
    complete_results[5,n_length,n_dif,"MAE","SLBDD"]=mean(c(AR5_0.7_0.7_0.7_0.7_0.7_SLBDD_mae,AR5_0.3_0.7_0.7_0.7_0.7_SLBDD_mae))
    
    complete_results[5,n_length,n_dif,"Acc","forecast"]=mean(c(AR5_0.7_0.7_0.7_0.7_0.7_forecast,AR5_0.3_0.7_0.7_0.7_0.7_forecast))
    complete_results[5,n_length,n_dif,"MSE","forecast"]=mean(c(AR5_0.7_0.7_0.7_0.7_0.7_forecast_mse,AR5_0.3_0.7_0.7_0.7_0.7_forecast_mse))
    complete_results[5,n_length,n_dif,"s2","forecast"]=mean(c(AR5_0.7_0.7_0.7_0.7_0.7_forecast_s2,AR5_0.3_0.7_0.7_0.7_0.7_forecast_s2))
    complete_results[5,n_length,n_dif,"MAE","forecast"]=mean(c(AR5_0.7_0.7_0.7_0.7_0.7_forecast_mae,AR5_0.3_0.7_0.7_0.7_0.7_forecast_mae))
  
    # MA models
    
    # MA(1) with G=0.3
    
    a=simulation(n_series, length_series, ma_coef=c(0.3),d=d,D=D)
    MA1_0.3_SLBDD <- mean(a[1,]==paste("MA(1)","_(",d,",",D,")",sep=""))
    MA1_0.3_SLBDD_mse <- mean(as.numeric(a[2,]))
    MA1_0.3_SLBDD_s2 <- mean(as.numeric(a[3,]))
    MA1_0.3_SLBDD_mae <- mean(as.numeric(a[4,]))
    MA1_0.3_forecast <- mean(a[5,]==paste("MA(1)","_(",d,",",D,")",sep=""))
    MA1_0.3_forecast_mse <- mean(as.numeric(a[6,]))
    MA1_0.3_forecast_s2 <- mean(as.numeric(a[7,]))
    MA1_0.3_forecast_mae <- mean(as.numeric(a[8,]))
    
    # MA(1) with G=0.7
    
    a=simulation(n_series, length_series, ar_coef=c(0.7),d=d,D=D)
    MA1_0.7_SLBDD <- mean(a[1,]==paste("MA(1)","_(",d,",",D,")",sep=""))
    MA1_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    MA1_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    MA1_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    MA1_0.7_forecast <- mean(a[5,]==paste("MA(1)","_(",d,",",D,")",sep=""))
    MA1_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    MA1_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    MA1_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # Joint MA(1) simulation 
    
    complete_results[6,n_length,n_dif,"Acc","SLBDD"]=mean(c(MA1_0.3_SLBDD,MA1_0.7_SLBDD))
    complete_results[6,n_length,n_dif,"MSE","SLBDD"]=mean(c(MA1_0.3_SLBDD_mse,MA1_0.7_SLBDD_mse))
    complete_results[6,n_length,n_dif,"s2","SLBDD"]=mean(c(MA1_0.3_SLBDD_s2,MA1_0.7_SLBDD_s2))
    complete_results[6,n_length,n_dif,"MAE","SLBDD"]=mean(c(MA1_0.3_SLBDD_mae,MA1_0.7_SLBDD_mae))
    
    complete_results[6,n_length,n_dif,"Acc","forecast"]=mean(c(MA1_0.3_forecast,MA1_0.7_forecast))
    complete_results[6,n_length,n_dif,"MSE","forecast"]=mean(c(MA1_0.3_forecast_mse,MA1_0.7_forecast_mse))
    complete_results[6,n_length,n_dif,"s2","forecast"]=mean(c(MA1_0.3_forecast_s2,MA1_0.7_forecast_s2))
    complete_results[6,n_length,n_dif,"MAE","forecast"]=mean(c(MA1_0.3_forecast_mae,MA1_0.7_forecast_mae))
    
    # MA(2) with G={0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)=(1-0.7*B)^2=(1-1.4*B+0.49*B^2)
    # > 1/polyroot(c(1,-1.4,0.49))
    # [1] 0.7-0i 0.7+0i
    
    a=simulation(n_series, length_series, ma_coef=c(1.4, -0.49),d=d,D=D)
    MA2_0.7_0.7_SLBDD <- mean(a[1,]==paste("MA(2)","_(",d,",",D,")",sep=""))
    MA2_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    MA2_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    MA2_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    MA2_0.7_0.7_forecast <- mean(a[5,]==paste("MA(2)","_(",d,",",D,")",sep=""))
    MA2_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    MA2_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    MA2_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # MA(2) with G={0.3, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)=(1-0.3*B)*(1-0.7*B)=(1-1*B+0.21*B^2)
    # > 1/polyroot(c(1,-1,0.21))
    # [1] 0.7-0i 0.3+0i
    
    a=simulation(n_series, length_series, ma_coef=c(1, -0.21),d=d,D=D)
    MA2_0.3_0.7_SLBDD <- mean(a[1,]==paste("MA(2)","_(",d,",",D,")",sep=""))
    MA2_0.3_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    MA2_0.3_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    MA2_0.3_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    MA2_0.3_0.7_forecast <- mean(a[5,]==paste("MA(2)","_(",d,",",D,")",sep=""))
    MA2_0.3_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    MA2_0.3_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    MA2_0.3_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # Joint MA(2) simulation
    
    complete_results[7,n_length,n_dif,"Acc","SLBDD"]=mean(c(MA2_0.7_0.7_SLBDD,MA2_0.3_0.7_SLBDD))
    complete_results[7,n_length,n_dif,"MSE","SLBDD"]=mean(c(MA2_0.7_0.7_SLBDD_mse,MA2_0.3_0.7_SLBDD_mse))
    complete_results[7,n_length,n_dif,"s2","SLBDD"]=mean(c(MA2_0.7_0.7_SLBDD_s2,MA2_0.3_0.7_SLBDD_s2))
    complete_results[7,n_length,n_dif,"MAE","SLBDD"]=mean(c(MA2_0.7_0.7_SLBDD_mae,MA2_0.3_0.7_SLBDD_mae))
    
    complete_results[7,n_length,n_dif,"Acc","forecast"]=mean(c(MA2_0.7_0.7_forecast,MA2_0.3_0.7_forecast))
    complete_results[7,n_length,n_dif,"MSE","forecast"]=mean(c(MA2_0.7_0.7_forecast_mse,MA2_0.3_0.7_forecast_mse))
    complete_results[7,n_length,n_dif,"s2","forecast"]=mean(c(MA2_0.7_0.7_forecast_s2,MA2_0.3_0.7_forecast_s2))
    complete_results[7,n_length,n_dif,"MAE","forecast"]=mean(c(MA2_0.7_0.7_forecast_mae,MA2_0.3_0.7_forecast_mae))
  
    # MA(3) with G={0.7, 0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)*(1-G3*B)=(1-0.7*B)^3=(1-2.1*B+1.47*B^2-0.343*B^3)
    # > 1/polyroot(c(1,-2.1,1.47,-0.343))
    # [1] 0.7+0i 0.7-0i 0.7-0i
    
    a=simulation(n_series, length_series, ma_coef=c(2.1, -1.47, 0.343),d=d,D=D)
    MA3_0.7_0.7_0.7_SLBDD <- mean(a[1,]==paste("MA(3)","_(",d,",",D,")",sep=""))
    MA3_0.7_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    MA3_0.7_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    MA3_0.7_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    MA3_0.7_0.7_0.7_forecast <- mean(a[5,]==paste("MA(3)","_(",d,",",D,")",sep=""))
    MA3_0.7_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    MA3_0.7_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    MA3_0.7_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # MA(3) with G={0.3, 0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)*(1-G3*B)=(1-0.3*B)*(1-0.7*B)^2=(1-1.7*B+0.91*B^2-0.147*B^3)
    # > 1/polyroot(c(1,-1.7,0.91,-0.147))
    # [1] 0.7+0i 0.7-0i 0.3-0i
    
    a=simulation(n_series, length_series, ma_coef=c(1.7,-0.91,0.147),d=d,D=D)
    MA3_0.3_0.7_0.7_SLBDD <- mean(a[1,]==paste("MA(3)","_(",d,",",D,")",sep=""))
    MA3_0.3_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    MA3_0.3_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    MA3_0.3_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    MA3_0.3_0.7_0.7_forecast <- mean(a[5,]==paste("MA(3)","_(",d,",",D,")",sep=""))
    MA3_0.3_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    MA3_0.3_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    MA3_0.3_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # Joint MA(3) simulation 
    
    complete_results[8,n_length,n_dif,"Acc","SLBDD"]=mean(c(MA3_0.7_0.7_0.7_SLBDD,MA3_0.3_0.7_0.7_SLBDD))
    complete_results[8,n_length,n_dif,"MSE","SLBDD"]=mean(c(MA3_0.7_0.7_0.7_SLBDD_mse,MA3_0.3_0.7_0.7_SLBDD_mse))
    complete_results[8,n_length,n_dif,"s2","SLBDD"]=mean(c(MA3_0.7_0.7_0.7_SLBDD_s2,MA3_0.3_0.7_0.7_SLBDD_s2))
    complete_results[8,n_length,n_dif,"MAE","SLBDD"]=mean(c(MA3_0.7_0.7_0.7_SLBDD_mae,MA3_0.3_0.7_0.7_SLBDD_mae))
    
    complete_results[8,n_length,n_dif,"Acc","forecast"]=mean(c(MA3_0.7_0.7_0.7_forecast,MA3_0.3_0.7_0.7_forecast))
    complete_results[8,n_length,n_dif,"MSE","forecast"]=mean(c(MA3_0.7_0.7_0.7_forecast_mse,MA3_0.3_0.7_0.7_forecast_mse))
    complete_results[8,n_length,n_dif,"s2","forecast"]=mean(c(MA3_0.7_0.7_0.7_forecast_s2,MA3_0.3_0.7_0.7_forecast_s2))
    complete_results[8,n_length,n_dif,"MAE","forecast"]=mean(c(MA3_0.7_0.7_0.7_forecast_mae,MA3_0.3_0.7_0.7_forecast_mae))
    
    # MA(4) with G={0.7, 0.7, 0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)*(1-G3*B)*(1-G4*B)=(1-0.7*B)^4=(1-2.8*B+2.94*B^2-1.372*B^3+0.2401*B^4)
    # > 1/polyroot(c(1,-2.8,2.94,-1.372,0.2401))
    # [1] 0.7-0i 0.7-0i 0.7+0i 0.7+0i
    
    a=simulation(n_series, length_series, ma_coef=c(2.8,-2.94,1.372,-0.2401),d=d,D=D)
    MA4_0.7_0.7_0.7_0.7_SLBDD <- mean(a[1,]==paste("MA(4)","_(",d,",",D,")",sep=""))
    MA4_0.7_0.7_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    MA4_0.7_0.7_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    MA4_0.7_0.7_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    MA4_0.7_0.7_0.7_0.7_forecast <- mean(a[5,]==paste("MA(4)","_(",d,",",D,")",sep=""))
    MA4_0.7_0.7_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    MA4_0.7_0.7_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    MA4_0.7_0.7_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # MA(4) with G={0.3, 0.7, 0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)*(1-G3*B)*(1-G4*B)=(1-0.3*B)*(1-0.7*B)^3=(1-2.4*B+2.1*B^2-0.784*B^3+0.1029*B^4)
    # > 1/polyroot(c(1,-2.4,2.1,-0.784,0.1029))
    # [1] 0.7+0i 0.7-0i 0.7-0i 0.3+0i
    
    a=simulation(n_series, length_series, ma_coef=c(2.4,-2.1,0.784,-0.1029),d=d,D=D)
    MA4_0.3_0.7_0.7_0.7_SLBDD <- mean(a[1,]==paste("MA(4)","_(",d,",",D,")",sep=""))
    MA4_0.3_0.7_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    MA4_0.3_0.7_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    MA4_0.3_0.7_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    MA4_0.3_0.7_0.7_0.7_forecast <- mean(a[5,]==paste("MA(4)","_(",d,",",D,")",sep=""))
    MA4_0.3_0.7_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    MA4_0.3_0.7_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    MA4_0.3_0.7_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # Joint MA(4) simulation
    
    complete_results[9,n_length,n_dif,"Acc","SLBDD"]=mean(c(MA4_0.7_0.7_0.7_0.7_SLBDD,MA4_0.3_0.7_0.7_0.7_SLBDD))
    complete_results[9,n_length,n_dif,"MSE","SLBDD"]=mean(c(MA4_0.7_0.7_0.7_0.7_SLBDD_mse,MA4_0.3_0.7_0.7_0.7_SLBDD_mse))
    complete_results[9,n_length,n_dif,"s2","SLBDD"]=mean(c(MA4_0.7_0.7_0.7_0.7_SLBDD_s2,MA4_0.3_0.7_0.7_0.7_SLBDD_s2))
    complete_results[9,n_length,n_dif,"MAE","SLBDD"]=mean(c(MA4_0.7_0.7_0.7_0.7_SLBDD_mae,MA4_0.3_0.7_0.7_0.7_SLBDD_mae))
    
    complete_results[9,n_length,n_dif,"Acc","forecast"]=mean(c(MA4_0.7_0.7_0.7_0.7_forecast,MA4_0.3_0.7_0.7_0.7_forecast))
    complete_results[9,n_length,n_dif,"MSE","forecast"]=mean(c(MA4_0.7_0.7_0.7_0.7_forecast_mse,MA4_0.3_0.7_0.7_0.7_forecast_mse))
    complete_results[9,n_length,n_dif,"s2","forecast"]=mean(c(MA4_0.7_0.7_0.7_0.7_forecast_s2,MA4_0.3_0.7_0.7_0.7_forecast_s2))
    complete_results[9,n_length,n_dif,"MAE","forecast"]=mean(c(MA4_0.7_0.7_0.7_0.7_forecast_mae,MA4_0.3_0.7_0.7_0.7_forecast_mae))
    
    # MA(5) with G={0.7, 0.7, 0.7, 0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)*(1-G3*B)*(1-G4*B)*(1-G5*B)=(1-0.7*B)^5=(1-3.5*B+4.9*B^2-3.43*B^3+1.2005*B^4-0.16807*B^5)
    # > 1/polyroot(c(1,-3.5,4.9,-3.43,1.2005,-0.16807))
    # [1] 0.7+0i 0.7-0i 0.7-0i 0.7+0i 0.7+0i
    
    a=simulation(n_series, length_series, ma_coef=c(3.5,-4.9,3.43,-1.2005,0.16807),d=d,D=D)
    MA5_0.7_0.7_0.7_0.7_0.7_SLBDD <- mean(a[1,]==paste("MA(5)","_(",d,",",D,")",sep=""))
    MA5_0.7_0.7_0.7_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    MA5_0.7_0.7_0.7_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    MA5_0.7_0.7_0.7_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    MA5_0.7_0.7_0.7_0.7_0.7_forecast <- mean(a[5,]==paste("MA(5)","_(",d,",",D,")",sep=""))
    MA5_0.7_0.7_0.7_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    MA5_0.7_0.7_0.7_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    MA5_0.7_0.7_0.7_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # MA(5) with G={0.3, 0.7, 0.7, 0.7, 0.7}
    
    # In this case (1-G1*B)*(1-G2*B)*(1-G3*B)*(1-G4*B)*(1-G5*B)=(1-0.3*B)*(1-0.7*B)^4=(1-3.1*B+3.78*B^2-2.254*B^3+0.6517*B^4-0.07203*B^5)
    # > 1/polyroot(c(1,-3.1,3.78,-2.254,0.6517,-0.07203))
    # [1] 0.3+0i 0.7-0i 0.7-0i 0.7+0i 0.7+0i
    
    a=simulation(n_series, length_series, ma_coef=c(3.1,-3.78,2.254,-0.6517,0.07203),d=d,D=D)
    MA5_0.3_0.7_0.7_0.7_0.7_SLBDD <- mean(a[1,]==paste("MA(5)","_(",d,",",D,")",sep=""))
    MA5_0.3_0.7_0.7_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    MA5_0.3_0.7_0.7_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    MA5_0.3_0.7_0.7_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    MA5_0.3_0.7_0.7_0.7_0.7_forecast <- mean(a[5,]==paste("MA(5)","_(",d,",",D,")",sep=""))
    MA5_0.3_0.7_0.7_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    MA5_0.3_0.7_0.7_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    MA5_0.3_0.7_0.7_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # Joint MA(5) simulation and graphics
    
    complete_results[10,n_length,n_dif,"Acc","SLBDD"]=mean(c(MA5_0.7_0.7_0.7_0.7_0.7_SLBDD,MA5_0.3_0.7_0.7_0.7_0.7_SLBDD))
    complete_results[10,n_length,n_dif,"MSE","SLBDD"]=mean(c(MA5_0.7_0.7_0.7_0.7_0.7_SLBDD_mse,MA5_0.3_0.7_0.7_0.7_0.7_SLBDD_mse))
    complete_results[10,n_length,n_dif,"s2","SLBDD"]=mean(c(MA5_0.7_0.7_0.7_0.7_0.7_SLBDD_s2,MA5_0.3_0.7_0.7_0.7_0.7_SLBDD_s2))
    complete_results[10,n_length,n_dif,"MAE","SLBDD"]=mean(c(MA5_0.7_0.7_0.7_0.7_0.7_SLBDD_mae,MA5_0.3_0.7_0.7_0.7_0.7_SLBDD_mae))
    
    complete_results[10,n_length,n_dif,"Acc","forecast"]=mean(c(MA5_0.7_0.7_0.7_0.7_0.7_forecast,MA5_0.3_0.7_0.7_0.7_0.7_forecast))
    complete_results[10,n_length,n_dif,"MSE","forecast"]=mean(c(MA5_0.7_0.7_0.7_0.7_0.7_forecast_mse,MA5_0.3_0.7_0.7_0.7_0.7_forecast_mse))
    complete_results[10,n_length,n_dif,"s2","forecast"]=mean(c(MA5_0.7_0.7_0.7_0.7_0.7_forecast_s2,MA5_0.3_0.7_0.7_0.7_0.7_forecast_s2))
    complete_results[10,n_length,n_dif,"MAE","forecast"]=mean(c(MA5_0.7_0.7_0.7_0.7_0.7_forecast_mae,MA5_0.3_0.7_0.7_0.7_0.7_forecast_mae))
    
    # ARMA models
    
    # ARMA(1,1) with G1={0.3} G2={0.7}
    
    a=simulation(n_series, length_series, ar_coef=c(0.3),ma_coef=c(0.7),d=d,D=D)
    ARMA11_0.3_0.7_SLBDD <- mean(a[1,]==paste("ARMA(1,1)","_(",d,",",D,")",sep=""))
    ARMA11_0.3_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    ARMA11_0.3_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    ARMA11_0.3_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    ARMA11_0.3_0.7_forecast <- mean(a[5,]==paste("ARMA(1,1)","_(",d,",",D,")",sep=""))
    ARMA11_0.3_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    ARMA11_0.3_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    ARMA11_0.3_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # ARMA(1,1) with G1={0.7} G2={0.3}
    
    a=simulation(n_series, length_series, ar_coef=c(0.7),ma_coef=c(0.3),d=d,D=D)
    ARMA11_0.7_0.3_SLBDD <- mean(a[1,]==paste("ARMA(1,1)","_(",d,",",D,")",sep=""))
    ARMA11_0.7_0.3_SLBDD_mse <- mean(as.numeric(a[2,]))
    ARMA11_0.7_0.3_SLBDD_s2 <- mean(as.numeric(a[3,]))
    ARMA11_0.7_0.3_SLBDD_mae <- mean(as.numeric(a[4,]))
    ARMA11_0.7_0.3_forecast <- mean(a[5,]==paste("ARMA(1,1)","_(",d,",",D,")",sep=""))
    ARMA11_0.7_0.3_forecast_mse <- mean(as.numeric(a[6,]))
    ARMA11_0.7_0.3_forecast_s2 <- mean(as.numeric(a[7,]))
    ARMA11_0.7_0.3_forecast_mae <- mean(as.numeric(a[8,]))
    
    # Joint ARMA(1,1) simulation and graphics
    
    complete_results[11,n_length,n_dif,"Acc","SLBDD"]=mean(c(ARMA11_0.3_0.7_SLBDD,ARMA11_0.7_0.3_SLBDD))
    complete_results[11,n_length,n_dif,"MSE","SLBDD"]=mean(c(ARMA11_0.3_0.7_SLBDD_mse,ARMA11_0.7_0.3_SLBDD_mse))
    complete_results[11,n_length,n_dif,"s2","SLBDD"]=mean(c(ARMA11_0.3_0.7_SLBDD_s2,ARMA11_0.7_0.3_SLBDD_s2))
    complete_results[11,n_length,n_dif,"MAE","SLBDD"]=mean(c(ARMA11_0.3_0.7_SLBDD_mae,ARMA11_0.7_0.3_SLBDD_mae))
    
    complete_results[11,n_length,n_dif,"Acc","forecast"]=mean(c(ARMA11_0.3_0.7_forecast,ARMA11_0.7_0.3_forecast))
    complete_results[11,n_length,n_dif,"MSE","forecast"]=mean(c(ARMA11_0.3_0.7_forecast_mse,ARMA11_0.7_0.3_forecast_mse))
    complete_results[11,n_length,n_dif,"s2","forecast"]=mean(c(ARMA11_0.3_0.7_forecast_s2,ARMA11_0.7_0.3_forecast_s2))
    complete_results[11,n_length,n_dif,"MAE","forecast"]=mean(c(ARMA11_0.3_0.7_forecast_mae,ARMA11_0.7_0.3_forecast_mae))
    
    # ARMA(2,1) with G1={0.3, 0.3} G2={0.7}
    
    # z_t*(1-phi_1*z_(t-1)-...-phi_p*z_(t-p))=(1-theta_1*z_(t-1)-...-theta_q*z_(t-q))*a_t 
    # => z_t*(1-G1*B)*(1-G2*B)=(1-G3*B)*a_t 
    # In this case (1-G1*B)*(1-G2*B)=(1-0.3*B)^2=(1-0.6*B+0.09*B^2)
    # > 1/polyroot(c(1,-0.6,0.09))
    # [1] 0.3+0i 0.3-0i
    
    a=simulation(n_series, length_series, ar_coef=c(0.6, -0.09),ma_coef=c(0.7),d=d,D=D)
    ARMA21_0.3_0.3_0.7_SLBDD <- mean(a[1,]==paste("ARMA(2,1)","_(",d,",",D,")",sep=""))
    ARMA21_0.3_0.3_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    ARMA21_0.3_0.3_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    ARMA21_0.3_0.3_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    ARMA21_0.3_0.3_0.7_forecast <- mean(a[5,]==paste("ARMA(2,1)","_(",d,",",D,")",sep=""))
    ARMA21_0.3_0.3_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    ARMA21_0.3_0.3_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    ARMA21_0.3_0.3_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # ARMA(2,1) with G1={0.7, 0.7} G2={0.3}
    
    # In this case (1-G1*B)*(1-G2*B)=(1-0.7*B)^2=(1-1.4*B+0.49*B^2)
    # > 1/polyroot(c(1,-1.4,0.49))
    # [1] 0.7-0i 0.7+0i
    
    a=simulation(n_series, length_series, ar_coef=c(1.4, -0.49),ma_coef=c(0.3),d=d,D=D)
    ARMA21_0.7_0.7_0.3_SLBDD <- mean(a[1,]==paste("ARMA(2,1)","_(",d,",",D,")",sep=""))
    ARMA21_0.7_0.7_0.3_SLBDD_mse <- mean(as.numeric(a[2,]))
    ARMA21_0.7_0.7_0.3_SLBDD_s2 <- mean(as.numeric(a[3,]))
    ARMA21_0.7_0.7_0.3_SLBDD_mae <- mean(as.numeric(a[4,]))
    ARMA21_0.7_0.7_0.3_forecast <- mean(a[5,]==paste("ARMA(2,1)","_(",d,",",D,")",sep=""))
    ARMA21_0.7_0.7_0.3_forecast_mse <- mean(as.numeric(a[6,]))
    ARMA21_0.7_0.7_0.3_forecast_s2 <- mean(as.numeric(a[7,]))
    ARMA21_0.7_0.7_0.3_forecast_mae <- mean(as.numeric(a[8,]))
    
    # Joint ARMA(2,1) simulation and graphics
  
    complete_results[12,n_length,n_dif,"Acc","SLBDD"]=mean(c(ARMA21_0.3_0.3_0.7_SLBDD,ARMA21_0.7_0.7_0.3_SLBDD))
    complete_results[12,n_length,n_dif,"MSE","SLBDD"]=mean(c(ARMA21_0.3_0.3_0.7_SLBDD_mse,ARMA21_0.7_0.7_0.3_SLBDD_mse))
    complete_results[12,n_length,n_dif,"s2","SLBDD"]=mean(c(ARMA21_0.3_0.3_0.7_SLBDD_s2,ARMA21_0.7_0.7_0.3_SLBDD_s2))
    complete_results[12,n_length,n_dif,"MAE","SLBDD"]=mean(c(ARMA21_0.3_0.3_0.7_SLBDD_mae,ARMA21_0.7_0.7_0.3_SLBDD_mae))
    
    complete_results[12,n_length,n_dif,"Acc","forecast"]=mean(c(ARMA21_0.3_0.3_0.7_forecast,ARMA21_0.7_0.7_0.3_forecast))
    complete_results[12,n_length,n_dif,"MSE","forecast"]=mean(c(ARMA21_0.3_0.3_0.7_forecast_mse,ARMA21_0.7_0.7_0.3_forecast_mse))
    complete_results[12,n_length,n_dif,"s2","forecast"]=mean(c(ARMA21_0.3_0.3_0.7_forecast_s2,ARMA21_0.7_0.7_0.3_forecast_s2))
    complete_results[12,n_length,n_dif,"MAE","forecast"]=mean(c(ARMA21_0.3_0.3_0.7_forecast_mae,ARMA21_0.7_0.7_0.3_forecast_mae))
    
    # ARMA(1,2) with G1={0.3} G2={0.7, 0.7}
    
    # z_t*(1-phi_1*z_(t-1)-...-phi_p*z_(t-p))=(1-theta_1*z_(t-1)-...-theta_q*z_(t-q))*a_t 
    # => z_t*(1-G1*B)=(1-G2*B)*(1-G3*B)*a_t 
    # In this case (1-G2*B)*(1-G3*B)=(1-0.7*B)^2=(1-1.4*B+0.49*B^2)
    # > 1/polyroot(c(1,-1.4,0.49))
    # [1] 0.7-0i 0.7+0i
    
    a=simulation(n_series, length_series, ar_coef=c(0.3),ma_coef=c(1.4, -0.49),d=d,D=D)
    ARMA12_0.3_0.7_0.7_SLBDD <- mean(a[1,]==paste("ARMA(1,2)","_(",d,",",D,")",sep=""))
    ARMA12_0.3_0.7_0.7_SLBDD_mse <- mean(as.numeric(a[2,]))
    ARMA12_0.3_0.7_0.7_SLBDD_s2 <- mean(as.numeric(a[3,]))
    ARMA12_0.3_0.7_0.7_SLBDD_mae <- mean(as.numeric(a[4,]))
    ARMA12_0.3_0.7_0.7_forecast <- mean(a[5,]==paste("ARMA(1,2)","_(",d,",",D,")",sep=""))
    ARMA12_0.3_0.7_0.7_forecast_mse <- mean(as.numeric(a[6,]))
    ARMA12_0.3_0.7_0.7_forecast_s2 <- mean(as.numeric(a[7,]))
    ARMA12_0.3_0.7_0.7_forecast_mae <- mean(as.numeric(a[8,]))
    
    # ARMA(1,2) with G1={0.7} G2={0.3, 0.3}
    
    # In this case (1-G2*B)*(1-G3*B)=(1-0.3*B)^2=(1-0.6*B+0.09*B^2)
    # > 1/polyroot(c(1,-0.6,0.09))
    # [1] 0.3+0i 0.3-0i
    
    a=simulation(n_series, length_series, ar_coef=c(0.7),ma_coef=c(0.6, -0.09),d=d,D=D)
    ARMA12_0.7_0.3_0.3_SLBDD <- mean(a[1,]==paste("ARMA(1,2)","_(",d,",",D,")",sep=""))
    ARMA12_0.7_0.3_0.3_SLBDD_mse <- mean(as.numeric(a[2,]))
    ARMA12_0.7_0.3_0.3_SLBDD_s2 <- mean(as.numeric(a[3,]))
    ARMA12_0.7_0.3_0.3_SLBDD_mae <- mean(as.numeric(a[4,]))
    ARMA12_0.7_0.3_0.3_forecast <- mean(a[5,]==paste("ARMA(1,2)","_(",d,",",D,")",sep=""))
    ARMA12_0.7_0.3_0.3_forecast_mse <- mean(as.numeric(a[6,]))
    ARMA12_0.7_0.3_0.3_forecast_s2 <- mean(as.numeric(a[7,]))
    ARMA12_0.7_0.3_0.3_forecast_mae <- mean(as.numeric(a[8,]))
    
    # Joint ARMA(1,2) simulation and graphics
    
    complete_results[13,n_length,n_dif,"Acc","SLBDD"]=mean(c(ARMA12_0.3_0.7_0.7_SLBDD,ARMA12_0.7_0.3_0.3_SLBDD))
    complete_results[13,n_length,n_dif,"MSE","SLBDD"]=mean(c(ARMA12_0.3_0.7_0.7_SLBDD_mse,ARMA12_0.7_0.3_0.3_SLBDD_mse))
    complete_results[13,n_length,n_dif,"s2","SLBDD"]=mean(c(ARMA12_0.3_0.7_0.7_SLBDD_s2,ARMA12_0.7_0.3_0.3_SLBDD_s2))
    complete_results[13,n_length,n_dif,"MAE","SLBDD"]=mean(c(ARMA12_0.3_0.7_0.7_SLBDD_mae,ARMA12_0.7_0.3_0.3_SLBDD_mae))
  
    complete_results[13,n_length,n_dif,"Acc","forecast"]=mean(c(ARMA12_0.3_0.7_0.7_forecast,ARMA12_0.7_0.3_0.3_forecast))
    complete_results[13,n_length,n_dif,"MSE","forecast"]=mean(c(ARMA12_0.3_0.7_0.7_forecast_mse,ARMA12_0.7_0.3_0.3_forecast_mse))
    complete_results[13,n_length,n_dif,"s2","forecast"]=mean(c(ARMA12_0.3_0.7_0.7_forecast_s2,ARMA12_0.7_0.3_0.3_forecast_s2))
    complete_results[13,n_length,n_dif,"MAE","forecast"]=mean(c(ARMA12_0.3_0.7_0.7_forecast_mae,ARMA12_0.7_0.3_0.3_forecast_mae))

  }
}

# Time

Time=toc()

# Save results

save(complete_results,file = "Simulation_results.RData")

# Load results

 load("Resultados.RData")

# Final graphics

datos_acc <- data.frame(
  Conjunto_de_datos = c("AR(1)", "AR(2)", "AR(3)","AR(4)","AR(5)","MA(1)", "MA(2)", "MA(3)","MA(4)","MA(5)",
                        "ARMA(1,1)", "ARMA(2,1)", "ARMA(1,2)"),
  SLBDD = apply(complete_results[,,,"Acc","SLBDD"],1,mean),
  forecast = apply(complete_results[,,,"Acc","forecast"],1,mean) 
)

datos_largos_acc <- gather(datos_acc, key = "Modelo", value = "Accuracy", -Conjunto_de_datos)
datos_largos_acc$Conjunto_de_datos <- factor(datos_largos_acc$Conjunto_de_datos, 
                                         levels = c("AR(1)", "AR(2)", "AR(3)","AR(4)","AR(5)","MA(1)", "MA(2)", "MA(3)","MA(4)","MA(5)",
                                                    "ARMA(1,1)", "ARMA(2,1)", "ARMA(1,2)"))

pAcc=ggplot(datos_largos_acc, aes(x = Conjunto_de_datos, y = Accuracy, fill = Modelo)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylim(0, 1) +
  labs(title = "SLBDD vs forecast",
       x = "",
       y = "Mean Accuracy",
       fill = "Library") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))
pAcc

datos_mse <- data.frame(
  Conjunto_de_datos = c("AR(1)", "AR(2)", "AR(3)","AR(4)","AR(5)","MA(1)", "MA(2)", "MA(3)","MA(4)","MA(5)",
                        "ARMA(1,1)", "ARMA(2,1)", "ARMA(1,2)"),
  SLBDD = apply(complete_results[,,,"MSE","SLBDD"],1,mean),
  forecast = apply(complete_results[,,,"MSE","forecast"],1,mean) 
)

datos_largos_mse <- gather(datos_mse, key = "Modelo", value = "Mean_MSE", -Conjunto_de_datos)
datos_largos_mse$Conjunto_de_datos <- factor(datos_largos_mse$Conjunto_de_datos, 
                                             levels = c("AR(1)", "AR(2)", "AR(3)","AR(4)","AR(5)","MA(1)", "MA(2)", "MA(3)","MA(4)","MA(5)",
                                                        "ARMA(1,1)", "ARMA(2,1)", "ARMA(1,2)"))

pMSE=ggplot(datos_largos_mse, aes(x = Conjunto_de_datos, y = Mean_MSE, fill = Modelo)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylim(0, 1) +
  labs(title = "SLBDD vs forecast",
       x = "",
       y = "Mean MSE",
       fill = "Library") +
  theme_minimal()+ ylim(c(0,max(datos_largos_mse$Mean_MSE)))+
  theme(plot.title = element_text(hjust = 0.5))
pMSE

datos_s2 <- data.frame(
  Conjunto_de_datos = c("AR(1)", "AR(2)", "AR(3)","AR(4)","AR(5)","MA(1)", "MA(2)", "MA(3)","MA(4)","MA(5)",
                        "ARMA(1,1)", "ARMA(2,1)", "ARMA(1,2)"),
  SLBDD = apply(complete_results[,,,"s2","SLBDD"],1,mean),
  forecast = apply(complete_results[,,,"s2","forecast"],1,mean) 
)

datos_largos_s2 <- gather(datos_s2, key = "Modelo", value = "s2", -Conjunto_de_datos)
datos_largos_s2$Conjunto_de_datos <- factor(datos_largos_s2$Conjunto_de_datos, 
                                             levels = c("AR(1)", "AR(2)", "AR(3)","AR(4)","AR(5)","MA(1)", "MA(2)", "MA(3)","MA(4)","MA(5)",
                                                        "ARMA(1,1)", "ARMA(2,1)", "ARMA(1,2)"))

ps2=ggplot(datos_largos_s2, aes(x = Conjunto_de_datos, y = s2, fill = Modelo)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylim(0, 1) +
  labs(title = "SLBDD vs forecast",
       x = "",
       y = expression(paste("Mean ",sigma^2)),
       fill = "Library") +ylim(c(0,max(datos_largos_s2$s2)))+
  theme_minimal()+ 
theme(plot.title = element_text(hjust = 0.5))
ps2

datos_MAE <- data.frame(
  Conjunto_de_datos = c("AR(1)", "AR(2)", "AR(3)","AR(4)","AR(5)","MA(1)", "MA(2)", "MA(3)","MA(4)","MA(5)",
                        "ARMA(1,1)", "ARMA(2,1)", "ARMA(1,2)"),
  SLBDD = apply(complete_results[,,,"MAE","SLBDD"],1,mean),
  forecast = apply(complete_results[,,,"MAE","forecast"],1,mean) 
)

datos_largos_MAE <- gather(datos_MAE, key = "Modelo", value = "MAE", -Conjunto_de_datos)
datos_largos_MAE$Conjunto_de_datos <- factor(datos_largos_MAE$Conjunto_de_datos, 
                                            levels = c("AR(1)", "AR(2)", "AR(3)","AR(4)","AR(5)","MA(1)", "MA(2)", "MA(3)","MA(4)","MA(5)",
                                                       "ARMA(1,1)", "ARMA(2,1)", "ARMA(1,2)"))

pMAE=ggplot(datos_largos_MAE, aes(x = Conjunto_de_datos, y = MAE, fill = Modelo)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylim(0, 1) +
  labs(title = "SLBDD vs forecast",
       x = "",
       y = "Mean MAE",
       fill = "Library") +ylim(c(0,max(datos_largos_MAE$MAE)))+
  theme_minimal()+ 
  theme(plot.title = element_text(hjust = 0.5))
pMAE

# Mean final results

## Accuracy

### n

SLBDD_acc_n=apply(complete_results[,,,"Acc","SLBDD"],2,mean)
SLBDD_acc_n

forecast_acc_n=apply(complete_results[,,,"Acc","forecast"],2,mean)
forecast_acc_n

plot(SLBDD_acc_n,type="l", ylab="Mean Accuracy", xlab="Length of the series",col="#00BFC4",lwd=2, xaxt="n",ylim=c(min(c(SLBDD_acc_n,forecast_acc_n)),max(c(SLBDD_acc_n,forecast_acc_n))))
axis(side=1,at=c(1,2,3),labels=length_series_vec)
lines(forecast_acc_n, type="l",col="#F8766D",lwd=2)
legend("topleft",legend=c("SLBDD","forecast"),fill=c("#00BFC4","#F8766D"))

### dD

SLBDD_acc_dD=apply(complete_results[,,,"Acc","SLBDD"],3,mean)
SLBDD_acc_dD

forecast_acc_dD=apply(complete_results[,,,"Acc","forecast"],3,mean)
forecast_acc_dD


plot(SLBDD_acc_dD,type="l", ylab="Mean Accuracy", xlab="(d, D)",col="#00BFC4",lwd=2, xaxt="n",ylim=c(min(c(SLBDD_acc_dD,forecast_acc_dD)),max(c(SLBDD_acc_dD,forecast_acc_dD))))
axis(side=1,at=c(1,2,3),labels=substring(pairs_dD_differences_vec,2))
lines(forecast_acc_dD, type="l",col="#F8766D",lwd=2)
legend("topleft",legend=c("SLBDD","forecast"),fill=c("#00BFC4","#F8766D"))

### Mean Accuracy

mean(complete_results[,,,"Acc","SLBDD"])
mean(complete_results[,,,"Acc","forecast"])

## MSE

### n

SLBDD_MSE_n=apply(complete_results[,,,"MSE","SLBDD"],2,mean)
SLBDD_MSE_n

forecast_MSE_n=apply(complete_results[,,,"MSE","forecast"],2,mean)
forecast_MSE_n

plot(SLBDD_MSE_n,type="l", ylab="Mean MSE", xlab="Length of the series",col="#00BFC4",lwd=2, xaxt="n",ylim=c(min(c(SLBDD_MSE_n,forecast_MSE_n)),max(c(SLBDD_MSE_n,forecast_MSE_n))))
axis(side=1,at=c(1,2,3),labels=length_series_vec)
lines(forecast_MSE_n, type="l",col="#F8766D",lwd=2)
legend("topright",legend=c("SLBDD","forecast"),fill=c("#00BFC4","#F8766D"))

### dD

SLBDD_MSE_dD=apply(complete_results[,,,"MSE","SLBDD"],3,mean)
SLBDD_MSE_dD

forecast_MSE_dD=apply(complete_results[,,,"MSE","forecast"],3,mean)
forecast_MSE_dD


plot(SLBDD_MSE_dD,type="l", ylab="Mean MSE", xlab="(d, D)",col="#00BFC4",lwd=2, xaxt="n",ylim=c(min(c(SLBDD_MSE_dD,forecast_MSE_dD)),max(c(SLBDD_MSE_dD,forecast_MSE_dD))))
axis(side=1,at=c(1,2,3),labels=substring(pairs_dD_differences_vec,2))
lines(forecast_MSE_dD, type="l",col="#F8766D",lwd=2)
legend("topleft",legend=c("SLBDD","forecast"),fill=c("#00BFC4","#F8766D"))

### Mean MSE

mean(complete_results[,,,"MSE","SLBDD"])
mean(complete_results[,,,"MSE","forecast"])

## s2

### n

SLBDD_s2_n=apply(complete_results[,,,"s2","SLBDD"],2,mean)
SLBDD_s2_n

forecast_s2_n=apply(complete_results[,,,"s2","forecast"],2,mean)
forecast_s2_n

plot(SLBDD_s2_n,type="l", ylab=expression(paste("Mean ",sigma^2)), xlab="Length of the series",col="#00BFC4",lwd=2, xaxt="n",ylim=c(min(c(SLBDD_s2_n,forecast_s2_n)),max(c(SLBDD_s2_n,forecast_s2_n))))
axis(side=1,at=c(1,2,3),labels=length_series_vec)
lines(forecast_s2_n, type="l",col="#F8766D",lwd=2)
legend("topright",legend=c("SLBDD","forecast"),fill=c("#00BFC4","#F8766D"))

### dD

SLBDD_s2_dD=apply(complete_results[,,,"s2","SLBDD"],3,mean)
SLBDD_s2_dD

forecast_s2_dD=apply(complete_results[,,,"s2","forecast"],3,mean)
forecast_s2_dD


plot(SLBDD_s2_dD,type="l", ylab=expression(paste("Mean ",sigma^2)), xlab="(d, D)",col="#00BFC4",lwd=2, xaxt="n",ylim=c(min(c(SLBDD_s2_dD,forecast_s2_dD)),max(c(SLBDD_s2_dD,forecast_s2_dD))))
axis(side=1,at=c(1,2,3),labels=substring(pairs_dD_differences_vec,2))
lines(forecast_s2_dD, type="l",col="#F8766D",lwd=2)
legend("topleft",legend=c("SLBDD","forecast"),fill=c("#00BFC4","#F8766D"))

### Mean s2

mean(complete_results[,,,"s2","SLBDD"])
mean(complete_results[,,,"s2","forecast"])

## MAE

### n

SLBDD_MAE_n=apply(complete_results[,,,"MAE","SLBDD"],2,mean)
SLBDD_MAE_n

forecast_MAE_n=apply(complete_results[,,,"MAE","forecast"],2,mean)
forecast_MAE_n

plot(SLBDD_MAE_n,type="l", ylab="Mean MAE", xlab="Length of the series",col="#00BFC4",lwd=2, xaxt="n",ylim=c(min(c(SLBDD_MAE_n,forecast_MAE_n)),max(c(SLBDD_MAE_n,forecast_MAE_n))))
axis(side=1,at=c(1,2,3),labels=length_series_vec)
lines(forecast_MAE_n, type="l",col="#F8766D",lwd=2)
legend("topright",legend=c("SLBDD","forecast"),fill=c("#00BFC4","#F8766D"))

### dD

SLBDD_MAE_dD=apply(complete_results[,,,"MAE","SLBDD"],3,mean)
SLBDD_MAE_dD

forecast_MAE_dD=apply(complete_results[,,,"MAE","forecast"],3,mean)
forecast_MAE_dD


plot(SLBDD_MAE_dD,type="l", ylab="Mean MAE", xlab="(d, D)",col="#00BFC4",lwd=2, xaxt="n",ylim=c(min(c(SLBDD_MAE_dD,forecast_MAE_dD)),max(c(SLBDD_MAE_dD,forecast_MAE_dD))))
axis(side=1,at=c(1,2,3),labels=substring(pairs_dD_differences_vec,2))
lines(forecast_MAE_dD, type="l",col="#F8766D",lwd=2)
legend("topleft",legend=c("SLBDD","forecast"),fill=c("#00BFC4","#F8766D"))

### Mean MAE

mean(complete_results[,,,"MAE","SLBDD"])
mean(complete_results[,,,"MAE","forecast"])

