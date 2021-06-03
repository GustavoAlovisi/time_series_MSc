library(naniar)
library(quantmod)
library(TTR)
library(rugarch)
library(stochvol)
library(summarytools)

## getting stock and index data
getSymbols(c('PETR4.SA','^BVSP','ITUB4.SA','VALE3.SA','BBAS3.SA'), src = 'yahoo', auto.assign = TRUE)

## binding in one multivariate xts object for plotting and NA handling
PETR4.SA <- PETR4.SA$PETR4.SA.Close
BVSP <- BVSP$BVSP.Close
ITUB4.SA <- ITUB4.SA$ITUB4.SA.Close
VALE3.SA <- VALE3.SA$VALE3.SA.Close
BBAS3.SA <- BBAS3.SA$BBAS3.SA.Close
stockPrices <- cbind(PETR4.SA, BVSP, ITUB4.SA, VALE3.SA, BBAS3.SA)

## viz miss 
naniar::vis_miss(as.data.frame(stockPrices))

## na.approx
stockPrices <- zoo::na.approx(stockPrices)

## plotting price series
xts::plot.xts(stockPrices, multi.panel = TRUE,  yaxis.same=FALSE, grid.col = "lightgray", main = NA)

## calculating log-Returns
stocksRet <- diff(log(stockPrices))
stocksRet <- stocksRet[-1,]

## plotting return series
xts::plot.xts(stocksRet, multi.panel = TRUE, yaxis.same=FALSE, grid.col = "lightgray", main = NA)

## descriptive stats of returns
summarytools::descr(stocksRet, round.digits = 3)


descr()
#### In sample GARCH and SV fit and volatility plot for BVSP

### GARCH(1,1)-t 
garchSpec = rugarch::ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1), variance.targeting = TRUE), 
                                mean.model=list(armaOrder=c(0,0)),  
                                distribution.model='std') ##vamos especificar uma 't'  ao inves de sstd

#in sample
inSampleGarch_t <- rugarch::ugarchfit(data = stocksRet[,2], spec = garchSpec)


### GARCH(1,1)-sstd
garchSpec = rugarch::ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1), variance.targeting = TRUE), 
                                mean.model=list(armaOrder=c(0,0)),  
                                distribution.model='sstd') ##vamos especificar uma 't'  ao inves de sstd

#in sample
inSampleGarch_sstd <- rugarch::ugarchfit(data = stocksRet[,2], spec = garchSpec)


###fitando um AR-SV-t

#in sample
inSampleSV <- svtsample(stocksRet[,2], draws = 50000, burnin = 1000, parallel = 'snow', n_cpus = 8)


### comparing in sample volatility 
plotInsample<- cbind(sigma(inSampleGarch_t),sigma(inSampleGarch_sstd),inSampleSV$summary$sd[,4])
names(plotInsample) <- c('GARCH(1,1)-t', 'GARCH(1,1)-sstd','SV-AR(1)-t')
xts::plot.xts(plotInsample, main = 'IBOV in sample vol', legend.loc = 'topleft', lwd = 1, grid.col = "lightgray", col = c(2,4,7))


### GARCH/StochVol rolling estimation with out-of-sample forecasts 
for (j in 2:5){
  ret <- stocksRet[,j]
  #nameret <- names(ret)
  garchSpec = rugarch::ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1),variance.targeting = TRUE), 
                                  mean.model=list(armaOrder=c(0,0)),  
                                  distribution.model='std')
  
  outSampleGarch_t <- 0
  outSampleGarch_t <- ugarchroll(garchSpec, data = ret, n.ahead = 1, 
                               n.start = 2100,  refit.every = 1, refit.window = "moving", 
                               solver = "hybrid", fit.control = list(),
                               calculate.VaR = TRUE, VaR.alpha = c(0.05),
                               keep.coef = TRUE)
  
  GarchVaR_t <- outSampleGarch_t@forecast$VaR$`alpha(5%)`
  RealizedRet <- outSampleGarch_t@forecast$VaR$realized
  
  ## out of sample GARCH(1,1)-sstd
  garchSpecSkew = rugarch::ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1),variance.targeting = TRUE), 
                                      mean.model=list(armaOrder=c(0,0)),  
                                      distribution.model='sstd') 
  
  outSampleGarch_sstd <- 0
  outSampleGarch_sstd <- ugarchroll(garchSpecSkew, data = ret, n.ahead = 1, 
                                   n.start = 2100,  refit.every = 1, refit.window = "moving", 
                                   solver = "hybrid", fit.control = list(),
                                   calculate.VaR = TRUE, VaR.alpha = c(0.05),
                                   keep.coef = TRUE)
  
  GarchVaR_sstd <- outSampleGarch_sstd@forecast$VaR$`alpha(5%)`
  

  ## out of sample SV-AR(1)
  stochVol <- 0
  stochVol <- stochvol::svtsample_roll(ret, n_ahead = 1, refit_window = 'moving', 
                                       n_start = 2101, parallel = 'snow', n_cpus = 8, 
                                       draws = 10000, refit_every = 1, calculate_quantile = 0.05,
                                       keep_draws = FALSE, calculate_predictive_likelihood = FALSE) ##VaR = 0.05 
  
  # salvando vetor com a prediction do VaR
  SVVaR <- 0
  for (i in 1:1466){
    SVVaR[i] <- stochVol[[i]]$predicted_quantile
  }
 # SVVaR <- SVVaR[-1]
  estimations <- 0
  estimations <- cbind(GarchVaR_t,GarchVaR_sstd, SVVaR, RealizedRet)
  saveRDS(estimations, paste('Estimation',j,'.RDS', sep = ""))
}




### Usando GAS para prever VaR: no momento nao vamos utilizar

#library(GAS)

#GASSpec_std  <- GAS::UniGASSpec(Dist = "std",  GASPar = list(scale = TRUE))

#cluster <- makeCluster(8)
#outOfSampleGAS <- GAS::UniGASRoll(petrRet,GASSpec = GASSpec_std, Nstart = 1000, 
#                                  RefitEvery = 1, RefitWindow = 'moving', cluster = cluster)

#GASVaR <- quantile(outOfSampleGAS, 0.05)
#GASVaR <- GASVaR[,1]



### VaR model comparison 

# Loading estimation data
PETR4_VaR <- readRDS('estimation1.RDS')

BVSP_VaR <- readRDS('estimation2.RDS')

ITUB_VaR <- readRDS('estimation3.RDS')

VALE_VaR <- readRDS('estimation4.RDS')

BBAS_VaR <- readRDS('estimation5.RDS')


### 5% VaR Forecasting Plot for BVSP

xaxis <- index(stocksRet)
xaxis <- xaxis[2101:3567]
xaxis <- c(xaxis[1], xaxis[750], xaxis[1467])
plot(BVSP_VaR[,4], type = "l", col = grey(0.4, 0.5), 
     ylab = "IBOV Log-Returns", main = "IBOV VaR (5%) Forecasting", xaxt = "n")
axis(1, at=c(1, 750, 1467), labels=xaxis)
lines(BVSP_VaR[,1], col = 2)
lines(BVSP_VaR[,2], col = 4)
lines(BVSP_VaR[,3], col = 7)
#lines(GASVaR, col = 3)
abline(h=0, col = 1)
legend('topright', c("GARCH(1,1)-t", 'GARCH(1,1)-sstd', "SV-AR(1)"), 
       lty=1, col=c(2,4, 7), bty='n', cex=.75)


## PETR4 comparison

Garch_tPerf <- GAS::BacktestVaR(PETR4_VaR[,4],PETR4_VaR[,1], alpha = 0.05)
Garch_sstdPerf <- GAS::BacktestVaR(PETR4_VaR[,4],PETR4_VaR[,2], alpha = 0.05)
SVVARPerf <- GAS::BacktestVaR(PETR4_VaR[,4],PETR4_VaR[,3], alpha = 0.05)

# Comparing the the three

# AD: Absolute Deviation
Garch_tPerf$AD
Garch_sstdPerf$AD
SVVARPerf$AD

# DQ Correct 5% VaR specification: H0: Incorrect
Garch_tPerf$DQ$pvalue
Garch_sstdPerf$DQ$pvalue
SVVARPerf$DQ$pvalue

# AE: < 1: conservative >1: underestimate risk
Garch_tPerf$AE
Garch_sstdPerf$AE
SVVARPerf$AE

# QL: Quantile Loss - model comparison
Garch_tPerf$Loss$Loss
Garch_sstdPerf$Loss$Loss
SVVARPerf$Loss$Loss

# LRcc 
Garch_tPerf$LRcc
Garch_sstdPerf$LRcc
SVVARPerf$LRcc

# LRuc
Garch_tPerf$LRuc
Garch_sstdPerf$LRuc
SVVARPerf$LRuc




## BVSP comparison

Garch_tPerf <- GAS::BacktestVaR(BVSP_VaR[,4],BVSP_VaR[,1], alpha = 0.05)
Garch_sstdPerf <- GAS::BacktestVaR(BVSP_VaR[,4],BVSP_VaR[,2], alpha = 0.05)
SVVARPerf <- GAS::BacktestVaR(BVSP_VaR[,4],BVSP_VaR[,3], alpha = 0.05)

# Comparing the the three

# AD: Absolute Deviation
Garch_tPerf$AD
Garch_sstdPerf$AD
SVVARPerf$AD

# DQ Correct 5% VaR specification: H0: Incorrect
Garch_tPerf$DQ$pvalue
Garch_sstdPerf$DQ$pvalue
SVVARPerf$DQ$pvalue

# AE: < 1: conservative >1: underestimate risk
Garch_tPerf$AE
Garch_sstdPerf$AE
SVVARPerf$AE

# QL: Quantile Loss - model comparison
Garch_tPerf$Loss$Loss
Garch_sstdPerf$Loss$Loss
SVVARPerf$Loss$Loss


# LRcc 
Garch_tPerf$LRcc
Garch_sstdPerf$LRcc
SVVARPerf$LRcc

# LRuc
Garch_tPerf$LRuc
Garch_sstdPerf$LRuc
SVVARPerf$LRuc





## ITUB comparison

Garch_tPerf <- GAS::BacktestVaR(ITUB_VaR[,4],ITUB_VaR[,1], alpha = 0.05)
Garch_sstdPerf <- GAS::BacktestVaR(ITUB_VaR[,4],ITUB_VaR[,2], alpha = 0.05)
SVVARPerf <- GAS::BacktestVaR(ITUB_VaR[,4],ITUB_VaR[,3], alpha = 0.05)

# Comparing the the three

# AD: Absolute Deviation
Garch_tPerf$AD
Garch_sstdPerf$AD
SVVARPerf$AD

# DQ Correct 5% VaR specification: H0: Incorrect
Garch_tPerf$DQ$pvalue
Garch_sstdPerf$DQ$pvalue
SVVARPerf$DQ$pvalue

# AE: < 1: conservative >1: underestimate risk
Garch_tPerf$AE
Garch_sstdPerf$AE
SVVARPerf$AE

# QL: Quantile Loss - model comparison
Garch_tPerf$Loss$Loss
Garch_sstdPerf$Loss$Loss
SVVARPerf$Loss$Loss


# LRcc 
Garch_tPerf$LRcc
Garch_sstdPerf$LRcc
SVVARPerf$LRcc

# LRuc
Garch_tPerf$LRuc
Garch_sstdPerf$LRuc
SVVARPerf$LRuc



## VALE comparison

Garch_tPerf <- GAS::BacktestVaR(VALE_VaR[,4],VALE_VaR[,1], alpha = 0.05)
Garch_sstdPerf <- GAS::BacktestVaR(VALE_VaR[,4],VALE_VaR[,2], alpha = 0.05)
SVVARPerf <- GAS::BacktestVaR(VALE_VaR[,4],VALE_VaR[,3], alpha = 0.05)

# Comparing the the three

# AD: Absolute Deviation
Garch_tPerf$AD
Garch_sstdPerf$AD
SVVARPerf$AD

# DQ Correct 5% VaR specification: H0: Incorrect
Garch_tPerf$DQ$pvalue
Garch_sstdPerf$DQ$pvalue
SVVARPerf$DQ$pvalue

# AE: < 1: conservative >1: underestimate risk
Garch_tPerf$AE
Garch_sstdPerf$AE
SVVARPerf$AE

# QL: Quantile Loss - model comparison
Garch_tPerf$Loss$Loss
Garch_sstdPerf$Loss$Loss
SVVARPerf$Loss$Loss


# LRcc 
Garch_tPerf$LRcc
Garch_sstdPerf$LRcc
SVVARPerf$LRcc

# LRuc
Garch_tPerf$LRuc
Garch_sstdPerf$LRuc
SVVARPerf$LRuc




## BBAS comparison

Garch_tPerf <- GAS::BacktestVaR(BBAS_VaR[,4],BBAS_VaR[,1], alpha = 0.05)
Garch_sstdPerf <- GAS::BacktestVaR(BBAS_VaR[,4],BBAS_VaR[,2], alpha = 0.05)
SVVARPerf <- GAS::BacktestVaR(BBAS_VaR[,4],BBAS_VaR[,3], alpha = 0.05)

# Comparing the the three

# AD: Absolute Deviation
Garch_tPerf$AD
Garch_sstdPerf$AD
SVVARPerf$AD

# DQ Correct 5% VaR specification: H0: Incorrect
Garch_tPerf$DQ$pvalue
Garch_sstdPerf$DQ$pvalue
SVVARPerf$DQ$pvalue

# AE: < 1: conservative >1: underestimate risk
Garch_tPerf$AE
Garch_sstdPerf$AE
SVVARPerf$AE

# QL: Quantile Loss - model comparison
Garch_tPerf$Loss$Loss
Garch_sstdPerf$Loss$Loss
SVVARPerf$Loss$Loss


# LRcc 
Garch_tPerf$LRcc
Garch_sstdPerf$LRcc
SVVARPerf$LRcc

# LRuc
Garch_tPerf$LRuc
Garch_sstdPerf$LRuc
SVVARPerf$LRuc







a <- read.csv('efetivos_val_aux.csv')