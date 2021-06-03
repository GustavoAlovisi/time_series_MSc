

library(naniar)
library(quantmod)
library(TTR)
library(rugarch)
library(stochvol)



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

### na.approx
stockPrices <- zoo::na.approx(stockPrices)

###calculando log-Returns
petrRet <- diff(log(petrPrice))
petrRet <- petrRet[-1,]

###plotando a série
xts::plot.xts(petrRet)

###fitando um garch(1,1): 
garchSpec = rugarch::ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), 
                  mean.model=list(armaOrder=c(0,0)),  
                  distribution.model='std') ##vamos especificar uma 't'  ao inves de sstd

#in sample
inSampleGarch <- rugarch::ugarchfit(data = petrRet, spec = garchSpec)
plot(sigma(inSampleGarch))

fitted(inSampleGarch)
mean(petrRet)

###fitando um AR-SV-t

#in sample
inSampleSV <- svtsample(petrRet, draws = 50000, burnin = 1000, parallel = 'snow', n_cpus = 8)

plot(inSampleSV)
summary(inSampleSV)
volplot(inSampleSV)

inSampleSV$summary$sd

###comparando in Sample garch-t e AR-SV-t
plotInsample<- cbind(sigma(inSampleGarch),inSampleSV$summary$sd[,4])
plot(plotInsample)





#######estimação out of sample do VaR: 1000 observações para treino, 1 step ahead forecast do VaR

## out of sample GARCH(1,1)-t
outSampleGarch <- ugarchroll(garchSpec, data = petrRet, n.ahead = 1, 
                             n.start = 1000,  refit.every = 1, refit.window = "moving", 
                             solver = "hybrid", fit.control = list(),
                             calculate.VaR = TRUE, VaR.alpha = c(0.05),
                             keep.coef = TRUE)

GarchVaR <- outSampleGarch@forecast$VaR$`alpha(5%)`
RealizedRet <- outSampleGarch@forecast$VaR$realized
#resume(outSampleGarch)

#plot(outSampleGarch@forecast$density$Sigma)

#plot(outSampleGarch@forecast$VaR$`alpha(5%)`)

saveRDS(outSampleGarch, 'garchEst.RDS')


## out of sample GARCH(1,1)-sstd
garchSpecSkew = rugarch::ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), 
                                mean.model=list(armaOrder=c(0,0)),  
                                distribution.model='sstd') 

outSampleGarchSkew <- ugarchroll(garchSpecSkew, data = petrRet, n.ahead = 1, 
                             n.start = 1000,  refit.every = 1, refit.window = "moving", 
                             solver = "hybrid", fit.control = list(),
                             calculate.VaR = TRUE, VaR.alpha = c(0.05),
                             keep.coef = TRUE)

GarchVaRSkew <- outSampleGarchSkew@forecast$VaR$`alpha(5%)`
saveRDS(GarchVaRSkew, 'garchEstSkew.RDS')

## out of sample SV-AR(1)
stochVol <- stochvol::svtsample_roll(petrRet, n_ahead = 1, refit_window = 'moving', 
                                     n_start = 1000, parallel = 'snow', n_cpus = 8, 
                                     draws = 10000, refit_every = 1, calculate_quantile = 0.05,
                                     keep_draws = FALSE, calculate_predictive_likelihood = FALSE) ##VaR = 0.05 
 
saveRDS(stochVol, 'SVolEst.RDS')
#VaR = fitted(filt) + sigma(filt)*qdist("sstd", p=0.05, mu = 0, sigma = 1,skew  = coef(fit)["skew"], shape=coef(fit)["shape"])

# salvando vetor com a prediction do VaR
SVVaR <- 0
for (i in 1:2565){
  SVVaR[i] <- stochVol[[i]]$predicted_quantile
}

SVVaR <- SVVaR[-1]


### Usando GAS para prever VaR: no momento nao vamos utilizar

#library(GAS)

#GASSpec_std  <- GAS::UniGASSpec(Dist = "std",  GASPar = list(scale = TRUE))

#cluster <- makeCluster(8)
#outOfSampleGAS <- GAS::UniGASRoll(petrRet,GASSpec = GASSpec_std, Nstart = 1000, 
#                                  RefitEvery = 1, RefitWindow = 'moving', cluster = cluster)

#GASVaR <- quantile(outOfSampleGAS, 0.05)
#GASVaR <- GASVaR[,1]




### 5% VaR Forecasting Plot

xaxis <- index(petrRet)
xaxis <- xaxis[1001:3564]
xaxis <- c(xaxis[1], xaxis[501], xaxis[1001], xaxis[1501],xaxis[2001],xaxis[2501])
plot((RealizedRet), type = "l", col = grey(0.4, 0.5), 
     ylab = "PETR4 Log Returns", main = "VaR (5%) Forecasting", xaxt = "n")
axis(1, at=c(1, 500, 1000, 1500, 2000, 2501), labels=xaxis)
lines(SVVaR, col = 2)
lines(GarchVaR, col = 4)
lines(GarchVaRSkew, col = 7)
#lines(GASVaR, col = 3)
abline(h=0, col = 1)
legend('topright', c("GARCH(1,1)-t",  "SV-AR(1)", 'GARCH(1,1)-sstd') , 
       lty=1, col=c(2,4, 7), bty='n', cex=.75)


### VaR model comparison 

OOSlnRet <- petrRet[1001:3564,]

# GARCH(1,1)
GarchVaRPerf <- GAS::BacktestVaR(OOSlnRet, GarchVaR, alpha = 0.05)

# SVAR(1)
SVVaRPerf <- GAS::BacktestVaR(OOSlnRet, SVVaR, alpha = 0.05)

## Comparing the two 

# AD: Absolute Deviation
GarchVaRPerf$AD
SVVaRPerf$AD

# DQ Correct 5% VaR specification: H0: Incorrect
GarchVaRPerf$DQ$pvalue
SVVaRPerf$DQ$pvalue

# AE: < 1: conservative >1: underestimate risk
GarchVaRPerf$AE
SVVaRPerf$AE

# QL: Quantile Loss
round(GarchVaRPerf$Loss$Loss/SVVaRPerf$Loss$Loss,2)
# GARCH outperforms SVAR by 1% wrt Quantile Loss

# LRcc 
GarchVaRPerf$LRcc
SVVaRPerf$LRcc

# LRuc
GarchVaRPerf$LRuc
SVVaRPerf$LRuc


