rm(list=ls())
library(quantmod)
library(moments)

# Retrieve 25 stock symbols for historical analysis since the Great Recession of 2007-2009
# These 25 stocks were chosen because they all more or less bounced off their 200 week moving average, 
#or never touched such average, which indicates to me a technically sound and attractive stock chart
# Also, the 50 week moving average stayed above 200 week moving average for most of its life, suggesting good technicals.

getSymbols('AAPL', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('MSFT', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('GOOGL', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('AMZN', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('BRK-B', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('NVDA', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('V', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('MA', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('WMT', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('HD', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('TMO', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('COST', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('ABT', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('ADBE', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('DHR', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('NKE', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('INTC', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('AMD', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('TMUS', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('AMGN', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('LMT', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('INTU', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('PGR', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('ICE', from='2009-6-1', to='2022-6-1', src='yahoo')
getSymbols('DG', from='2009-6-1', to='2022-6-1', src='yahoo')


# Calculate the yearly log returns of each of the 25 stock symbols
AAPL.Logyearly = periodReturn(AAPL$AAPL.Close, period='yearly', type='log')
MSFT.Logyearly = periodReturn(MSFT$MSFT.Close, period='yearly', type='log')
GOOGL.Logyearly = periodReturn(GOOGL$GOOGL.Close, period='yearly', type='log')
AMZN.Logyearly = periodReturn(AMZN$AMZN.Close, period='yearly', type='log')
BRK.B.Logyearly = periodReturn(`BRK-B`$`BRK-B.Close`, period='yearly', type='log')
NVDA.Logyearly = periodReturn(NVDA$NVDA.Close, period='yearly', type='log')
V.Logyearly = periodReturn(V$V.Close, period='yearly', type='log')
MA.Logyearly = periodReturn(MA$MA.Close, period='yearly', type='log')
WMT.Logyearly = periodReturn(WMT$WMT.Close, period='yearly', type='log')
HD.Logyearly = periodReturn(HD$HD.Close, period='yearly', type='log')
TMO.Logyearly = periodReturn(TMO$TMO.Close, period='yearly', type='log')
COST.Logyearly = periodReturn(COST$COST.Close, period='yearly', type='log')
ABT.Logyearly = periodReturn(ABT$ABT.Close, period='yearly', type='log')
ADBE.Logyearly = periodReturn(ADBE$ADBE.Close, period='yearly', type='log')
DHR.Logyearly = periodReturn(DHR$DHR.Close, period='yearly', type='log')
NKE.Logyearly = periodReturn(NKE$NKE.Close, period='yearly', type='log')
INTC.Logyearly = periodReturn(INTC$INTC.Close, period='yearly', type='log')
AMD.Logyearly = periodReturn(AMD$AMD.Close, period='yearly', type='log')
TMUS.Logyearly = periodReturn(TMUS$TMUS.Close, period='yearly', type='log')
AMGN.Logyearly = periodReturn(AMGN$AMGN.Close, period='yearly', type='log')
LMT.Logyearly = periodReturn(LMT$LMT.Close, period='yearly', type='log')
INTU.Logyearly = periodReturn(INTU$INTU.Close, period='yearly', type='log')
PGR.Logyearly = periodReturn(PGR$PGR.Close, period='yearly', type='log')
ICE.Logyearly = periodReturn(ICE$ICE.Close, period='yearly', type='log')
DG.Logyearly = periodReturn(DG$DG.Close, period='yearly', type='log')


# The average ten-year risk-free treasury rate from 2009 to 2022 is 2.28% (annual)
#https://www.macrotrends.net/2016/10-year-treasury-bond-rate-yield-chart

# I make my own R function
sharpe<-function(logyearly){
  yearly.risk.free = 0.0228
  sharpe.result = (mean(logyearly)-yearly.risk.free)/sd(logyearly)
  return(sharpe.result)
  
}


# compute the annual average return based on 13 years of historical returns for each stock:
mean.logreturns.list = cbind(mean(AAPL.Logyearly), mean(MSFT.Logyearly), mean(GOOGL.Logyearly), mean(AMZN.Logyearly),
                             mean(BRK.B.Logyearly), mean(NVDA.Logyearly), mean(V.Logyearly), mean(MA.Logyearly),
                             mean(WMT.Logyearly), mean(HD.Logyearly), mean(TMO.Logyearly), mean(COST.Logyearly),
                             mean(ABT.Logyearly), mean(ADBE.Logyearly), mean(DHR.Logyearly), mean(NKE.Logyearly),
                             mean(INTC.Logyearly), mean(AMD.Logyearly), mean(TMUS.Logyearly), mean(AMGN.Logyearly),
                             mean(LMT.Logyearly), mean(INTU.Logyearly), mean(PGR.Logyearly), mean(ICE.Logyearly), mean(DG.Logyearly)
)



sharpe.list = cbind(sharpe(AAPL.Logyearly), sharpe(MSFT.Logyearly), sharpe(GOOGL.Logyearly), sharpe(AMZN.Logyearly),
                    sharpe(BRK.B.Logyearly), sharpe(NVDA.Logyearly), sharpe(V.Logyearly), sharpe(MA.Logyearly),
                    sharpe(WMT.Logyearly), sharpe(HD.Logyearly), sharpe(TMO.Logyearly), sharpe(COST.Logyearly),
                    sharpe(ABT.Logyearly), sharpe(ADBE.Logyearly), sharpe(DHR.Logyearly), sharpe(NKE.Logyearly),
                    sharpe(INTC.Logyearly), sharpe(AMD.Logyearly), sharpe(TMUS.Logyearly), sharpe(AMGN.Logyearly),
                    sharpe(LMT.Logyearly), sharpe(INTU.Logyearly), sharpe(PGR.Logyearly), sharpe(ICE.Logyearly), sharpe(DG.Logyearly)
                    )


kurtosis.list = cbind(kurtosis(AAPL.Logyearly), kurtosis(MSFT.Logyearly), kurtosis(GOOGL.Logyearly), kurtosis(AMZN.Logyearly),
                      kurtosis(BRK.B.Logyearly), kurtosis(NVDA.Logyearly), kurtosis(V.Logyearly), kurtosis(MA.Logyearly),
                      kurtosis(WMT.Logyearly), kurtosis(HD.Logyearly), kurtosis(TMO.Logyearly), kurtosis(COST.Logyearly),
                      kurtosis(ABT.Logyearly), kurtosis(ADBE.Logyearly), kurtosis(DHR.Logyearly), kurtosis(NKE.Logyearly),
                      kurtosis(INTC.Logyearly), kurtosis(AMD.Logyearly), kurtosis(TMUS.Logyearly), kurtosis(AMGN.Logyearly),
                      kurtosis(LMT.Logyearly), kurtosis(INTU.Logyearly), kurtosis(PGR.Logyearly), kurtosis(ICE.Logyearly), kurtosis(DG.Logyearly)
)


df<-t(rbind.data.frame(mean.logreturns.list, sharpe.list,kurtosis.list))


colnames(df) = c('13-Year Average of Log Yearly Returns','Sharpe Ratio', 'Excess Kurtosis (normal=0, thinner tails<0, fatter tails>0')
rownames(df) = c('AAPL', 'MSFT', 'GOOGL', 'AMZN', 'BRK.B', 'NVDA', 'V', 'MA', 'WMT', 'HD', 'TMO', 
                          'COST', 'ABT', 'ADBE', 'DHR', 'NKE', 'INTC', 'AMD', 'TMUS', 'AMGN', 'LMT', 'INTU',
                          'PGR', 'ICE', 'DG')


View(df)
# I am going to pick the top 10 sharpe ratios




# The below ten stocks are the top 10 in sharpe ratio (risk-adjusted returns)

# Calculate the vol of the portfolio of the ten selected stocks
# each symbol gets an equal weight of 1/10

#combine ten return series into a matrix

weight=matrix(c(rep(1/10, times=10)),nrow=1)
log.returns = cbind(
  DG.Logyearly,
  COST.Logyearly,
  MA.Logyearly,
  AAPL.Logyearly,
  PGR.Logyearly,
  DHR.Logyearly,
  V.Logyearly,
  HD.Logyearly,
  BRK.B.Logyearly,
  TMO.Logyearly
  
)

#__________________________________________________________________________________________-

# 10 rows for log.returns (matrix 1) times 10 columns transposed of weight (matrix 2)
portfolio.returns = log.returns%*%t(weight) #matrix multiplication; t() is transpose


# *** equal-weighted portfolio returned on average about 17.64% each year from 2009 to 2022
mean(portfolio.returns)

# *** find volatility of portfolio of selected 10 stocks
# portfolio vol is a function of the correlation among the 10 stock components
annualvol_portfolio = sqrt(1)*sd(portfolio.returns) #sqrt(periods per year)*std(portfolio returns) # sqrt(252) for daily data, sqrt(52) for weekly data

# annual volatility of portfolio is only 0.1199
annualvol_portfolio


# Note: 17.64% at 0.1199 equal-weighted portfolio is NOT on the envelope 
#(this can be determined by plotting the dot on the chart and seeing if it lies on envelope convex curve)
# this means that there is another portfolio with higher return for the same variance at 0.1199

vol = cbind(sqrt(1)*sd(DG.Logyearly),sqrt(1)*sd(COST.Logyearly),
            sqrt(1)*sd(MA.Logyearly),sqrt(1)*sd(AAPL.Logyearly), 
            sqrt(1)*sd(PGR.Logyearly), sqrt(1)*sd(DHR.Logyearly),
            sqrt(1)*sd(V.Logyearly),sqrt(1)*sd(HD.Logyearly),
            sqrt(1)*sd(BRK.B.Logyearly),sqrt(1)*sd(TMO.Logyearly),
            annualvol_portfolio
            )

colnames(vol) = c('DG','COST','MA','AAPL','PGR', 'DHR', 'V', 'HD', 'BRK.B', 'TMO', 'My Portfolio')


View(vol)
# Individual names all have higher volatility than portfolio volatility



#___________________________________________________________________________________________

# compare portfolio returns with SP500 via market-based ETF called SPY
getSymbols('SPY', from='2009-6-1', to='2022-6-1', src='yahoo')
SPY.Logyearly = periodReturn(SPY$SPY.Close, period='yearly', type='log')

#average return of SP500 from 2009 to 2022 is 10.51%
mean(SPY.Logyearly)
# volatility (standard deviation) of such log returns is 0.1204
sd(SPY.Logyearly)
#So, my portfolio returns is 1.67 times more than the SP500, 
#and my portfolio has slightly lower volatility than the SP500 (my portfolio vol is 0.1199 vs SPY vol at 0.1204)

# Remember, to get 95% confidence interval: 
#mean(SPY.Logyearly) +/- 1.96 * sd(SPY.Logyearly)/(sqrt(13)) , where 13 years is the number of data points 
#______________________________________________________________________________



# Construct the efficient frontier (the set of points about the global min variance point)
# And find the Global Minimum Variance Portfolio



# Find the mean return of each stock in the portfolio:
mean_vec = apply(log.returns,2,mean, nrow=1)
# mean(mean_vec) is just the average portfolio return of 17.64%
# 2 in the apply function just means column (1 means row)

Sigma = matrix(var(log.returns), nrow=10, byrow=TRUE) #variance-covariance matrix




# to find the two envelopes, you need to first find v and w (see lecture 5, example 1.3, page 12-13)

# find w, with arbitrary c:

# c represents a tangent line on the envelope
# changing the constant c will change the portfolio w weights and hence portfolio w's envelope location... 
# But portfolio w will still lie somewhere on the envelope 
# and regardless of what c is, the SAME envelope is constructed

c=0.0 #choose an arbitrary constant c (0 is often chosen)
mean_vec
Sigma

# minimize the covariance-variance matrix by using solve():
z = solve(Sigma)%*%(mean_vec-c)
w = matrix(z/sum(z),nrow=1) #envelope portfolio
w  # weights (both long stock and short stock) of portfolio w

mean_w = w%*%(mean_vec)
sd_w = sqrt(w%*%Sigma%*%t(w))
mean_w
sd_w


# do same process for v: find v, with arbitrary c:
c=0.04 #choose a different arbitrary constant c
mean_vec
Sigma
z = solve(Sigma)%*%(mean_vec-c)
v = matrix(z/sum(z),nrow=1) #envelope portfolio
v # the weights of portfolio v

mean_v = v%*%(mean_vec)
sd_v = sqrt(v%*%Sigma%*%t(v))
mean_v
sd_v





#####construct efficient frontier using two envelope portfolios that you just found, v and w

# v and w are of the same 10 stocks but different weights for each stock 
#(e.g., portfolio v can have aapl at 56% and portfolio w can have aapl at -16% short weight)


####mean return and standard deviation of w
a = seq(-2.5, 4,by=0.1)  #this is the range of the y axis
envelope_weight = matrix(NA,nrow=length(a),ncol=length(w))
envelope_return=c()
envelope_sd=c()
for (i in 1:length(a)){
  envelope_weight[i,] = a[i]*w+(1-a[i])*v  # a is some constant; w and v are vectors of stock weights
  envelope_return[i] = envelope_weight[i,]%*%(mean_vec)
  envelope_sd[i] =  sqrt(envelope_weight[i,]%*%Sigma%*%t(matrix(envelope_weight[i,],nrow=1)))
}

plot(envelope_sd,envelope_return,xlab='Portfolio Standard Deviation',ylab='Portfolio Return', type='l',lwd=2,col='blue')

points(envelope_sd,envelope_return,pch=16,col='green')

points(sd_w,mean_w,pch=17,col='red')
text(1.1*sd_w,mean_w,'Portfolio w')

points(sd_v,mean_v,pch=15,col='orange')
text(1.02*sd_v,0.95*mean_v,'Portfolio v')
#1.02 and 1.1 just means adjust the chart label to the right of the point
#0.95 just means adjust the chart label to the left of the point


#####find the GMVP
######method 1: use optimizer
library(quadprog)
number_assets = length(w) #number of assets in the portfolio
Dmat = Sigma
dvec = rep(0,number_assets)
Amat = matrix(rep(1,number_assets),ncol=1)
bvec = 1
GMVP_weight = solve.QP(Dmat,dvec,Amat,bvec,meq=1)$solution
GMVP_weight
GMVPweight_vec = matrix(GMVP_weight,nrow=1)


#compute the portfolio mean and standard deviation for GMVP
GMVP_mean = mean_vec%*%t(GMVPweight_vec)
GMVP_sd = sqrt(GMVPweight_vec%*%Sigma%*%t(GMVPweight_vec))

GMVP_mean
GMVP_sd




###plot the portfolio standard deviation versus mean return
#plot(envelope_sd,envelope_return,xlab='Portfolio Standard Deviation',ylab='Portfolio Return', type='l',lwd=2,col='blue',ylim=c(0, 0.16),xlim=c(0.2,0.7))


points(GMVP_sd,GMVP_mean,pch=c('@'),col='purple')
text(1.08*GMVP_sd,GMVP_mean,'GMVP')
points(sqrt(diag(Sigma)),mean_vec,pch=16,col='green')


# you can't actually see below text because they are far to the right of the chart
# the portfolio sd is low and each stock component is high standard deviation
#text(1.1*sqrt(Sigma[1,1]),mean_vec[1],'Asset 1')
#text(1.05*sqrt(Sigma[2,2]),mean_vec[2],'Asset 2')
#text(1.04*sqrt(Sigma[3,3]),mean_vec[3],'Asset 3')
#text(0.98*sqrt(Sigma[4,4]),0.96*mean_vec[4],'Asset 4')


lines(envelope_sd[envelope_return >= as.numeric(GMVP_mean)],envelope_return[envelope_return >= as.numeric(GMVP_mean)],col='red',lwd=2)
legend('bottomright',c('Efficient Frontier'), col=c('red'),lty=1,lwd=2)



# Global Variance Minimum Portfolio 
# The least variance of all portfolios on the envelope

GMVP_weight
GMVP_mean
GMVP_sd


GMVP_table<-(rbind.data.frame(GMVP_weight))
rownames(GMVP_table) = c('Weights for GMVP')
colnames(GMVP_table) = c('DG', 'COST', 'MA', 'AAPL', 'PGR', 'DHR', 'V', 'HD', 'BRK.B', 'TMO')
View(GMVP_table)

