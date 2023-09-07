# Load the necessary packages:
library(tseries);library(dplyr);library(ggplot2);library(fGarch)
library(zoo);library(rugarch);library(forecast);library(Ecdat)
library(ggpubr); library(PerformanceAnalytics)
# Load the data
df <- get.hist.quote(instrument="^FTMC", start="2000-01-01",end="2022-02-15",
                     quote="Adj",retclass = 'zoo')
# Save the dates
Date <- index(df)
# Convert data into a data frame
if (is.data.frame(df) == FALSE) {
  df = as.data.frame(df)
}
# Create a data frame with Dates and Prices
df <- data.frame('Date'=Date,'Price'=df$Adjusted)
head(df,5)
# Omit the NA values from the data
df   <- df %>% filter(!is.na(Price) & !is.na(Date))
Date <- df %>% select(Date)
# Summary of the data:
df %>% select(Price) %>% summary() 
#Examination of Stationarity 
# Plot of the FTSE 250 stock
g1 <- ggplot(df,aes(x=Date,y=Price)) +
       geom_line(color='darkred') +
        labs(x='Time',y='FTSE 250',
          title='Plot of FTSE 250 stock during T = 01-01-2000 - 15-02-2022')+
            theme_bw();g1 # This is like a trend.
# Convert to stationarity data  
adf.test(df$Price,alternative = 'stationary',k=10)$p.val
#Transformations for Stationarity 
# Square root
sq <- apply(df %>% select(Price),2,sqrt) #sqroo
sq <- diff(sq); adf.test(sq,k=15)$p.val
g2 <- ggplot(data.frame('Date'=Date[-1,],'Price'=sq),aes(x=Date,y=Price)) +
    geom_line(color='darkred') +
    labs(x='Time',y='Square root differences',
    title='Plot of square root differences of prices during T')+
    theme_bw() # This is like a trend.

# Cube root
cb <- apply(df %>% select(Price),2,function(x){x^(1/3)})
cb <- diff(cb); adf.test(cb,k=15)$p.val 
g3 <- ggplot(data.frame('Date'=Date[-1,],'Price'=cb),aes(x=Date,y=Price)) +
    geom_line(color='darkred') +
    labs(x='Time',y='Cube root differences',
    title='Plot of cube root differences of prices during T')+
    theme_bw() # This is like a trend.
ggarrange(g2,g3,nrow=2,ncol=1)

# Log transformation
rt <- rep(NA,dim(df)[1]-1)
for(i in 1:(dim(df)[1]-1)){
  rt[i] <- log(df[i+1,2]) - log(df[i,2]) 
}
rt <- data.frame('Return'=rt, 'Square.return' = (rt)^2)
rownames(rt) <- Date[-1,];summary(rt %>% select(1)) # Summary of log-returns
# Plot of the log returns, squared returns and absolute returns of FTSE 250:
g4 <-  ggplot(data.frame('Date' = Date[-1,],'Return'=rt$Return),
      aes(x=Date,y=Return)) +
      geom_line(color='darkred') +
      labs(x='Time',y='log-returns FTSE 250',
      title='Plot of log-returns during T')+
      theme_bw()
g5 <-  ggplot(data.frame('Date' = Date[-1,],'Square.return'=rt$Square.return),
      aes(x=Date,y=Square.return)) +
      geom_line(color='darkred') +
      labs(x='Time',y='squared log-returns FTSE 250',
      title='Plot of squared log returns during T')+
      theme_bw()
g6 <-  ggplot(data.frame('Date' = Date[-1,],'abs.return'=abs(rt$Return)),
      aes(x=Date,y=abs.return)) +
      geom_line(color='darkred') +
      labs(x='Time',y='absolute log-returns FTSE 250',
      title='Plot of absolute log-returns during T')+
      theme_bw()
ggarrange(g4,g5,g6,nrow = 3,ncol = 1) # Combination of ggplots
# Histogram of the Log-returns with the N(0,1) density curve:
g7 <- rt %>% 
      ggplot() +
      geom_histogram(aes(x=Return,y=..density..),color='black',fill='white') +
      stat_function(fun=dnorm, args = list(mean = mean(rt$Return), 
      sd = sd(rt$Return)), size=1,col='darkred')+ labs(
      x='Log-returns of FTSE 250',
      title='Histogram of log-returns during T')+
      theme_bw();g7
#Fitting ARMA-GARCH MODEL 
# log returns,square log returns, and absolute log returns:
par (mfrow=c(2,2))
acf (rt  %>% select(1),main='ACF of FTSE 250 stock log-returns')
pacf(rt %>% select(1),main='PACF of FTSE 250 stock log-returns')
acf (rt  %>% select(2),main='ACF of FTSE 250 stock squared log-returns')
acf (abs(rt %>% select(1)),main='ACF of FTSE 250 stock absolute log-returns')
par (mfrow=c(1,1))
# MODEL SELECTION: (We suppose 12*12 models)
garch_order <- matrix(NA,12,2)
for(i in 1:3){
  garch_order[(4*i-3):(4*i),1] <- rep(i,4) 
  garch_order[(4*i-3):(4*i),2] <- seq(0,3,by=1)
             }
garch_order  ;arma_order <- garch_order

bayes_crit <- matrix(NA,12,12)
for(j in 1:dim(garch_order)[1]) {
  for(i in 1:dim(arma_order)[1])  {
    fit.model <- ugarchfit(spec = ugarchspec(variance.model=list(
              model="sGARCH", garchOrder=garch_order[j,]),
              mean.model=list(armaOrder=arma_order[i,]),
              distribution.model="norm"),data = rt$Return,solver = 'hybrid')
     bayes_crit[i,j] <- infocriteria(fit.model)[2] # With the BIC criterion
                                  }
                                }
index.garch_order <- 
  which(apply(bayes_crit,2,min)==min(apply(bayes_crit,2,min)))
index.arma_order  <- 
  which(bayes_crit[,index.garch_order]==min(apply(bayes_crit,2,min)))
# Best model:
garch.fit <- garch_order[index.garch_order,]
arma.fit  <- arma_order[index.arma_order,]
# We conclude that according to BIC criterion the model is:
# GARCH(1,1)-ARMA(1,0)
ftse.model <- ugarchfit(spec = ugarchspec(variance.model=list(
              model="sGARCH",garchOrder=garch.fit),
              mean.model=list(armaOrder=arma.fit),distribution.model="norm"),
              data = rt$Return,solver = 'hybrid')
ftse.model;coef(ftse.model) # Coefficients of the ARMA(1,0)-GARCH(1,1)

#Testing for auto correlations of the residuals and square residuals
### Residuals
# Box.test(*,type='Ljung',fitdf=1+0):
res <- ftse.model@fit$z
Box.test(res,type = 'Ljung',lag = 10,fitdf = 1)$p.val
Box.test(res,type = 'Ljung',lag = 20,fitdf = 1)$p.val
Box.test(res,type = 'Ljung',lag = 40,fitdf = 1)$p.val
### Square Residuals:
# Box.test(*,type='Ljung',fitdf=1+1):
Box.test(res^2,type = 'Ljung',lag = 10,fitdf = 2)$p.val
Box.test(res^2,type = 'Ljung',lag = 20,fitdf = 2)$p.val
Box.test(res^2,type = 'Ljung',lag = 40,fitdf = 2)$p.val

#Testing the Normality in residuals
### Histogram of standardized residuals
g8 <- ggplot(data.frame('Residuals' = res)) +
    geom_histogram(aes(x=Residuals,y=..density..),
    color='black',fill='white') +
    stat_function(fun=dnorm, args = list(mean = mean(res),sd = sd(res)),
    size=1,col='darkred') + labs(x='Residuals')+  
    ggtitle('Density of Standardized Residuals and N(0,1) curve') +
    theme_bw();g8 
### QQ-plot of standardized residuals:
qqnorm(res,pch=20,
main = expression("Q-Q plot for" ~~ {'Standardized residuals'}))
qqline(res,col='darkred',lwd=2);mtext("qqline(*, dist = qnorm)")

### Tests for normality:
# kolmogorov-smirnov, Jarque-bera test
fun <- function(x){
  c(round(ks.test(x,'pnorm',mean=mean(res),sd = sd(res))$p.val,4),# ks.test
    round(jarque.bera.test(x)$p.val,4))# jarque-bera test
};res %>% fun 

#Fitting t-distribution for residuals
par(mfrow=c(1,2)) # We will use the function stdFit:
res.t.fit <- stdFit(res)
dof <- as.numeric(res.t.fit$par[3])
qqplot(res,qt(ppoints(10000),df=dof), ylab = 'Sample Quantiles',
  xlab = 'Theoretical Quantiles',
  main = expression("Q-Q plot for" ~~ {t}[nu = 8.061278 ]),pch=20)
qqline(res,distribution = function(y,df=dof) qt(y,df=dof), 
  col = 'darkred',lwd=2);mtext("qqline(*, dist = qt(., df=8.061278)")

#Fitting skew t-distribution for residuals
res.st.fit <- sstdFit(res)
param <- round(res.st.fit$estimate,3)
set.seed(1235)
qqplot(x=rsstd(10000,mean = param[1],
  sd = param[2],nu = param[3],xi = param[4]),
  y = res, main = expression("Q-Q plot for" ~~ {st}[nu = (8.527)]),pch=20,
  ylab = 'Sample Quantiles', xlab = 'Theoretical Quantiles')
qqline(res,
 distribution = function(p) qsstd(p,mean=param[1],sd=param[2],nu=param[3], 
                                  xi=param[4]),lwd=2,col='darkred')
mtext("qqline(*, dist = qsstd(., mean=-0.044,sd=0.998,nu=8.527,xi=0.860)")
par(mfrow=c(1,1))

# Skew t distribution is a better fit than normal and t distribution.
ftse.model <- ugarchfit(spec = ugarchspec(variance.model=list(
              model="sGARCH", garchOrder=garch.fit),
              mean.model=list(armaOrder=arma.fit),distribution.model="sstd"),
              data = rt %>% select(1),solver = 'hybrid')
ftse.model                                    
BIC.sgarch <- infocriteria(ftse.model)[2]  
# Normal Vs Skew-Student for Residuals.
plot(ftse.model,which=8)
coef(ftse.model)
# Extensions of Garch processes
BIC.egarch <- infocriteria(ugarchfit(spec = ugarchspec(variance.model=list(
              model="eGARCH",garchOrder=garch.fit),
              mean.model=list(armaOrder=arma.fit),distribution.model="sstd"),
              data = rt %>% select(1),solver = 'hybrid'))[2]
BIC.aparch <- infocriteria(ugarchfit(spec = ugarchspec(variance.model=list(
              model="apARCH",garchOrder=garch.fit),
              mean.model=list(armaOrder=arma.fit),distribution.model="sstd"),
              data = rt %>% select(1),solver = 'hybrid'))[2]
BIC.igarch <- infocriteria(ugarchfit(spec = ugarchspec(variance.model=list(
              model="iGARCH",garchOrder=garch.fit),
              mean.model=list(armaOrder=arma.fit),distribution.model="sstd"),
              data = rt %>% select(1),solver = 'hybrid'))[2]
name_model <- c('sGARCH','iGARCH','eGARCH','apARCH')
method     <- c(BIC.sgarch,BIC.igarch,BIC.egarch,BIC.aparch)
# Create a matrix with the BIC values according to different Garch processes:
MODEL_SELECTION <- matrix(method,nrow=1,ncol=4)
colnames(MODEL_SELECTION)  <- name_model
row.names(MODEL_SELECTION) <- 'BIC'
MODEL_SELECTION
# Find the model with the lowest BIC value:
min(MODEL_SELECTION)
# So, we have a better fit for the problem which is the apARCH(1,1):
name_model[which(MODEL_SELECTION == min(MODEL_SELECTION))]
# This suggests that the AR(1)-APARCH(1,1) model is the one that
# fits the FTSE 250 time series best.
best.fit <- ugarchfit(spec = ugarchspec(variance.model=list(
            model="apARCH",garchOrder=garch.fit),
            mean.model=list(armaOrder=arma.fit),distribution.model="sstd"),
            data = rt %>% select(1),solver = 'hybrid');best.fit
# Leverage effect
# Method:
col_1 <- rep(NA,dim(rt)[1]); col_2 <- rep(NA,dim(rt)[1])
for(i in 1:dim(rt)[1]) {
  col_1[i] <- rt$Return[i] 
  col_2[i] <- rt$Square.return[i+1] }
leverage <- data.frame('ret'=col_1[-dim(rt)[1]],'s.ret'=col_2[-dim(rt)[1]])
cor(leverage$ret,leverage$s.ret)

# Another method:
r.neg <- abs(rt$Return[which(rt$Return < 0)]) # Absolute negative log-returns
r.pos <-     rt$Return[which(rt$Return > 0)]  # Positive log-returns
# Plot of the relationship between +ve and absolute -ve log-returns and 
# Leverage effect plot(QQ-plots between +ve and absolute -ve log-returns). 
set.seed(565)
par(mfrow=c(1,2))
plot(x = r.neg, y = sample(r.pos,length(r.neg)), pch = 20,
    xlab = '-ve absolute log Returns',ylab='+ve log Returns',
    main='Plot of +ve and -ve absolute log Returns')
qqplot(x = sample(r.pos,length(r.neg)), y = r.neg,
       xlab = 'quantiles of log returns after +ve returns',
       ylab = 'quantiles of absolute log returns after -ve returns',
      pch = 20, lwd = 2, col = 'black', main = 'Leverage effect plot')
abline(lm(sort(r.neg) ~ sort(sample(r.pos,length(r.neg)))), 
       col = 'darkred', lwd = 2) ;par(mfrow=c(1,1))
# Conditional Volatility plot
sigma(best.fit) # Sigma values
head(sigma(best.fit),20)
# Plot of the Volatility 
plot(sigma(best.fit), ylab="sigma(t)", col="darkred", 
        main = "Conditional Variance of ARMA(1,0)-apARCH(1,1) Model")
# -------------------------------------------
# END: Vasileios Diplas
# -------------------------------------------