#Packages needed for the simulation to run
install.packages("MASS")
library(MASS)
install.packages("plyr")
library(plyr)

#Generating the data
set.seed(98)
data <- mvrnorm(1000, mu = c(0,0,0,0), Sigma = matrix(c(1,0,0.5,0.5,0,1,0,0,0.5,0,1,0,0.5,0,0,1), ncol = 4),empirical = TRUE)
X1 = data[,1]
X2 = data[,2]
X3 = data[,3]
Z = data[,4]
e = rnorm(1000, 0, 10)
Y = 2*X1 + X2 + X3 + e
orig = data.frame(Y,X1,X2,X3,Z,e)
head(orig)

#Sanity check to make sure the correlations are what we intended for them to be
cor(X1, X2)
cor(X1, X3)
cor(X1, Z)
cor(X2, X3)
cor(X2, Z)
cor(X3, Z)

#Functions for our simulation to run
#this function removes the data at a specific rate (what percent of entries we want removed from X1) for a given specified type
#this will allow us to first remove the data in a certain way to then refill it in with the three imputation techniques
remove.data = function(orig,type="MCAR",rate=0.25){
  if(type=="MCAR"){
    inds = sample(1000,rate*1000)
    orig$X1[inds]=NA
  }
  if(type=="MAR"){
    rand=rbinom(1000,prob=rate*(2/5)+rate*(6/5)*(orig$X3 >= median(orig$X3)),1)
    orig$X1[rand==1]=NA
  }
  if(type=="MNAR"){
    rand=rbinom(1000,prob=rate*(2/5)+rate*(6/5)*(orig$Z >= median(orig$Z)),1)
    orig$X1[rand==1]=NA
  }
  return(newdata=orig)
}

#this function specifies which type of imputation we would like to refill the dataset with
impute.data=function(newdata,type="mean"){
  if(type=="mean"){
    newdata$X1[is.na(newdata$X1)] = mean(newdata$X1, na.rm = TRUE)
  }
  if(type=="hotdeck"){
    newdata$X1[is.na(newdata$X1)] = sample(newdata$X1[!is.na(newdata$X1)], 1, replace=TRUE)
  }
  if(type=="regress"){
    newdata.full = newdata[complete.cases(newdata),]
    newdata.miss= newdata[!complete.cases(newdata),]
    regress = lm(X1~X2+X3, data=newdata.full)
    newdata$X1[is.na(newdata$X1)]=predict(regress,newdata=newdata.miss)
  }
  return(imputed.data=newdata)
}

#we will run 10,000 simulations so that we can analyze the results beta1 values
nsims=10000

#creating tables to store the beta values according to type of missingness and imputation technique used
mcar.results=matrix(NA,nrow=nsims,ncol=3)
colnames(mcar.results) = c('mean', 'hotdeck', 'regress')

mar.results=matrix(NA,nrow=nsims,ncol=3)
colnames(mar.results) = c('mean', 'hotdeck', 'regress')

mnar.results=matrix(NA,nrow=nsims,ncol=3)
colnames(mnar.results) = c('mean', 'hotdeck', 'regress')

#creating a confidence interval dataframe to store the intervals for each simulation
mcar.confint=matrix(NA,nrow=nsims,ncol=6)
colnames(mcar.confint) = c('meanCIlow','meanCIhigh', 'hotdeckCIlow','hotdeckCIhigh', 'regressCIlow','regressCIhigh')

mar.confint=matrix(NA,nrow=nsims,ncol=6)
colnames(mar.confint) = c('meanCIlow','meanCIhigh', 'hotdeckCIlow','hotdeckCIhigh', 'regressCIlow','regressCIhigh')

mnar.confint=matrix(NA,nrow=nsims,ncol=6)
colnames(mnar.confint) = c('meanCIlow','meanCIhigh', 'hotdeckCIlow','hotdeckCIhigh', 'regressCIlow','regressCIhigh')

#We will now run the simulation for each type of missingness, in which the different imputation techniques will all be employed

#MCAR
for(i in 1:nsims){
  missing.data=remove.data(orig,"MCAR",rate=0.2)
  filled.data1=impute.data(missing.data,type="mean")
  fit1=lm(Y~X1+X2+X3, data=filled.data1)
  filled.data2=impute.data(missing.data,type="hotdeck")
  fit2=lm(Y~X1+X2+X3, data=filled.data2)
  filled.data3=impute.data(missing.data,type="regress")
  fit3=lm(Y~X1+X2+X3, data=filled.data3)
  mcar.results[i,1]=coef(fit1)["X1"]
  mcar.results[i,2]=coef(fit2)["X1"]
  mcar.results[i,3]=coef(fit3)["X1"]
  mcar.confint[i,1]=confint(fit1)[2,][1]
  mcar.confint[i,2]=confint(fit1)[2,][2]
  mcar.confint[i,3]=confint(fit2)[2,][1]
  mcar.confint[i,4]=confint(fit2)[2,][2]
  mcar.confint[i,5]=confint(fit3)[2,][1]
  mcar.confint[i,6]=confint(fit3)[2,][2]
}
mcar.test=cbind(mcar.results,mcar.confint)

#MAR
for(i in 1:nsims){
  missing.data=remove.data(orig,"MAR",rate=0.2)
  filled.data1=impute.data(missing.data,type="mean")
  fit1=lm(Y~X1+X2+X3, data=filled.data1)
  filled.data2=impute.data(missing.data,type="hotdeck")
  fit2=lm(Y~X1+X2+X3, data=filled.data2)
  filled.data3=impute.data(missing.data,type="regress")
  fit3=lm(Y~X1+X2+X3, data=filled.data3)
  mar.results[i,1]=coef(fit1)["X1"]
  mar.results[i,2]=coef(fit2)["X1"]
  mar.results[i,3]=coef(fit3)["X1"]
  mar.confint[i,1]=confint(fit1)[2,][1]
  mar.confint[i,2]=confint(fit1)[2,][2]
  mar.confint[i,3]=confint(fit2)[2,][1]
  mar.confint[i,4]=confint(fit2)[2,][2]
  mar.confint[i,5]=confint(fit3)[2,][1]
  mar.confint[i,6]=confint(fit3)[2,][2]
}
mar.test=cbind(mar.results,mar.confint)

##MNAR
for(i in 1:nsims){
  missing.data=remove.data(orig,"MNAR",rate=0.2)
  filled.data1=impute.data(missing.data,type="mean")
  fit1=lm(Y~X1+X2+X3, data=filled.data1)
  filled.data2=impute.data(missing.data,type="hotdeck")
  fit2=lm(Y~X1+X2+X3, data=filled.data2)
  filled.data3=impute.data(missing.data,type="regress")
  fit3=lm(Y~X1+X2+X3, data=filled.data3)
  mnar.results[i,1]=coef(fit1)["X1"]
  mnar.results[i,2]=coef(fit2)["X1"]
  mnar.results[i,3]=coef(fit3)["X1"]
  mnar.confint[i,1]=confint(fit1)[2,][1]
  mnar.confint[i,2]=confint(fit1)[2,][2]
  mnar.confint[i,3]=confint(fit2)[2,][1]
  mnar.confint[i,4]=confint(fit2)[2,][2]
  mnar.confint[i,5]=confint(fit3)[2,][1]
  mnar.confint[i,6]=confint(fit3)[2,][2]
}
mnar.test=cbind(mnar.results,mnar.confint)

## Analysis of the beta1 values

#Variance
var.mcar=c(var(mcar.results[,1]), var(mcar.results[,2]),var(mcar.results[,3]))
var.mar=c(var(mar.results[,1]),var(mar.results[,2]),var(mar.results[,3]))
var.mnar=c(var(mnar.results[,1]),var(mnar.results[,2]),var(mnar.results[,3]))

variances = data.frame(var.mcar,var.mar,var.mnar)
colnames(variances) = c('MCAR', 'MAR', 'MNAR')
round(variances,3)

#Bias
bias.mcar=mcar.results-2
bias.mnar=mnar.results-2
bias.mar=mar.results-2

mcar.bias=abs(round(c(mean(bias.mcar[,1]),mean(bias.mcar[,2]),mean(bias.mcar[,3]))/2*100,2))
mar.bias=abs(round(c(mean(bias.mar[,1]),mean(bias.mar[,2]),mean(bias.mar[,3]))/2*100,2))
mnar.bias=abs(round(c(mean(bias.mnar[,1]),mean(bias.mnar[,2]),mean(bias.mnar[,3]))/2*100,2))

bias.percent=rbind(mcar.bias,mar.bias,mnar.bias)
colnames(bias.percent) = c('mean', 'hotdeck', 'regress')
rownames(bias.percent) = c('MCAR', 'MAR', 'MNAR')

#Visual representation: barplot for percent biased
barplot(bias.percent,beside=TRUE,legend.text=TRUE,
        main="20% of X1 entries removed",xlab="Imputation Techniques", 
        ylab="Percent Biased",col=c("purple","lavender","blue"), 
        args.legend = list(x=ncol(coverageRatio)+12,y=50),xlim=c(0,14),ylim=c(0,100))


#Creating confidence intervals for each beta1, and then calculating coverage for each type of missingness and imputation

#MCAR confidence intervals
ci.mcar.mean=cbind(mcar.test[,4],mcar.test[,5])
ci.mcar.hotdeck=cbind(mcar.test[,6],mcar.test[,7])
ci.mcar.regress=cbind(mcar.test[,8],mcar.test[,9])

#Calculationg coverage for MCAR:
#creating variables to store results
coverageMCAR.mean = rep(NA,length(ci.mcar.mean[,1]))
coverageMCAR.hotdeck = rep(NA,length(ci.mcar.hotdeck[,1]))
coverageMCAR.regress = rep(NA,length(ci.mcar.regress[,1]))

#Checking each interval for each imputation technique to see if 2 falls within it
for (i in 1:length(ci.mcar.mean[,1])) {
  if (ci.mcar.mean[i,1] > 2 || ci.mcar.mean[i,2] < 2) {
    coverageMCAR.mean[i]=F
  }
  else coverageMCAR.mean[i]=T 
}

for (i in 1:length(ci.mcar.hotdeck[,1])) {
  if (ci.mcar.hotdeck[i,1] > 2 || ci.mcar.hotdeck[i,2] < 2) {
    coverageMCAR.hotdeck[i]=F
  }
  else coverageMCAR.hotdeck[i]=T 
}

for (i in 1:length(ci.mcar.regress[,1])) {
  if (ci.mcar.regress[i,1] > 2 || ci.mcar.regress[i,2] < 2) {
    coverageMCAR.regress[i]=F
  }
  else coverageMCAR.regress[i]=T 
}

#Calculating the coverage ratio for each imputation technique
cov.ratio.mcar.mean = sum(coverageMCAR.mean)/length(coverageMCAR.mean)
cov.ratio.mcar.hotdeck = sum(coverageMCAR.hotdeck)/length(coverageMCAR.hotdeck)
cov.ratio.mcar.regress = sum(coverageMCAR.regress)/length(coverageMCAR.regress)

#MAR Confidence Intervals
ci.mar.mean=cbind(mar.test[,4],mar.test[,5])
ci.mar.hotdeck=cbind(mar.test[,6],mar.test[,7])
ci.mar.regress=cbind(mar.test[,8],mar.test[,9])

#Calculationg coverage for MAR:
#creating variables to store results
coverageMAR.mean = rep(NA,length(ci.mar.mean[,1]))
coverageMAR.hotdeck = rep(NA,length(ci.mar.hotdeck[,1]))
coverageMAR.regress = rep(NA,length(ci.mar.regress[,1]))

#Checking each interval for each imputation technique to see if 2 falls within it
for (i in 1:length(ci.mar.mean[,1])) {
  if (ci.mar.mean[i,1] > 2 || ci.mar.mean[i,2] < 2) {
    coverageMAR.mean[i]=F
  }
  else coverageMAR.mean[i]=T 
}

for (i in 1:length(ci.mar.hotdeck[,1])) {
  if (ci.mar.hotdeck[i,1] > 2 || ci.mar.hotdeck[i,2] < 2) {
    coverageMAR.hotdeck[i]=F
  }
  else coverageMAR.hotdeck[i]=T 
}

for (i in 1:length(ci.mar.regress[,1])) {
  if (ci.mar.regress[i,1] > 2 || ci.mar.regress[i,2] < 2) {
    coverageMAR.regress[i]=F
  }
  else coverageMAR.regress[i]=T 
}

#Calculating the coverage ratio for each imputation technique
cov.ratio.mar.mean = sum(coverageMAR.mean)/length(coverageMAR.mean)
cov.ratio.mar.hotdeck = sum(coverageMAR.hotdeck)/length(coverageMAR.hotdeck) 
cov.ratio.mar.regress = sum(coverageMAR.regress)/length(coverageMAR.regress)

#MNAR confidence intervals
ci.mnar.mean=cbind(mnar.test[,4],mnar.test[,5])
ci.mnar.hotdeck=cbind(mnar.test[,6],mnar.test[,7])
ci.mnar.regress=cbind(mnar.test[,8],mnar.test[,9])

#Calculationg coverage for MNAR:
#creating variables to store results
coverageMNAR.mean = rep(NA,length(ci.mnar.mean[,1]))
coverageMNAR.hotdeck = rep(NA,length(ci.mnar.hotdeck[,1]))
coverageMNAR.regress = rep(NA,length(ci.mnar.regress[,1]))

#Checking each interval for each imputation technique to see if 2 falls within it
for (i in 1:length(ci.mnar.mean[,1])) {
  if (ci.mnar.mean[i,1] > 2 || ci.mnar.mean[i,2] < 2) {
    coverageMNAR.mean[i]=F
  }
  else coverageMNAR.mean[i]=T 
}

for (i in 1:length(ci.mnar.hotdeck[,1])) {
  if (ci.mnar.hotdeck[i,1] > 2 || ci.mnar.hotdeck[i,2] < 2) {
    coverageMNAR.hotdeck[i]=F
  }
  else coverageMNAR.hotdeck[i]=T 
}

for (i in 1:length(ci.mnar.regress[,1])) {
  if (ci.mnar.regress[i,1] > 2 || ci.mnar.regress[i,2] < 2) {
    coverageMNAR.regress[i]=F
  }
  else coverageMNAR.regress[i]=T 
}

#Calculating the coverage ratio for each imputation technique
cov.ratio.mnar.mean = sum(coverageMNAR.mean)/length(coverageMNAR.mean)
cov.ratio.mnar.hotdeck = sum(coverageMNAR.hotdeck)/length(coverageMNAR.hotdeck) 
cov.ratio.mnar.regress = sum(coverageMNAR.regress)/length(coverageMNAR.regress)

#Creating one large table with all the coverage ratios
coverageRatio = rbind(c(cov.ratio.mcar.mean,cov.ratio.mcar.hotdeck,cov.ratio.mcar.regress),
                      c(cov.ratio.mar.mean,cov.ratio.mar.hotdeck,cov.ratio.mcar.regress),
                      c(cov.ratio.mnar.mean,cov.ratio.mnar.hotdeck,cov.ratio.mnar.regress))
colnames(coverageRatio) = c('mean', 'hotdeck', 'regress')
rownames(coverageRatio) = c('MCAR', 'MAR', 'MNAR')
coverageRatio

#Visual Representation: barplot for coverage
barplot(coverageRatio,beside=TRUE,legend.text=TRUE,
        main="20% of X1 entries removed",xlab="Imputation Techniques", 
        ylab="Coverage",col=c("purple","lavender","blue"), 
        args.legend = list(x=ncol(coverageRatio)+12,y=0.5),xlim=c(0,14),ylim=c(0,1))

#remember must change rate and titles on barplots each time simulation is run                                                                                                                     