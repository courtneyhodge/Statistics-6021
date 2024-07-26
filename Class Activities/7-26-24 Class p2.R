#######################################################
## File:	Regreg.R
## Purpose: Explore regularized regression 
##
## Author: 	J. Blume
## Date:	July 26, 2025
#######################################################

#######################################################
## Packages
#######################################################

library(glmnet)
library(corrplot)
library(lars)

#######################################################
## Diabetes Data
#######################################################

# Diabetes data (10 variables, 442 measurements) as used in the study of Efron et al. (2004). 
# The data is standardized such that the means of all variables are zero, and all variances are equal to one. 
# https://search.r-project.org/CRAN/refmans/care/html/efron2004.html

data(diabetes)
help(diabetes)

y <- diabetes$y ; length(y)

x <- diabetes$x ; dim(x); 
head(x)
colnames(x)

x2 <- diabetes$x2 ; dim(x2); colnames(x2)

rhomat=cor(cbind("y"=c(y),x))
rhomat2=cor(cbind("y"=c(y),x2))

## Correlation Plot
corrplot(rhomat, method="circle",tl.srt=45,tl.col="black",diag=F) #,addCoef.col = "black")
corrplot(rhomat, method="number",tl.srt=45,tl.col="black",diag=F)

corrplot(rhomat2, method="circle",tl.srt=45,tl.col="black",diag=F) #,addCoef.col = "black")

#######################################################
## Regression
#######################################################

# Main effects model (no interactions)
mod.me <- lm(y~x)
summary(mod.me)
summary(mod.me)$r.squared
summary(mod.me)$sigma

# Check MSE
yhat <- predict(mod.me, newdata=x)
sum((y-yhat)^2)/(length(y)-ncol(x)-1)
(summary(mod.me)$sigma)^2

# Big model
mod.bg <- lm(y~x2)
summary(mod.bg)
summary(mod.bg)$r.squared
summary(mod.bg)$sigma

# check MSE
yhat2 <- predict(mod.bg, newdata=x2)
sum((y-yhat2)^2)/(length(y)-ncol(x2)-1)
(summary(mod.bg)$sigma)^2

# Careful with the predict statement (Example for one at a time)
# predict(mod.me, newdata=data.frame(x=t(x[2,])))

#######################################################
## Ridge regression
#######################################################

lasso.0 <- glmnet(x, y, alpha=0, family="gaussian")

plot(lasso.0, label=TRUE, xvar="lambda")
plot(lasso.0, label=TRUE, xvar="norm") 
plot(lasso.0, label=TRUE, xvar="dev") 

summary(lasso.0)
lasso.0

lnlam=2
coef(lasso.0, s=exp(lnlam))
tmp_coeffs = coef(lasso.0, s=exp(lnlam))

betas = data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], loglam=lnlam, coefficient = tmp_coeffs@x)

plot(lasso.0, label=TRUE, xvar="lambda")
points(betas[-1,2:3], col='purple')

## Where are the OLS parameters
points(x=rep(1.5,10),y=coef(mod.me)[-1], col='blue', col=1:10)

cbind(ridge=coef(lasso.0, s=exp(2)),ols=coef(mod.me))

#### Fun fun -- Add variable names in right margin 
plot(lasso.0)
vn=names(x[1,])
vnat=coef(lasso.0) 
vnat=vnat[-1,ncol(vnat)] # remove the intercept, and get the coefficients at the end of the path 
axis(4, at=vnat,line=-.5,label=vn,las=1,tick=FALSE, cex.axis=0.8)

#######################################################
## Lasso regression
#######################################################

lasso.1 <- glmnet(x, y, alpha=1, family="gaussian")

plot(lasso.1, label=TRUE, xvar="lambda")
plot(lasso.1, label=TRUE, xvar="norm") 
plot(lasso.1, label=TRUE, xvar="dev") 

summary(lasso.1)
lasso.1

lnlam=2
coef(lasso.1, s=exp(lnlam))
tmp_coeffs = coef(lasso.1, s=exp(lnlam))

betas = data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], loglam=lnlam, coefficient = tmp_coeffs@x)

plot(lasso.1, label=TRUE, xvar="lambda")
points(betas[-1,2:3], col='purple')

#######################################################
## Elastic net
#######################################################

lasso.5 <- glmnet(x, y, alpha=0.5, family="gaussian")	

plot(lasso.5, label=TRUE, xvar="lambda")
plot(lasso.5, label=TRUE, xvar="norm") 
plot(lasso.5, label=TRUE, xvar="dev") 

summary(lasso.5)
lasso.5

lnlam=2
coef(lasso.5, s=exp(lnlam))
tmp_coeffs = coef(lasso.5, s=exp(lnlam))

betas = data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], loglam=lnlam, coefficient = tmp_coeffs@x)

plot(lasso.5, label=TRUE, xvar="lambda")
points(betas[-1,2:3], col='purple')

## Change alpha and see how smoother the lines get

#######################################################
## Force variable(s) to remain in model 
#######################################################

## variable order
cbind(names(x[1,]),seq(dim(x)[[2]]))

## Set penalty vector (1=normal lasso, 0=Keep, huge=Drop)
keep=rep(1,dim(x)[[2]])
keep[3]=0					# Force BMI [3rd column in X] in model

lasso.bmi <- glmnet(x, y, alpha=1, penalty.factor=keep)	

plot(lasso.bmi, label=TRUE, xvar="lambda", main="Force BMI in model")

dev.new()
plot(lasso.1, label=TRUE, xvar="lambda", main="Just Lasso")
dev.off()

## Check Lasso against OLS
coef(lasso.bmi,s=max(lasso.bmi$lambda))
lm(y~x[,3])

coef(lasso.bmi,s=min(lasso.bmi$lambda))
lm(y~x) 

#######################################################
## Fully Relaxed Lasso (Lasso for selection; OLS for Fit)
#######################################################

b <- coef(lasso.1, s=exp(2))
colnames(b) <- "lasso.1"
subset <- colnames(x)[which(b!=0)-1]
dropped <- setdiff(colnames(x),subset)

ols.1 <- lm(y~x[,subset])
summary(ols.1)		# this is a "Fully Relaxed Lasso" 

b.ols <- as(c(coef(ols.1),rep(0,length(dropped))),"sparseMatrix")
rownames(b.ols)=c("(Intercept)",subset,dropped)

cbind(b,'ols'=b.ols[match(rownames(b),rownames(b.ols))])

#######################################################
## Crossvalidation 
#######################################################

## Lasso
l1.cv <- cv.glmnet(x, y, alpha=1, family="gaussian")
summary(l1.cv)

plot(l1.cv)  # top axis is number of variables in model

plot(lasso.1,label=TRUE,xvar="lambda")
abline(v=log(l1.cv$lambda.min),lty=2,lwd=0.5)
abline(v=log(l1.cv$lambda.1se),lty=2,lwd=0.5)

cbind(coef(l1.cv,s="lambda.min"),coef(l1.cv,s="lambda.1se"))

## Fully Relaxed Lasso (Cross-validated)
b <- coef(l1.cv,s="lambda.1se")
colnames(b) <- "lasso.1"
subset <- colnames(x)[which(b!=0)-1]
dropped <- setdiff(colnames(x),subset)

ols.1 <- lm(y~x[,subset])
summary(ols.1)		# this is a "Fully Relaxed Lasso" 

#######################################################
## Bootstrap Lasso
#######################################################

## Examine the "stability" of the regularized model

boots=500
b.betas <- matrix(NA,nrow=boots,ncol=dim(x)[[2]]+1)

for (i in 1:boots) {
  select <- sample(length(y),replace=TRUE)
  l1.cv.boot <- cv.glmnet(x[select,],y[select],alpha=1,family="gaussian")
  b.betas[i,]=as.vector(coef(l1.cv.boot,s="lambda.min"))
  if(i %% 50==0) {cat(paste0("iteration: ", i, "\n"))}
}

colMeans(b.betas)
hist(b.betas[,2])
hist(b.betas[,2],breaks=80)	

#######################################################
## Exercise
#######################################################

1. Repeat for X2 instead of X

2. Repeart using Boston Housing Data. (predict 'medv' = median value of owner-occupied homes in USD 1000s )
## https://search.r-project.org/CRAN/refmans/mlbench/html/BostonHousing.html

data("BostonHousing", package = "mlbench")
summary(BostonHousing)

data("BostonHousing2", package = "mlbench")
summary(BostonHousing2)


####
###
##
#