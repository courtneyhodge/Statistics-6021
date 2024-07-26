#######################################################
## File:	Regex.R
## Purpose: Regression exploration
##
## Author: 	J. Blume
## Date:	July 26, 2025
#######################################################

#### Data !
n=1000
b0=5
b1=1

x = rnorm(n, mean=20 , sd=4) 
#x.seq=seq(min(x), max(x), 0.1)

y = b0 + b1*x + rnorm(n, mean=0, sd=6) #adds noise 

#mean of zero above is an assumption

mod.slr = lm(y ~ x)
coef(mod.slr)

#### Plot

plot(x, y, type="n")
points(x, y, col="black", pch=20)
abline(a=b0,b=b1,lty=1,lwd=2, col="black") #true population line
abline(lm(y~x),lty=2, lwd=2, col="red") #estimator model


#### Compute Mean Square Error

yhat = coef(mod.slr)[1] + coef(mod.slr)[2]*x
mse = sum((y-yhat)^2)/length(y)

find.mse = function(y, x, a, b){ #a and b are estimates
  yhat = a + b*x
  mse = sum((y-yhat)^2)/length(y)
  return(mse)
}

find.mse(y=y, x=x, a=coef(mod.slr)[1], b=coef(mod.slr)[2])
find.mse(y=y, x=x, a=coef(mod.slr)[1]+0.5, b=coef(mod.slr)[2]+0.5)  #added a half to intercept and slope just because

## start simply
a=seq( -20, 20, 1) #yintercept
b=seq( -20, 20, 1) #slope


#a=seq( -20, 20, 0.05)
#b=seq( -20, 20, 0.05)

## Find minimum

object.plain = matrix(NA, nrow=length(a), ncol=length(b)) # a big ole matrix, where each row is a column is a slope
object.ridge = matrix(NA, nrow=length(a), ncol=length(b))
object.lasso = matrix(NA, nrow=length(a), ncol=length(b))

for (i in 1:length(a)) {
  for (j in 1:length(b)) {
    object.plain[i,j] = find.mse(y=y, x=x, a=a[i], b=b[j]) 
    object.ridge[i,j] = find.mse(y=y, x=x, a=a[i], b=b[j]) + 0.2*(a[i]^2 + b[j]^2)
    object.lasso[i,j] = find.mse(y=y, x=x, a=a[i], b=b[j]) + 0.2*(abs(a[i]) + abs(b[j]))
    if((i %% 50==0) & (j==1)) {cat(paste0("iteration: ", i, "\n"))}
  }
}

## Standard
#object.plain
min(object.plain)
soln.plain=which(object.plain == min(object.plain), arr.ind = TRUE)
soln.plain
a[soln.plain[1]] ; b[soln.plain[2]]

## Ridge
#object.ridge
min(object.ridge)
soln.ridge=which(object.ridge == min(object.ridge), arr.ind = TRUE)
soln.ridge
a[soln.ridge[1]] ; b[soln.ridge[2]]

## Lasso
#object.lasso
min(object.lasso)
soln.lasso=which(object.lasso == min(object.lasso), arr.ind = TRUE)
soln.lasso
a[soln.lasso[1]] ; b[soln.lasso[2]]


#### Updated Plot
plot(x, y, type="n")
points(x, y, col="black", pch=20)
abline(a=b0,b=b1,lty=1,lwd=2, col="red")

abline(lm(y~x),lty=2, lwd=2, col="black")
abline(a=a[soln.plain[1]],b=b[soln.plain[2]], lty=1, lwd=2, col="purple") #plain old regression, matches the original black line
abline(a=a[soln.ridge[1]],b=b[soln.ridge[2]], lty=1, lwd=2, col="blue") #ridge regression
abline(a=a[soln.lasso[1]],b=b[soln.lasso[2]], lty=1, lwd=2, col="forestgreen") #lasso regression


## Exercise -->
## Do this for Elastic Net!
## Keep the intercept; only shrink the slope.
## Draw the plots, and interpret. What do you see?
a=seq(-20, 20, 0.05) #start ranges
b=seq(-20, 20, 0.05) #start ranges

object.elasticnet = matrix(NA, nrow=length(a), ncol=length(b))

## Keep the intercept; only shrink the slope
for (i in 1:length(a)) {
  for (j in 1:length(b)) {
    object.elasticnet[i,j] = find.mse(y=y, x=x, a=a[i], b=b[j]) + 0.2*(b[j]^2) + 0.2*(abs(b[j]))
    #object.ridge[i,j] = find.mse(y=y, x=x, a=a[i], b=b[j]) + 0.2*(a[i]^2 + b[j]^2)
    #object.lasso[i,j] = find.mse(y=y, x=x, a=a[i], b=b[j]) + 0.2*(abs(a[i]) + abs(b[j]))
    if((i %% 50==0) & (j==1)) {cat(paste0("iteration: ", i, "\n"))}
  }
}

## Elasticnet
#object.elasticnet
min(object.elasticnet)
soln.elasticnet=which(object.elasticnet == min(object.elasticnet), arr.ind = TRUE)
soln.elasticnet
a[soln.elasticnet[1]] ; b[soln.elasticnet[2]]

plot(x, y, type="n")
points(x, y, col="black", pch=20)
abline(a=b0,b=b1,lty=1,lwd=2, col="red")

abline(lm(y~x),lty=2, lwd=2, col="black")
abline(a=a[soln.elasticnet[1]],b=b[soln.elasticnet[2]], lty=1, lwd=2, col="purple")



###
##
#