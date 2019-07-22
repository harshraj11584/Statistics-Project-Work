## Beetles data set
# conc = CS2 concentration
# y= number of beetles killed
# n= number of beetles exposed
# rep = Replicate number (1 or 2)
beetles <- read.table("http://statacumen.com/teach/SC1/SC1_11_beetles.dat", header = TRUE)

library("caTools")
library("car")

beetles
# create data variables: m, y, X
n <- nrow(beetles)
m <- beetles$n
y <- beetles$y
conc <- beetles$conc
conc2<-conc^2

#USING LINEAR LEAST SQUARES FOR FITTING
linear_LS_model <- lm(beetles$y/beetles$n ~ beetles$conc)
linear_LS_model$fitted.values
p_fit1 <- linear_LS_model$fitted.values
plot(beetles$conc,beetles$y/beetles$n)
lines(beetles$conc,p_fit1)

#USING QUADRATIC LEAST SQUARES
quad_LS_model <- lm(beetles$y/beetles$n ~ conc+conc2)
quad_LS_model$fitted.values
p_fit2 <- quad_LS_model$fitted.values
plot(beetles$conc,beetles$y/beetles$n)
lines(beetles$conc,p_fit2)



#Logistic Regression

#Fitting Quadratic Logistic Regression Model

glm_lr <- glm(formula = cbind(y,n-y)~conc+conc2, family=binomial(link="logit"), data=beetles)
summary(glm_lr)
glm_lr$coefficients
glm_lr$fitted.values

p_hat <- glm_lr$fitted.values
p<-beetles$y/beetles$n
plot(beetles$conc,p,legend=TRUE)
points(beetles$conc,p_hat)
lines(conc,p_hat)
p_hat
beetles$conc
