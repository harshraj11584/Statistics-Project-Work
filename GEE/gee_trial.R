library(geepack)
data(respiratory)
summary(respiratory)



#trying glm approach
patients <- respiratory
m.glm <- glm(outcome ~ baseline + center + sex + treat + age, data=patients, family = binomial)
summary(m.glm)
#want to find correlation in residuals within a patient manually
m.glm$residuals
patients$residuals <- m.glm$residuals
residuals_all <- matrix( rep( 0, len=111*4), nrow = 111)
#residuals_all
for (i in 1:444)
  residuals_all[floor((i-1)/4 +1),patients$visit[i]] <- patients$residuals[i]
cor(residuals_all)






#trying gee approach

patients <- respiratory
m.ex <- geeglm(outcome ~ baseline + center + sex + treat + age, data=patients, family=binomial, id=interaction(center,id), corstr="exchangeable",std.err="fij")
summary(m.ex)






#GEE Approach for Data from Center 1 from table 2.36 

dataset1 <- read.csv(file='/home/harsh/Desktop/Statistics-Project-Work/GEE/Data2x2_Table2.36 - Center1.csv',header = TRUE)
summary(dataset1)
#View(datase)
model1 <- geeglm(outcome ~ period_effect + treatment_effect + carryover_effect, data=dataset1, id = id, family=binomial, corstr="exchangeable",std.err="fij")
summary(model1)

#want to find correlation in residuals within a patient manually
#model1$residuals
dataset1$residuals <- model1$residuals
residuals_all <- matrix( rep( 0, len=33*2), nrow = 33)
#residuals_all
for (i in c(1:66))
{
  id1  <-  dataset1$id[i] 
  cnum <- dataset1$period_effect[i]+1
  residuals_all[id1,cnum] <- dataset1$residuals[i]
}
cor(residuals_all)






#GEE Approach for Data from Center 2 from table 2.36 

dataset1 <- read.csv(file='/home/harsh/Desktop/Statistics-Project-Work/GEE/Data2x2_Table2.36 - Center2.csv',header = TRUE)
summary(dataset1)
#View(datase)
model1 <- geeglm(outcome ~ period_effect + treatment_effect + carryover_effect, data=dataset1, id = id, family=binomial, corstr="exchangeable",std.err="fij")
summary(model1)

#want to find correlation in residuals within a patient manually
#model1$residuals
dataset1$residuals <- model1$residuals
residuals_all <- matrix( rep( 0, len=67*2), nrow = 67)
#residuals_all
for (i in c(1:134))
{
  id1  <-  dataset1$id[i] 
  cnum <- dataset1$period_effect[i]+1
  residuals_all[id1,cnum] <- dataset1$residuals[i]
}
cor(residuals_all)







#GEE Approach for Combined Data from both Centers from table 2.36 

dataset1 <- read.csv(file='/home/harsh/Desktop/Statistics-Project-Work/GEE/Data2x2_Table2.36 - FullData.csv',header = TRUE)
summary(dataset1)
#View(datase)
model1 <- geeglm(outcome ~ period_effect + treatment_effect + carryover_effect, data=dataset1, id = interaction(center,id), family=binomial, corstr="exchangeable",std.err="fij")
summary(model1)

#want to find correlation in residuals within a patient manually
#model1$residuals
dataset1$residuals <- model1$residuals
residuals_all <- matrix( rep( 0, len=100*2), nrow = 100)
#residuals_all
for (i in c(1:200))
{
  id1  <-  dataset1$id[i] 
  center1 <- dataset1$center[i]
  rnum<-0
  if (center1==1)
    rnum <- as.integer(id1)
  else
    rnum <- as.integer(id1+33)
  cnum <- dataset1$period_effect[i]+1
  residuals_all[rnum,cnum] <- dataset1$residuals[i]
}
cor(residuals_all)
