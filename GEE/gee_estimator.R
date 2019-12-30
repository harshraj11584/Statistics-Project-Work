library(geepack)

# WITH CARRYOVER EFFECTS
# GEE Approach for Combined Data from both Centers from table 2.36

dataset1 <- read.csv(file='/home/harsh/Desktop/Statistics-Project-Work/GEE/Data2x2_Table2.36 - FullData.csv',header = TRUE)
summary(dataset1)
#View(dataset1)
model1 <- geeglm(outcome ~ period_effect + treatment_effect + carryover_effect, data=dataset1, id = id, family=binomial, corstr="exchangeable",std.err="fij")
summary(model1)


# WITHOUT CARRYOVER EFFECTS
# GEE Approach for Combined Data from both Centers from table 2.36

dataset2 <- read.csv(file='/home/harsh/Desktop/Statistics-Project-Work/GEE/Data2x2_Table2.36 - FullData.csv',header = TRUE)
dataset2 <- subset(dataset2,select=-c(carryover_effect))
summary(dataset2)
#View(dataset2)
model2 <- geeglm(outcome ~ period_effect + treatment_effect, data=dataset2, id = id, family=binomial, corstr="exchangeable",std.err="fij")
summary(model2)

