## Leukemia white blood cell types example
# ntotal = number of patients with IAG and WBC combination
# nres = number surviving at least one year
# ag = 1 for AG+, 0 for AG-
# wbc  = white cell blood count
# lwbc = log white cell blood count
# p.emp = Emperical Probability
leuk <- read.table("http://statacumen.com/teach/SC1/SC1_11_leuk.dat", header = TRUE)
leuk$lwbc <- log(leuk$wbc)
leuk$p.emp <- leuk$nres / leuk$ntotal
leuk$ag  


glm_lr = glm(formula = cbind(leuk$nres,leuk$ntotal-leuk$nres) ~ leuk$ag+leuk$lwbc, family = binomial, data = leuk)  
summary(glm_lr)

leuk$p.MLE <- glm_lr$fitted.values

plot(leuk[(leuk$ag==1),]$lwbc, leuk[(leuk$ag==1),]$p.emp)
points(leuk[(leuk$ag==0),]$lwbc, leuk[(leuk$ag==0),]$p.emp)
points(leuk[(leuk$ag==1),]$lwbc, leuk[(leuk$ag==1),]$p.MLE,col='red')
points(leuk[(leuk$ag==0),]$lwbc, leuk[(leuk$ag==0),]$p.MLE,col='red')
