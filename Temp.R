rm(list=ls())
library(nlstools)

# Determination de l'influence de la temperature :

data=read.table("temp.txt",h=T)
Act = (data$Act*4*10^(-3))/(15*13600)*10^6
# Representation graphique
sunflowerplot(data$T,Act,main="Activité de l'enzyme en fonction
     de la température",xlab="Température (°C)",ylab="Activité (µmol/min)",
              cex.main=1.2, cex.lab=1, cex.sub=1)

# Loi d'Arrhenius
v = Act[1:21]*10^(-6)/(60*4*10^(-3))
k = v/(7*10^(-4))
invT = 1/(8.31*(data$T[1:21]+273))
plot(invT,log(k),main="ln(k)=f(1/RT)",xlab="1/RT (T en K)",ylab="log(k)")
lm1 = lm(log(k)~invT)
par(mfrow=c(2,2))
plot(lm1)
summary(lm1)
abline(lm1,col="red")
shapiro.test(residuals(lm1))


logQ10 = 17110*10/(8.31*(5+273)*(15+273))
Q10 = exp(logQ10)
