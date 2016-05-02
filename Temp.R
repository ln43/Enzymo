rm(list=ls())

#Determination de l'influence de la temperature : 

data=read.table("temp.txt",h=T)
sunflowerplot(data$T,data$Act,main="Activité de l'enzyme en fonction
     de la température",xlab="Température (°C)",ylab="Activité")
# Température optimale : environ 35 °C
data
#Reste a verifier que la loi d'Arrhenius s'applique
Ea = 8.31*(5-273)*(10-273)/5*log(mean(data$Act[1:3])/mean(data$Act[4:6]))
logQ10 = Ea*10/(8.31*(5-273)*(15-273))
Q10 = exp(logQ10)
