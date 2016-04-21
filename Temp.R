rm(list=ls())

#Determination de l'influence de la temperature : 

data=read.table("temp.txt",h=T)
sunflowerplot(data$T,data$Act,main="Activité de l'enzyme en fonction
     de la température",xlab="Température (°C)",ylab="Activité")
# Température optimale : environ 35 °C