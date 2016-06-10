rm(list=ls())
library(nlstools)


#===============================================================================================
#============================== Regression Linéaire ============================================
#===============================================================================================

rm(list=ls())
library(nlstools)

# Import des données 
d=read.table('EnzymoReg.txt', header=T)

S=d$"S"
V=d$"V"

# Visualtisation des données
plot(S[1:10],V[1:10],main="Vitesse de l'activité enzymatique de l'éserine \n en fonction de la concentration en ATC", xlab="[ATC] en M",ylab="V en µmol/min.mg",col="black",pch=20,ylim=c(0,500))
points(S[11:20],V[11:20], col= "red",pch=20)
points(S[21:30],V[21:30], col= "blue",pch=20)
points(S[31:40],V[31:40], col= "green",pch=20)
legend("bottomright",col = c("black","red","blue", "green"),legend=c("pas d'ésérine","[ésérine] = 2","[ésérine] = 4","[ésérine] = 8"),pch=20)

#Le modèle
michaelisI<-as.formula("V ~ ((Vm * S)/((1+ ( I / Ki))*( Km + S)))")
valinit<-list(Km=0.000124,Ki=2.145342e-07,Vm=553.4692)
mod=nls(formula=michaelisI,data=d,start=valinit)

plotfit(mod,main="Vitesse de l'activité enzymatique de l'éserine \n en fonction de la concentration en ATC", xlab="[ATC] en M",ylab="V en µmol/min.mg")
overview(mod)

#Vérification des hypothèses du modèle

#Residus
resbar=nlsResiduals(mod)
plot(resbar)

qqnorm(residuals(mod))
qqline(residuals(mod))

shapiro.test(residuals(mod))

#Regions de confiance
beale=nlsConfRegions(mod)
plot(beale)

#Contour de la Somme des Carres des ecarts residuels 
cont=nlsContourRSS(mod,lseq=50)
plot(cont,nlev=10)
