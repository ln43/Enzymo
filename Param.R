rm(list=ls())

#library(nlstools)

#Détermination des paramètres cinétiques avec 3 substrats de l'ACE : 
# Modèle de Michaelis-Menten


#ATC
d1=read.table("ATC.txt",h=T)
d1$"S"=d1$S*10^(-4)
d1=d1[1:21,] #Probable erreur sur les deux dernières concentrations, à vérifier.
plot(d1)

mod1<-michaelis
vi1<-list(Vmax=1.5,Km=1e-04)
aj1<-nls(mod1,d1,vi1)
plotfit(aj1,smooth=T) 
overview(aj1)

#Regions de confiance
beale=nlsConfRegions(aj1)
plot(beale)

#residus
resbar=nlsResiduals(aj1)
plot(resbar) 

shapiro.test(residuals(aj1))

#Contour de la Somme des Carrés des écarts résiduels (rouge : région de confiance selon Beale)
contour1=nlsContourRSS(aj1,lseq=50)
plot(contour1,nlev=10)

##############################

#PTC
d2=read.table("PTC.txt",h=T)
d2$"S"=d2$S*10^(-4)
plot(d2)

mod2<-michaelis
vi2<-list(Vmax=1.5,Km=1e-04)
aj2<-nls(mod2,d2,vi2)
plotfit(aj2,smooth=T) 
overview(aj2)

#Regions de confiance
beale=nlsConfRegions(aj2)
plot(beale)

#residus
resbar=nlsResiduals(aj2)
plot(resbar) 

shapiro.test(residuals(aj2))

#Contour de la Somme des Carrés des écarts résiduels (rouge : région de confiance selon Beale)
contour2=nlsContourRSS(aj2,lseq=50)
plot(contour2,nlev=10)

##############################

#BTC --> ???
d3=read.table("BTC.txt",h=T)
d3$"S"=d3$S*10^(-4)
plot(d3)

mod3<-michaelis
vi3<-list(Vmax=1.5,Km=1e-04)
aj3<-nls(mod3,d3,vi3)
plotfit(aj3,smooth=T) 
overview(aj3)

#Regions de confiance
beale=nlsConfRegions(aj3)
plot(beale)

#residus
resbar=nlsResiduals(aj3)
plot(resbar) 

shapiro.test(residuals(aj3))

#Contour de la Somme des Carrés des écarts résiduels (rouge : région de confiance selon Beale)
contour3=nlsContourRSS(aj3,lseq=50)
plot(contour3,nlev=10)