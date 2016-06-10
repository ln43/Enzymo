rm(list=ls())

library(nlstools)

#Détermination des paramètres cinétiques avec 3 substrats de l'ACE : 
# Modèle de Michaelis-Menten


#ATC
#Representations et determination des params
d1=read.table("ATC.txt",h=T)
d1$"S"=d1$S*10^(-4)
sunflowerplot(d1, xlab = "[S]", main = "V = f([S])") # V = f[S]
invv = 1/d1$v
invs = 1/d1$S
sunflowerplot(invs,invv,xlab="1/[S]",ylab="1/V",main="1/V = f(1/[S])") # Lineweaver et Burk
summary(lm(invv~invs))
abline(lm(invv~invs),col=2)
sv = d1$S/d1$v
sunflowerplot(d1$S,sv,xlab="[S]",ylab="[S]/V",main="[S]/V = f([S])") # Hanes
summary(lm(sv~d1$"S"))
abline(lm(sv~d1$"S"),col=2)

mod1<-michaelis
vi1<-list(Vmax=600,Km=1e-04)
aj1<-nls(mod1,d1,vi1)
plotfit(aj1,smooth=T) 
overview(aj1)

#Regions de confiance
beale=nlsConfRegions(aj1)
plot(beale)

#residus
resbar=nlsResiduals(aj1)
plot(resbar)

shapiro.test(residuals(aj1))#p-value <0.05!

#Contour de la Somme des Carres des ecarts residuels 
# (rouge : region de confiance selon Beale)
contour1=nlsContourRSS(aj1,lseq=50)
plot(contour1,nlev=10)

##############################

#PTC
d2=read.table("PTC.txt",h=T)
d2$"S"=d2$S*10^(-4)
sunflowerplot(d2, xlab = "[S]", main = "V = f([S])") # V = f[S]
invv = 1/d2$v
invs = 1/d2$S
sunflowerplot(invs,invv,xlab="1/[S]",ylab="1/V",main="1/V = f(1/[S])") # Lineweaver et Burk
summary(lm(invv~invs))
abline(lm(invv~invs),col=2)
sv = d2$S/d2$v
sunflowerplot(d2$S,sv,xlab="[S]",ylab="[S]/V",main="[S]/V = f([S])") # Hanes
summary(lm(sv~d2$"S"))
abline(lm(sv~d2$"S"),col=2)


mod2<-michaelis
vi2<-list(Vmax=300,Km=1e-04)
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

#Contour de la Somme des Carres des ecarts residuels 
# (rouge : region de confiance selon Beale)
contour2=nlsContourRSS(aj2,lseq=50)
plot(contour2,nlev=10)

##############################

#BTC
d3=read.table("BTC.txt",h=T)
d3$"S"=d3$S*10^(-4)
sunflowerplot(d3, xlab = "[S]", main = "V = f([S])") # V = f[S]
invv = 1/d3$v
invs = 1/d3$S
sunflowerplot(invs,invv,xlab="1/[S]",ylab="1/V",main="1/V = f(1/[S])") # Lineweaver et Burk
summary(lm(invv~invs))
abline(lm(invv~invs),col=2)
sv = d3$S/d3$v
sunflowerplot(d3$S,sv,xlab="[S]",ylab="[S]/V",main="[S]/V = f([S])") # Hanes
summary(lm(sv~d3$"S"))
abline(lm(sv~d3$"S"),col=2)
