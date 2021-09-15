#On importe les packages nécessaires

install.packages("zoo")
install.packages("tseries")
library(zoo)
library(tseries)
require(zoo)
require(tseries)
library(fUnitRoots)
require(forecast)

#On importe les données
#lien des données = " https://www.insee.fr/fr/statistiques/serie/001562489#Tableau "
data <- read.csv("/Users/Rayane/Documents/valeurs_mensuelles.csv",sep=";")

View(data)

#On prépare la base de données afin de travailler dessus
dates_char <- as.character(data$dates)
dates_char[1];dates_char[length(dates_char)]
dates<-as.yearmon(seq(from=1990+1/12,to=2012+12/12,by=1/12))
savon <- zoo(data$savon,order.by=dates)
plot(savon)

#Graphiquement elle ne semble pas stationnaire. Elle semble comporter une tendance linéaire croissante
#On vérifie ca avec une régression linéaire de savon sur le temps

lt<-lm(savon~dates)

r<-lt$residuals
summary(lt)

#On différencie une fois pour essayer de rendre notre STL stationnaire
dsavon <- diff(savon,1)
plot(cbind(savon,dsavon))

#On veut vérifier la stationnarité, donc on effectue différents tests

#On réalise tout d'abord un pp-test

pp.test(dsavon)

#On réalise ensuite un ADF test, nc car pas de constante ni de trend


adf<-adfTest(dsavon,lag=1,type="nc")
adf


#On vérifie que les résidus du modèle d'autoregressions sont bien non autocorrélés
Qtests <- function(series, k, fitdf=0){
  pvals <- apply(matrix(1:k), 1, FUN=function(l){
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))})
  return(t(pvals))}

t(Qtests(adf@test$lm$residuals[1:24],24))

#On réalise un kpss-test
kpss.test(dsavon)


#On affiche les autocorrélations et les autocorrélations partielles
par(mfrow=c(1,2))
acf(dsavon);pacf(dsavon)

#On calibre notre modèle ARMA

#affichage du test de significativité des coefficients
signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

#On estime des modèles arima et on en vérifie l'ajustement et la validité

modelchoice <- function(p,q,data=dsavon, k=24){
  estim <- try(arima(data, c(p,0,q),optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
  
  arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05
  masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
  resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
  checks <- c(arsignif,masignif,resnocorr)
  print(checks)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
}

#On teste tout les arma pour trouver les bonnes valeurs de p et q
armamodelchoice <- function(pmax,qmax){
  pqs <- expand.grid(0:pmax,0:qmax)
  t(apply(matrix(1:dim(pqs)[1]),1,function(row) {
    p <- pqs[row,1]; q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",",q,") \n"))
    modelchoice(p,q)
  }))
}
#On a fixé pmax et qmax grâce à l'acf et le pacf
pmax=11; qmax=13
T <- length(dsavon)


# On regarde les modèles bien ajustés et valides
armamodels <- armamodelchoice(pmax,qmax) #estime tous les arima
selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),] #modèles bien ajustés et valides
selec

#On regarde tout les p et q candidats pour le modèlet.
pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2]))) 
names(pqs) <- paste0("arma(",selec[,1],",",selec[,2],")")
models <- lapply(pqs, function(pq) arima(r,c(pq[["p"]],0,pq[["q"]]))) #crée une liste des modèles
vapply(models, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m)))

#On détermine le modèle que l'on conservera en minimisant le R2 de nos modèles
adj_r2 <- function(model){
  ssres <- sum(model$residuals^2) #somme des résidus au carré
  p <- model$arma[1]
  q <- model$arma[2]
  sstot <- sum(dsavon[-c(1:max(p,q))]^2) #somme des observations de l' ??echantillon au carré
  n <- model$nobs-max(p,q) #taille de l'échantillon
  adj_r2 <- 1-(ssres/(n-p-q-1))/(sstot/(n-1)) #r2 ajusté
  return(adj_r2)}
adj_r2(arma603)
adj_r2(arma0013)

#Fonction pour implémenter les modèles arima
arimafit <- function(estim){
  adjust <- round(signif(estim),3)
  pvals <- Qtests(estim$residuals,36,length(estim$coef)-1)
  pvals <- matrix(apply(matrix(1:36,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)
  colnames(pvals) <- rep(c("lag", "pval"),6)
  cat("tests de nullité des coefficients :\n")
  print(adjust)
  cat("\n tests d'absence d'autocorrélation des résidus : \n")
  print(pvals)
}

#On implémente notre modèle et arima

arma603 <- arima(dsavon,c(6,0,3)); arimafit(arma603)
arma0013<- arima(dsavon,c(0,0,13));arimafit(arma0013)

arma603
arma0013


#On trace la représentation de la prévision à 95%
autoplot(forecast(arma603, level = c(95)))
autoplot(forecast(arma0013,level=c(95)))

arima613 <- arima(savon,c(6,1,3),method='ML')
arima0113 <-arima(savon,c(0,1,13),method='ML')

#On teste l'absence d'autocorrélation des résidus.
t(Qtests(arima613$residuals,24))
Qtests(arima0113$residuals,24)

autoplot(forecast(arima613),level=c(95))
autoplot(forecast(arima0113),level=c(95))



