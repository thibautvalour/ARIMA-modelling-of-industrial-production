#On importe les packages n�cessaires

install.packages("zoo")
install.packages("tseries")
library(zoo)
library(tseries)
require(zoo)
require(tseries)
library(fUnitRoots)
require(forecast)

#On importe les donn�es
#lien des donn�es = " https://www.insee.fr/fr/statistiques/serie/001562489#Tableau "
data <- read.csv("/Users/Rayane/Documents/valeurs_mensuelles.csv",sep=";")

View(data)

#On pr�pare la base de donn�es afin de travailler dessus
dates_char <- as.character(data$dates)
dates_char[1];dates_char[length(dates_char)]
dates<-as.yearmon(seq(from=1990+1/12,to=2012+12/12,by=1/12))
savon <- zoo(data$savon,order.by=dates)
plot(savon)

#Graphiquement elle ne semble pas stationnaire. Elle semble comporter une tendance lin�aire croissante
#On v�rifie ca avec une r�gression lin�aire de savon sur le temps

lt<-lm(savon~dates)

r<-lt$residuals
summary(lt)

#On diff�rencie une fois pour essayer de rendre notre STL stationnaire
dsavon <- diff(savon,1)
plot(cbind(savon,dsavon))

#On veut v�rifier la stationnarit�, donc on effectue diff�rents tests

#On r�alise tout d'abord un pp-test

pp.test(dsavon)

#On r�alise ensuite un ADF test, nc car pas de constante ni de trend


adf<-adfTest(dsavon,lag=1,type="nc")
adf


#On v�rifie que les r�sidus du mod�le d'autoregressions sont bien non autocorr�l�s
Qtests <- function(series, k, fitdf=0){
  pvals <- apply(matrix(1:k), 1, FUN=function(l){
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))})
  return(t(pvals))}

t(Qtests(adf@test$lm$residuals[1:24],24))

#On r�alise un kpss-test
kpss.test(dsavon)


#On affiche les autocorr�lations et les autocorr�lations partielles
par(mfrow=c(1,2))
acf(dsavon);pacf(dsavon)

#On calibre notre mod�le ARMA

#affichage du test de significativit� des coefficients
signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

#On estime des mod�les arima et on en v�rifie l'ajustement et la validit�

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
#On a fix� pmax et qmax gr�ce � l'acf et le pacf
pmax=11; qmax=13
T <- length(dsavon)


# On regarde les mod�les bien ajust�s et valides
armamodels <- armamodelchoice(pmax,qmax) #estime tous les arima
selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),] #mod�les bien ajust�s et valides
selec

#On regarde tout les p et q candidats pour le mod�let.
pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2]))) 
names(pqs) <- paste0("arma(",selec[,1],",",selec[,2],")")
models <- lapply(pqs, function(pq) arima(r,c(pq[["p"]],0,pq[["q"]]))) #cr�e une liste des mod�les
vapply(models, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m)))

#On d�termine le mod�le que l'on conservera en minimisant le R2 de nos mod�les
adj_r2 <- function(model){
  ssres <- sum(model$residuals^2) #somme des r�sidus au carr�
  p <- model$arma[1]
  q <- model$arma[2]
  sstot <- sum(dsavon[-c(1:max(p,q))]^2) #somme des observations de l' ??echantillon au carr�
  n <- model$nobs-max(p,q) #taille de l'�chantillon
  adj_r2 <- 1-(ssres/(n-p-q-1))/(sstot/(n-1)) #r2 ajust�
  return(adj_r2)}
adj_r2(arma603)
adj_r2(arma0013)

#Fonction pour impl�menter les mod�les arima
arimafit <- function(estim){
  adjust <- round(signif(estim),3)
  pvals <- Qtests(estim$residuals,36,length(estim$coef)-1)
  pvals <- matrix(apply(matrix(1:36,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)
  colnames(pvals) <- rep(c("lag", "pval"),6)
  cat("tests de nullit� des coefficients :\n")
  print(adjust)
  cat("\n tests d'absence d'autocorr�lation des r�sidus : \n")
  print(pvals)
}

#On impl�mente notre mod�le et arima

arma603 <- arima(dsavon,c(6,0,3)); arimafit(arma603)
arma0013<- arima(dsavon,c(0,0,13));arimafit(arma0013)

arma603
arma0013


#On trace la repr�sentation de la pr�vision � 95%
autoplot(forecast(arma603, level = c(95)))
autoplot(forecast(arma0013,level=c(95)))

arima613 <- arima(savon,c(6,1,3),method='ML')
arima0113 <-arima(savon,c(0,1,13),method='ML')

#On teste l'absence d'autocorr�lation des r�sidus.
t(Qtests(arima613$residuals,24))
Qtests(arima0113$residuals,24)

autoplot(forecast(arima613),level=c(95))
autoplot(forecast(arima0113),level=c(95))



