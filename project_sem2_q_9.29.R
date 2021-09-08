install.packages('psych')
library(psych)
library(ggplot2)
library("car")

my_theme <- theme(panel.background = element_rect(fill = "white")) +
  theme(plot.title = element_text(hjust = 0.5 , face = "bold", color = "cadetblue" )) +
  theme(panel.grid = element_blank())+
  theme(plot.subtitle = element_text(hjust = 0.5 , face = "bold", color = "magenta" ))+
  theme(axis.title.x = element_text(hjust = 0.5 , face = "bold", color = "yellow" ))+
  theme(axis.title.y = element_text(hjust = 0.5 , face = "bold", color = "yellow" ))+
  theme(axis.text.x = element_text(hjust = 0.5 , face = "bold", color = "green" ))+
  theme(axis.text.y = element_text(hjust = 0.5 , face = "bold", color = "green" ))+
  theme(plot.background = element_rect(fill = "blue")) +
  theme(legend.position = 'bottom',legend.background = element_rect(fill = 'yellow'))

# Getting the data ready

trackdata = read.csv('data.csv')
trackdata
sp100 <- 100/trackdata$X100m
sp200 <- 200/trackdata$X200m
sp400 <- 400/trackdata$X400m
sp800 <- 800/(trackdata$X800m*60)
sp1500 <- 1500/(trackdata$X1500m*60)
sp5000 <- 5000/(trackdata$X5000m*60)
sp10000 <- 10000/(trackdata$X10000m*60)
spmar <- 42195/(trackdata$Marathon*60)
spmar
dataspeed <- data.frame(sp100,sp200,sp400,sp800,sp1500,sp5000,sp10000,spmar)
datamat <- as.matrix(dataspeed); datamat



n <- nrow(datamat); n
p <- ncol(datamat); p
summary(datamat)
dataMean <- colMeans(datamat); dataMean
sqvar <- sqrt((diag(var(datamat))))
datavar <- var(datamat); datavar
datacor <- cor(datamat); datacor

centered_data<- datamat - matrix(rep(dataMean,n),n,p,byrow=TRUE)
stddata<-centered_data/matrix(rep(sqvar,n),n,p,byrow=TRUE)

ggplot(trackdata, aes(x=sp100, y=spmar)) + geom_point() #ggplot scatterplot

scatterplotMatrix(datamat)

ggplot(trackdata,aes(y=reorder(Country,-X100m),x=X100m))+geom_col()+ylab("Countries")+ggtitle("Speed of Various Countries in 100m race")
ggplot(trackdata,aes(y=reorder(Country,-X200m),x=X200m))+geom_col()+ylab("Countries")+ggtitle("Speed of Various Countries in 200m race")
ggplot(trackdata,aes(y=reorder(Country,-X400m),x=X400m))+geom_col()+ylab("Countries")+ggtitle("Speed of Various Countries in 400m race")
ggplot(trackdata,aes(y=reorder(Country,-X800m),x=X800m))+geom_col()+ylab("Countries")+ggtitle("Speed of Various Countries in 800m race")
ggplot(trackdata,aes(y=reorder(Country,-Marathon),x=Marathon))+geom_col()+ylab("Countries")+ggtitle("Speed of Various Countries in Marathon race")

corPlot(datacor);title(main="Correlation Plot of National Track Records for Men",sub="Source: IAAF/ATES Track and Field Statistics Handbook for the Helinski 2005 Olympics.Courtesy of Ottavio Castellini",col.main="blue",col.sub="red")

boxplot(datamat,notch=TRUE,names=c("100m","200m","400m","800m","1500m","5000m","10000m","Marathon(42195m)"));title(main="Boxplot of National Track Records for Men",sub="Source: IAAF/ATES Track and Field Statistics Handbook for the Helinski 2005 Olympics.Courtesy of Ottavio Castellini",xlab="Length of Track (in m)",ylab="Speed (in m/s)",col.main="blue",col.sub="red",col.lab="magenta")

# =======================================================
# Factor Analysis by PCA                                
# =======================================================

#=========================================================================
# Principal Component Analysis by "prcomp" method on S (covariance matrix)
#=========================================================================
#####
#pca

dataPCAS <- prcomp(datamat,center = FALSE,scale. = FALSE); dataPCAS[[2]][,1:2]
summary(dataPCAS)
plot(dataPCAS,type="lines", main="Screeplot")
###biplot(dataPCAS)
dataPCAS
# factor loading
dataFactLoadS <- dataPCAS[[2]][,1:2]*t(matrix(rep(dataPCAS[[1]][1:2],p),2,p))  #[\sqrt(\lambda_1)e_1:...:\sqrt(\lambda_p)e_p]
print(dataFactLoadS)  #Loadings on factors, i.e., full L matrix

# rotated factor loading
varimax(dataFactLoadS)$loadings
# communalities
commS <- cbind(dataFactLoadS,varimax(dataFactLoadS)$loadings,'communalities'=dataFactLoadS[,1]^2 + dataFactLoadS[,2]^2)
commS
# specific variance
specvarS <- diag(stddata) - commS[,5]; specvarS
#####
#screeplot
eigenvalues<-eigen(datavar)[[1]]
plot(eigenvalues,type="lines");title(main="Screeplot",sub="Index refers to the Factors")

pca_S <- prcomp(datamat)[[2]]; pca_S
sd_S <- prcomp(datamat)[[1]]; sd_S
fac_S <- pca_S * t(matrix(rep(sd_S,p),p,p)) ; fac_S # factor loading
varimax(fac_S)$loadings # rotated factor loading
(varimax(fac_S)$loadings)*10
#table
tableS1 <- cbind(fac_S[,1:2],varimax(fac_S)$loadings[,1:2],'communalities'=fac_S[,1]^2 + fac_S[,2]^2,specvarS[,2])
colnames(tableS1) <- c('F1','F2','rotated F1','rotated F2','communalities','specific variance')
tableS1
specvarS <- diag(datavar) - t(apply(fac_S^2,1,cumsum)) #specific variances
t_var_S <- (cumsum(sd_S^2))/sum(sd_S^2) ; t_var_S  #total explained variability
exp_var_S<-sd_S^2/sum(sd_S^2);exp_var_S


# factor scores
fact_scoreS <- t(fac_S)%*%solve(datavar)%*%t(centered_data)
print(fact_scoreS)
plot(fact_scoreS[1,],fact_scoreS[2,],xlab = 'Factor 1',ylab='Factor 2',main='Factor scores')

# outliers
outliers <- trackdata$Country[fact_scoreS[1,]< (-2)]
outliers



#==========================================================================
# Principal Component Analysis by "prcomp" method on R (correlation matrix)
#==========================================================================
#screeplot
eigenvalues<-eigen(datacor)[[1]]
plot(eigenvalues,type="lines");title(main="Screeplot",sub="Index refers to the Factors")


# pca
dataPCAR <- prcomp(datamat,center = TRUE,scale. = TRUE); dataPCAR[[2]][,1:2]
summary(dataPCAR)

dataPCAR1<-list(dataPCAR[[1]],dataPCAR[[2]]);dataPCAR1

# factor loading
dataFactLoadR <- dataPCAR[[2]][,1:2]*t(matrix(rep(dataPCAR[[1]][1:2],p),2,p))  #[\sqrt(\lambda_1)e_1:...:\sqrt(\lambda_p)e_p]
print(dataFactLoadR)  #Loadings on factors,i.e., full L matrix

dataSpecVarR<-1-t(apply(dataFactLoadR^2,1,cumsum));dataSpecVarR
# rotated factor loading
varimax(dataFactLoadR)$loadings
# communalities
commR <- cbind(dataFactLoadR,varimax(dataFactLoadR)$loadings,'communalities'=dataFactLoadR[,1]^2 + dataFactLoadR[,2]^2)
commR

# specific variance
specvarR <- 1 - commR[,5]; specvarR
tableR <- cbind(commR,specvarR); tableR
colnames(tableR) <- c('F1','F2','rotated F1','rotated F2','communalities','specific variance')
tableR


fac <- fa(datacor,nfactors = 2,rotate = 'varimax',fm='mle')
fa.diagram(fac)

# Factor Score by Regression L_z'R^{-1}z'
fact_scoreR <- t(dataFactLoadR)%*%solve(datacor)%*%t(stddata)
fact_scoreR
plot(fact_scoreR[1,],fact_scoreR[2,],xlab = 'Factor 1',ylab='Factor 2',main='Factor scores')

# outliers
outliers <- trackdata$Country[fact_scoreR[1,] < (-2)]
outliers


# residual matrix on 1 factor
Res1=datacor-dataFactLoadR[,1]%*%t(dataFactLoadR[,1])-diag(dataSpecVarR[,1])
#Res1=datacor-dataFactLoadR[,1]%*%t(dataFactLoadR[,1])-diag(specvarR)
print(Res1) #R - LL' -\Psi, when one factor is considered
sum(Res1*Res1)

# residual matrix on 2 factors
Res2=datacor-dataFactLoadR[,c(1,2)]%*%t(dataFactLoadR[,c(1,2)])-diag(specvarR)
print(Res2) #R - LL' -\Psi, when two factors are considered
sum(Res2*Res2)





#===================#
#       MLE         #
#===================#

factmle <- factanal(datamat,factors = 2,scores = "none")
print(factmle)
mlespvr <- 1-t(apply(factmle$loadings^2,1,cumsum))
print(diag(mlespvr[,2]))
# R-LL'-\shi
Res2M<-datacor-factmle$loadings[,c(1,2)]%*%t(factmle$loadings[,c(1,2)])-diag(mlespvr[,2])
print(Res2M)
sum(Res2M*Res2M)


