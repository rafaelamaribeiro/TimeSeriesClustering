rm(list=ls())

##############Packages#######################

#install.packages("TSdist")
library(ggplot2)
library(lubridate)
library(stats)
library(timeDate)
library(clue)
library(fpc)
library(clv)
library(factoextra)
library(cluster)
#library(TSclust)
library(TSclust)
library(forecast)

#############Opening the data###############

#setwd("C:\\Users\\Rafaela\\Downloads\\2022\\TFC\\codigos eol")
#setwd("/home/rafaela/Documents/PUC/2022.1/TFC/codigos/Fernando")
#setwd("C:\\Users\\Rafaela\\Downloads")
getwd()

fernando <- read.csv2("hweol_FERNANDO_RN.csv")
head(fernando)

#fernando$dates_new <- as.POSIXct(fernando$dates, format = "%Y-%m-%d %H:%M")
fernando$dates_new <- as.Date(fernando$dates, format = "%Y-%m-%d %H:%M")

maximo = max(fernando$FERNANDO_RN)
minimo = min(fernando$FERNANDO_RN)
#lim <- ceiling(maximo)

fernando$months <- month(fernando$dates_new)
fernando$months1 <- as.character(fernando$months)

head(fernando)

#PlotFernandoBoxplot <- ggplot(fernando, 
#                              aes(y = FERNANDO_RN, x = reorder(months1, months))) + 
#  geom_boxplot(outlier.colour="black", outlier.shape=16,
#               outlier.size=2) +
#  xlab("Months") + ylab("Wind Speed (m/s)") + 
#  theme(legend.position = "none") + theme_classic()
#PlotFernandoBoxplot

#for (i in 1:length(fernando$FERNANDO_RN))
#{fernando$wind_speed_norm[i] = (fernando$FERNANDO_RN[i]-minimo)/(maximo-minimo)}

#PlotFernando <- ggplot(fernando) + 
#  geom_line(aes(x=dates_new, y=FERNANDO_RN)) +
#  xlab("Hours") + ylab("Wind Speed (m/s)") + 
#  theme(legend.position = "none") + theme_classic() + ylim(c(0,lim))
#PlotFernando

max(fernando$wind_speed_norm)

############Data Preparation#######################

#vec_norm = as.vector(fernando$wind_speed_norm)
vec = as.vector(fernando$FERNANDO_RN)
length(vec)

CreatePsi <- function (vec) {
  Nh = 24
  Ni = (365*41)
  Nc = 1
  matriz = matrix(0, nrow =  Ni, ncol = Nc * Nh)
  m = 1
  a = 1
  for (i in 1:length(vec)) {
    matriz[a, m] = vec[i]
    m = m + 1
    if (m > Nh) {
      m = 1
      a = a + 1
    }
  }
  return(matriz)
}

Psi = CreatePsi(vec)
Psi
#Psi_norm = CreatePsi(vec_norm)
#Psi_norm
head(Psi)
tail(Psi)

cols = ncol(Psi)
rows = nrow(Psi)

d <- diss(Psi, METHOD = "AR.LPC.CEPS", k = 10, order = NULL)






v <- matrix(ncol = 3, nrow = rows)

for (i in 1:rows) {
  aa <- auto.arima(Psi[i,], method = "ML");
  v[i,] <- aa$arma[1:3]}

m <- matrix(ncol = rows, nrow = rows)
for (r1 in 1:rows) {
  for (r2 in 1:rows) {
    m[r1,r2] <- diss.AR.PIC(Psi[r1,], Psi[r2,], order.x = v[r1,],
                            order.y = v[r2,]);
    m[r2,r1] <- m[r1,r2]
  }
}

m <- matrix(ncol = rows, nrow = rows)
for (r1 in 1:rows) {
  for (r2 in 1:rows) {
    m[r1,r2] <- diss.AR.PIC(Psi[r1,], Psi[r2,], order.x = NULL,
                            order.y = NULL);
    m[r2,r1] <- m[r1,r2]
  }
}






start_time <- Sys.time()
d <- diss(Psi, METHOD = "AR.PIC", order = v)
end_time <- Sys.time()
time <- end_time - start_time

###############ARIMA###############

#diss.AR.PIC(Psi[1,], Psi[2,])
#start_time <- Sys.time()
#m <- matrix(ncol = rows, nrow = rows)
#for (r1 in 1:rows) {
#  for (r2 in 1:rows) {
#    m[r1,r2] <- diss.AR.PIC(Psi[r1,], Psi[r2,], order.x = v[r1,],
#                            order.y = v[r2,]);
#    m[r2,r1] <- m[r1,r2]
#  }
#}
#end_time <- Sys.time()

#a <- matrix(rep(c(1,0,1), rows), ncol = 3, byrow = T)
#a

#plot(Psi[1,])
#acf(Psi[1,])
#pacf(Psi[1,])
#fit = arima(Psi[1,], order=c(1,0,0),
#            seasonal = list(order = c(1, 0, 0), period = NA))

#fit = auto.arima(Psi[48,])
#fit$arma[1:3]

#install.packages("forecast")
#library(forecast)
#v <- matrix(ncol = 3, nrow = 50)
#dPsi <- diff(Psi)


#aa <- auto.arima(Psi[500,])
#aa$model$Z
#fit = arima(Psi[1,], order=c(2,0,0))

#diff(Psi)
#plot(diff(Psi[1,], lag = 12))

#diss.AR.PIC(Psi[49,], Psi[1,], order.x = v[49,], order.y = v[1,])

for (i in 1:rows) {
  aa <- auto.arima(Psi[i,], method = "ML")
  if (aa$arma[1] == 0 && aa$arma[3] == 0) {
    v[i,] <- c(1,0,0)}
  else {v[i,] <- aa$arma[1:3]}}


#d = diss(Psi_norm[1:50,], METHOD="AR.PIC", order = v[1:50,])
#d = diss(Psi[1:50,], METHOD = "AR.PIC", order = NULL)
#d = diss(Psi[1:50,], METHOD = "AR.MAH") # Demora mais
#pvalues.clust(d$p_value, significance = 0.1) # escolhe o num de cluster automaticamente
start_time <- Sys.time()
d <- diss(Psi, METHOD = "AR.PIC", order = NULL)
end_time <- Sys.time()
time <- end_time - start_time

###############Silhouette, CH, DB pam##################

n_rows <- 12
mat = matrix(0, nrow = n_rows)
mat2 = matrix(0, nrow = n_rows)
mat3 = matrix(0, nrow = n_rows)
for (i in 2:n_rows){
  set.seed(321)
  clust = pam(d, i, diss = TRUE)
  sil = silhouette(clust$cluster, d)
  mat[i] = mean(as.matrix(sil)[,3])
  CH = calinhara(Psi_norm, clust$cluster)
  mat2[i] = as.matrix(CH)
  obj <- cls.scatt.data(Psi_norm, clust$cluster, dist = "euclidean")
  mat3[i] <- clv.Davies.Bouldin(obj, c("average"), c("average"))
}
colnames(mat) <- c("Avg_Silhouette_Value")
colnames(mat2) <- c("CH_Value")
max(mat2)
colnames(mat3) <- c("DB_Value")
min(mat3[2:12])

plot(mat)



# teste
c<-pam(d, 2, diss = TRUE)
summary(c)
sil <- silhouette(c$cluster, d)
mean(as.matrix(sil)[,3])
ch <- calinhara(Psi[1:100,], c$cluster)
obj <- cls.scatt.data(Psi[1:100,], c$cluster, dist = "euclidean")
db <- clv.Davies.Bouldin(obj, c("average"), c("average"))


