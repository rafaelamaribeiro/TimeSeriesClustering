rm(list=ls())

##############Packages#######################

#install.packages("dendextend")
library(fpc)
library(ggplot2)
library(lubridate)
library(stats)
library(timeDate)
library(clue)
library(cluster)
library(factoextra)
library(clv)
library(dendextend)

#############Opening the data###############

#setwd("C:\\Users\\Rafaela\\Downloads\\2022\\TFC\\codigos eol")
#setwd("/home/rafaela/Documents/PUC/2022.1/TFC/codigos")
getwd()

fernando <- read.csv2("hweol_FERNANDO_RN.csv")
head(fernando)

#fernando$dates_new <- as.POSIXct(fernando$dates, format = "%Y-%m-%d %H:%M")
fernando$dates_new <- as.Date(fernando$dates, format = "%Y-%m-%d %H:%M")

#maximo = max(fernando$FERNANDO_RN)
#minimo = min(fernando$FERNANDO_RN)
#lim <- ceiling(maximo)

fernando$months <- month(fernando$dates_new)
fernando$months1 <- as.character(fernando$months)

head(fernando)

PlotFernandoBoxplot <- ggplot(fernando, 
                             aes(y = FERNANDO_RN, x = reorder(months1, months))) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2) +
  xlab("Months") + ylab("Wind Speed (m/s)") + 
  theme(legend.position = "none") + theme_classic()
PlotFernandoBoxplot

for (i in 1:length(fernando$FERNANDO_RN))
	{fernando$wind_speed_norm[i] = (fernando$FERNANDO_RN[i]-minimo)/(maximo-minimo)}

#PlotFernando <- ggplot(fernando) + 
#  geom_line(aes(x=dates_new, y=FERNANDO_RN)) +
#  xlab("Hours") + ylab("Wind Speed (m/s)") + 
#  theme(legend.position = "none") + theme_classic() + ylim(c(0,lim))
#PlotFernando

min(fernando$wind_speed_norm)

############Data Preparation#######################

vec_norm = as.vector(fernando$wind_speed_norm)
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
Psi_norm = CreatePsi(vec_norm)
Psi_norm
head(Psi)
tail(Psi)

cols = ncol(Psi)
rows = nrow(Psi)

d = dist(Psi_norm, method = "euclidean")

###############Silhouette, CH, DB K-Means##################

n_rows <- 12
mat = matrix(0, nrow = n_rows)
mat2 = matrix(0, nrow = n_rows)
mat3 = matrix(0, nrow = n_rows)
for (i in 2:n_rows){
   set.seed(321)
   clust = kmeans(Psi_norm, i, nstart = 25)
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
c <- pam(d, 2, diss = TRUE)
sil <- silhouette(c$cluster, d)
ch <- calinhara(Psi_norm, c$cluster)
obj <- cls.scatt.data(Psi_norm, clust$cluster, dist = "euclidean")
db <- clv.Davies.Bouldin(obj, c("average"), c("average"))

###############Silhouette, CH, DB agnes##################

n_rows <- 12
mat = matrix(0, nrow = n_rows)
mat2 = matrix(0, nrow = n_rows)
mat3 = matrix(0, nrow = n_rows)
set.seed(321)
clust = agnes(d, method = "ward")
# teste
sil <- silhouette(cutree(clust, k=2), d)
ch <- calinhara(Psi_norm, cutree(clust, k=2))
obj <- cls.scatt.data(Psi_norm, cutree(clust, k=2), dist = "euclidean")
db <- clv.Davies.Bouldin(obj, c("average"), c("average"))

for (i in 2:n_rows){
  sil = silhouette(cutree(clust, k=i), d)
  mat[i] = mean(as.matrix(sil)[,3])
  CH = calinhara(Psi_norm, cutree(clust, k=i))
  mat2[i] = as.matrix(CH)
  obj <- cls.scatt.data(Psi_norm, cutree(clust, k=i), dist = "euclidean")
  mat3[i] <- clv.Davies.Bouldin(obj, c("average"), c("average"))
}
colnames(mat) <- c("Avg_Silhouette_Value")
colnames(mat2) <- c("CH_Value")
max(mat2)
colnames(mat3) <- c("DB_Value")
min(mat3[2:12])

plot(mat)

dend <- as.dendrogram(clust, leaf_label = "none")
gg1 <- as.ggdend(dend)
ggplot(gg1, labels = FALSE, offset_labels = TRUE, theme = theme_classic()) +  
  ggtitle("Dendrogram") + ylab("Height") + xlab("Observations")





###############K-Means Best Silhouette##################

set.seed(321)
system.time(means.res <- kmeans(Psi_norm, 2))
means.res
install.packages("clv")
library(clv)
library(clusterSim)

#aaaa <- cls.scatt.data(Psi, means.res$cluster, dist = "euclidean")
#clv.Davies.Bouldin(aaaa, c("centroid", "complete", "average"), c("centroid", "complete", "average"))

#index.DB(Psi, means.res$cluster, centrotypes = "centroids")
#index.DB(Psi, means.res$cluster, d, centrotypes = "medoids")

print(means.res)
aggregate(Psi, by=list(cluster=means.res$cluster), mean)

sil = silhouette(means.res$cluster, d)
summary(sil)

plot(sil)



data(data_ratio)
cl1 <- pam(data_ratio, 4)
d<-dist(data_ratio)
print(index.DB(data_ratio, cl1$clustering,d, centrotypes="medoids"))
