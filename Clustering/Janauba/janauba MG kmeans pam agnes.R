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

setwd("C:\\Users\\Rafaela\\Downloads\\2022\\TFC\\codigos eol\\Janauba")
#setwd("/home/rafaela/Documents/PUC/2022.1/TFC/codigos")
getwd()

janauba <- read.csv("JANAUBA_MG_UCT.csv")
janauba <- janauba[1:359160,]
head(janauba)
tail(janauba)

#janauba$dates_new <- as.POSIXct(janauba$dates_UCT, format = "%Y/%m/%d %H:%M")
janauba$dates_new <- as.Date(janauba$dates_UCT, format = "%Y/%m/%d %H:%M")

maximo = max(janauba$SWGDN)
minimo = min(janauba$SWGDN)
#lim <- ceiling(maximo)

janauba$months <- month(janauba$dates_new)
janauba$months1 <- as.character(janauba$months)

head(janauba)

library(dplyr)
janauba_pos <- filter(janauba, SWGDN>0.1)
length(janauba_pos$SWGDN)

PlotjanaubaBoxplot <- ggplot(janauba_pos, 
                           aes(y = SWGDN, x = reorder(months1, months))) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2) +
  xlab("Months") + ylab("Global Horizontal Irradaytion (Wh/m�)") + 
  theme(legend.position = "none") + theme_classic()
PlotjanaubaBoxplot

for (i in 1:length(janauba$SWGDN))
{janauba$SWGDN_norm[i] = (janauba$SWGDN[i]-minimo)/(maximo-minimo)}

#Plotjanauba <- ggplot(janauba) + 
#  geom_line(aes(x=dates_new, y=SWGDN)) +
#  xlab("Hours") + ylab("Irradaynce Normal (kM/m2)") + 
#  theme(legend.position = "none") + theme_classic() + ylim(c(0,lim))
#Plotjanauba

min(janauba$SWGDN_norm)

############Data Preparation#######################

vec_norm = as.vector(janauba$SWGDN_norm)
vec = as.vector(janauba$SWGDN)
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

###############Silhouette, CH, DB agnes##################

n_rows <- 12
mat = matrix(0, nrow = n_rows)
mat2 = matrix(0, nrow = n_rows)
mat3 = matrix(0, nrow = n_rows)
set.seed(321)
clust = agnes(d, method = "ward")
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
print(means.res)
aggregate(Psi, by=list(cluster=means.res$cluster), mean)

sil = silhouette(means.res$cluster, d)
summary(sil)

plot(sil)
