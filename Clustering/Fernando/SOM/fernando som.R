rm(list=ls())

##############Packages#######################

#install.packages("dendextend")
#library(ggplot2)
library(lubridate)
#library(stats)
#library(timeDate)
#library(clue)
library(fpc)
library(clv)
#library(factoextra)
library(cluster)
library(philentropy)
library(Rcpp)
library(kohonen)

#############Opening the data###############

setwd("C:\\Users\\Rafaela\\Downloads\\2022\\TFC\\codigos eol\\Fernando")
#setwd("/home/rafaela/Documents/PUC/2022.1/TFC/codigos/Fernando")
#setwd("C:\\Users\\Rafaela\\Downloads")
getwd()

fernando <- read.csv2("hweol_FERNANDO_RN.csv")
head(fernando)

#fernando$dates_new <- as.POSIXct(fernando$dates, format = "%Y-%m-%d %H:%M")
fernando$dates_new <- as.Date(fernando$dates, format = "%Y-%m-%d %H:%M")

fernando$months <- month(fernando$dates_new)
fernando$months1 <- as.character(fernando$months)

head(fernando)
maximo = max(fernando$FERNANDO_RN)
minimo = min(fernando$FERNANDO_RN)

for (i in 1:length(fernando$FERNANDO_RN))
{fernando$wind_speed_norm[i] = (fernando$FERNANDO_RN[i]-minimo)/(maximo-minimo)}

min(fernando$FERNANDO_RN)

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

#Psi_norm_df = as.data.frame(Psi_norm)
#head(Psi_norm_df)

########## SOM ##########

#PsiTrain <- Psi[1:1000,]
#PsiTest <- Psi[1:1000,]

## and the same for rectangular maps
#n <- ceiling(rows ^ (1/2.5))
n <- 20
set.seed(321)
sommap <- supersom(Psi, grid = somgrid(n, n, "rectangular"), 
                   keep.data = TRUE, rlen = 100)

som.clust1 <- sommap$unit.classif
max(sommap$unit.classif)

plot(sommap, type="changes") # progress over iterations
plot(sommap, type="dist.neighbours", main = "SOM neighbour distances") # U-Matrix
plot(sommap, type="count", main="Node Counts") # how many samples are mapped to each node on the map
plot(sommap, type="codes")
#plot(sommap, type = "property", property = getCodes(sommap)[,4], 
#     main=colnames(getCodes(sommap))[4])

## use hierarchical clustering to cluster the codebook vectors
#som.hc <- cutree(hclust(object.distances(sommap, "codes")), 3)
#add.cluster.boundaries(sommap, som.hc)

som.events <- sommap$codes[[1]]

set.seed(321)
#som_cluster <- pam(som.events, k = 3, diss = FALSE)
som.cluster.k <- pam(som.events, k = 2, diss = FALSE)$cluster # k-means
plot(sommap, type="dist.neighbours", main = "SOM neighbour distances")
add.cluster.boundaries(sommap, som.cluster.k)

# get vector with cluster value for each original data sample
#cluster_assignment <- som_cluster[sommap$unit.classif]
# for each of analysis, add the assignment as a column in the original data:
#data$cluster <- cluster_assignment


n_rows <- 12
mat = matrix(0, nrow = n_rows)
mat2 = matrix(0, nrow = n_rows)
mat3 = matrix(0, nrow = n_rows)
for (i in 2:n_rows){
  set.seed(321)
  clust = pam(som.events, i, diss = FALSE)
  str(sil <- silhouette(clust))
  (ssi <- summary(sil))
  mat[i] = mean(as.matrix(sil)[,3])
  CH = calinhara(som.events, clust$cluster)#Psi_norm
  mat2[i] = as.matrix(CH)
  obj <- cls.scatt.data(som.events, clust$cluster, dist = "euclidean")#Psi_norm
  mat3[i] <- clv.Davies.Bouldin(obj, c("average"), c("average"))
}
colnames(mat) <- c("Avg_Silhouette_Value")
colnames(mat2) <- c("CH_Value")
max(mat2)
colnames(mat3) <- c("DB_Value")
min(mat3[2:12])

write.csv(som.cluster.k, "2k_som.csv")
write.csv(som.clust1, "clust_som.csv")
