rm(list=ls())

##############Packages#######################

library(lubridate)
library(fpc)
library(clv)
library(cluster)
library(philentropy)
library(Rcpp)
library(kohonen)

#############Opening the data###############

setwd("C:\\Users\\Rafaela\\Downloads\\2022\\TFC\\codigos eol\\Cidreira")
#setwd("/home/rafaela/Documents/PUC/2022.1/TFC/codigos/Cidreira")
#setwd("C:\\Users\\Rafaela\\Downloads")
getwd()

cidreira <- read.csv2("hweol_CIDREIRA_RS.csv")
head(cidreira)

#cidreira$dates_new <- as.POSIXct(cidreira$dates, format = "%Y-%m-%d %H:%M")
cidreira$dates_new <- as.Date(cidreira$dates, format = "%Y-%m-%d %H:%M")

cidreira$months <- month(cidreira$dates_new)
cidreira$months1 <- as.character(cidreira$months)

head(cidreira)
maximo = max(cidreira$CIDREIRA_RS)
minimo = min(cidreira$CIDREIRA_RS)

for (i in 1:length(cidreira$CIDREIRA_RS))
{cidreira$wind_speed_norm[i] = (cidreira$CIDREIRA_RS[i]-minimo)/(maximo-minimo)}

min(cidreira$wind_speed_norm)

############Data Preparation#######################

vec_norm = as.vector(cidreira$wind_speed_norm)
vec = as.vector(cidreira$CIDREIRA_RS)
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
n <- 4
set.seed(321)
sommap <- supersom(Psi_norm, grid = somgrid(n, n, "rectangular"), 
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
som.cluster.k <- pam(som.events, k = 15, diss = FALSE)$cluster # pam
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
