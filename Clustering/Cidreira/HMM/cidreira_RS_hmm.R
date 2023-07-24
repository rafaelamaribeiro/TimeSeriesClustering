rm(list=ls())

##############Packages#######################

#install.packages("dendextend")
library(ggplot2)
library(lubridate)
library(stats)
library(timeDate)
library(clue)
library(fpc)
library(clv)
library(factoextra)
library(LaplacesDemon)
library("depmixS4")
library(cluster)
library(philentropy)
library(Rcpp)

#############Opening the data###############

#setwd("C:\\Users\\Rafaela\\Downloads\\2022\\TFC\\codigos eol\\Cidreira")
#setwd("/home/rafaela/Documents/PUC/2022.1/TFC/codigos/Cidreira")
#setwd("C:\\Users\\Rafaela\\Downloads")
getwd()

cidreira <- read.csv2("hweol_CIDREIRA_RS.csv")
head(cidreira)

#cidreira$dates_new <- as.POSIXct(cidreira$dates, format = "%Y-%m-%d %H:%M")
cidreira$dates_new <- as.Date(cidreira$dates, format = "%Y-%m-%d %H:%M")

maximo = max(cidreira$CIDREIRA_RS)
minimo = min(cidreira$CIDREIRA_RS)
#lim <- ceiling(maximo)

cidreira$months <- month(cidreira$dates_new)
cidreira$months1 <- as.character(cidreira$months)

head(cidreira)

#PlotcidreiraBoxplot <- ggplot(cidreira, 
#                              aes(y = cidreira_RS, x = reorder(months1, months))) + 
#  geom_boxplot(outlier.colour="black", outlier.shape=16,
#               outlier.size=2) +
#  xlab("Months") + ylab("Wind Speed (m/s)") + 
#  theme(legend.position = "none") + theme_classic()
#PlotcidreiraBoxplot

for (i in 1:length(cidreira$CIDREIRA_RS))
{cidreira$wind_speed_norm[i] = (cidreira$CIDREIRA_RS[i]-minimo)/(maximo-minimo)}

#Plotcidreira <- ggplot(cidreira) + 
#  geom_line(aes(x=dates_new, y=cidreira_RS)) +
#  xlab("Hours") + ylab("Wind Speed (m/s)") + 
#  theme(legend.position = "none") + theme_classic() + ylim(c(0,lim))
#Plotcidreira

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

Psi_norm_df = as.data.frame(Psi_norm)
head(Psi_norm_df)

###############hidden states##################
# change nstates from 2 to 5

#response = list(V1~1, V2~1, V3~1, V4~1, V5~1,
#V6~1, V7~1, V8~1, V9~1, V10~1,
#V11~1, V12~1, V13~1, V14~1, V15~1,
#V16~1, V17~1, V18~1, V19~1, V20~1,
#V21~1, V22~1, V23~1, V24~1)
#family = list(gaussian(), gaussian(), gaussian(), gaussian(),
#              gaussian(), gaussian(), gaussian(), gaussian(),
#              gaussian(), gaussian(), gaussian(), gaussian(),
#              gaussian(), gaussian(), gaussian(), gaussian(),
#              gaussian(), gaussian(), gaussian(), gaussian(),
#              gaussian(), gaussian(), gaussian(), gaussian())

calc_mod_fit <- function(r, s) {
  d <- t(Psi_norm_df[r,])
  colnames(d) <- "day"
  d <- as.data.frame(d)
  mod <- depmix(response = day~1,
                data = d, nstates = s, 
                family = gaussian())

  #tryCatch(fit(mod), finally = print("aaaaa"))
  log_100 <- vector()
  bic_100 <- vector()
  aic_100 <- vector()
  iter_100 <- vector()
  funcionaram <- vector()
  set.seed(123)
  iter_100 <- replicate(25, try(fit(mod, verbose = F, 
                                     emcontrol = em.control(tol = 1e-3))), 
                        simplify = "array")
  # https://stackoverflow.com/questions/25363871/na-nan-inf-error-when-fitting-hmm-using-depmixs4-in-r
  iter_100_array <- as.array(iter_100)
  funcionaram <- iter_100_array[grepl("log", iter_100_array)]
  for (i in 1:length(funcionaram)) {
    bic_100[i] <- BIC(funcionaram[[i]]); 
    log_100[i] <- logLik(funcionaram[[i]]);
    aic_100[i] <- AIC(funcionaram[[i]])
  }
  mean(log_100)
  mean(bic_100)
  mean(aic_100)
  #fm <- fit(mod)
  #summary(fm)
  #BIC(fm)
  #logLik(fm)[1]
  return(list(mean(log_100), mean(bic_100), mean(aic_100)))
}

#test <- calc_mod_fit(65,3)

vector_log_bic <- vector()
vector_log_aic <- vector()
vector_bic <- vector()
vector_aic <- vector()
start <- Sys.time()
#vec_log_2st <- vector()
for (r in 11932:rows) {
  #print(r)
  max_bic <- Inf
  max_aic <- Inf
  #r<-65
  #r<-484
  #r<- 1665
  #r<-7004
  #r<-11931
  #s<-2
  for (s in 2:3) {
    log_bic_aic <- calc_mod_fit(r,s)
    if (log_bic_aic[[2]] < max_bic) {
      vector_log_bic[r] <- log_bic_aic[[1]];
      #vec_log_2st[r] <- log_bic_aic[[1]]
      vector_bic[r] <- s;
      max_bic <- log_bic_aic[[2]]  
    }
    if (log_bic_aic[[3]] < max_aic) {
      vector_log_aic[r] <- log_bic_aic[[1]];
      vector_aic[r] <- s;
      max_aic <- log_bic_aic[[3]]
    }
  }
}
end <- Sys.time()
#end-start
#tail(vector_aic)

#length(vec_log_2st)

write.csv(vector_log_aic, "log_aic_cidreira.csv")
write.csv(vector_log_bic, "log_bic_cidreira.csv")

vector_aic == vector_bic

matriz <- cbind(vector_log_bic, vector_log_aic)
d <- KL(matriz, test.na = TRUE, unit = "log2", est.prob = "empirical", epsilon = 1e-05)
isSymmetric(d)

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
c <- pam(as.dist.default(KLDfinal$mean.KLD), 2, diss = TRUE)
sil <- silhouette(c$cluster, d)
ch <- calinhara(Psi_norm, c$cluster)
obj <- cls.scatt.data(Psi_norm, clust$cluster, dist = "euclidean")
db <- clv.Davies.Bouldin(obj, c("average"), c("average"))

# Os centros

set.seed(321)
clust = pam(d, 2, diss = TRUE)

summary(clust)

dates <- unique(cidreira$dates_new)
tail(dates)
meses <- month(dates)
length(meses)

a1 = matrix(0, clust$clusinfo[1,1])
aa1 = matrix(0, clust$clusinfo[1,1], 24)
ii = 1
for (j in 1:rows) {
  if (clust$clustering[j] == 1) {
    a1[ii] = j
    aa1[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a2 = matrix(0, clust$clusinfo[2,1])
aa2 = matrix(0, clust$clusinfo[2,1], 24)
ii = 1
for (j in 1:rows) {
  if (clust$clustering[j] == 2) {
    a2[ii] = j
    aa2[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a3 = matrix(0, clust$clusinfo[3,1])
aa3 = matrix(0, clust$clusinfo[3,1], 24)
ii = 1
for (j in 1:rows) {
  if (clust$clustering[j] == 3) {
    a3[ii] = j
    aa3[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a4 = matrix(0, clust$clusinfo[4,1])
aa4 = matrix(0, clust$clusinfo[4,1], 24)
ii = 1
for (j in 1:rows) {
  if (clust$clustering[j] == 4) {
    a4[ii] = j
    aa4[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a5 = matrix(0, clust$clusinfo[5,1])
aa5 = matrix(0, clust$clusinfo[5,1], 24)
ii = 1
for (j in 1:rows) {
  if (clust$clustering[j] == 5) {
    a5[ii] = j
    aa5[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a6 = matrix(0, clust$clusinfo[6,1])
aa6 = matrix(0, clust$clusinfo[6,1], 24)
ii = 1
for (j in 1:rows) {
  if (clust$clustering[j] == 6) {
    a6[ii] = j
    aa6[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a7 = matrix(0, clust$clusinfo[7,1])
aa7 = matrix(0, clust$clusinfo[7,1], 24)
ii = 1
for (j in 1:rows) {
  if (clust$clustering[j] == 7) {
    a7[ii] = j
    aa7[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a8 = matrix(0, clust$clusinfo[8,1])
aa8 = matrix(0, clust$clusinfo[8,1], 24)
ii = 1
for (j in 1:rows) {
  if (clust$clustering[j] == 8) {
    a8[ii] = j
    aa8[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a9 = matrix(0, clust$clusinfo[9,1])
aa9 = matrix(0, clust$clusinfo[9,1], 24)
ii = 1
for (j in 1:rows) {
  if (clust$clustering[j] == 9) {
    a9[ii] = j
    aa9[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a10 = matrix(0, clust$clusinfo[10,1])
aa10 = matrix(0, clust$clusinfo[10,1], 24)
ii = 1
for (j in 1:rows) {
  if (clust$clustering[j] == 10) {
    a10[ii] = j
    aa10[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a11 = matrix(0, clust$clusinfo[11,1])
aa11 = matrix(0, clust$clusinfo[11,1], 24)
ii = 1
for (j in 1:rows) {
  if (clust$clustering[j] == 11) {
    a11[ii] = j
    aa11[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a12 = matrix(0, clust$clusinfo[12,1])
aa12 = matrix(0, clust$clusinfo[12,1], 24)
ii = 1
for (j in 1:rows) {
  if (clust$clustering[j] == 12) {
    a12[ii] = j
    aa12[ii,] = Psi[j,]
    ii = ii + 1
  }
}

# january
jan1 = 0
jan2 = 0
jan3 = 0
jan4 = 0
jan5 = 0
jan6 = 0
jan7 = 0
jan8 = 0
jan9 = 0
jan10 = 0
jan11 = 0
jan12 = 0

# february
feb1 = 0
feb2 = 0
feb3 = 0
feb4 = 0
feb5 = 0
feb6 = 0
feb7 = 0
feb8 = 0
feb9 = 0
feb10 = 0
feb11 = 0
feb12 = 0

# march
mar1 = 0
mar2 = 0
mar3 = 0
mar4 = 0
mar5 = 0
mar6 = 0
mar7 = 0
mar8 = 0
mar9 = 0
mar10 = 0
mar11 = 0
mar12 = 0

# april
apr1 = 0
apr2 = 0
apr3 = 0
apr4 = 0
apr5 = 0
apr6 = 0
apr7 = 0
apr8 = 0
apr9 = 0
apr10 = 0
apr11 = 0
apr12 = 0

# may
may1 = 0
may2 = 0
may3 = 0
may4 = 0
may5 = 0
may6 = 0
may7 = 0
may8 = 0
may9 = 0
may10 = 0
may11 = 0
may12 = 0

# june
jun1 = 0
jun2 = 0
jun3 = 0
jun4 = 0
jun5 = 0
jun6 = 0
jun7 = 0
jun8 = 0
jun9 = 0
jun10 = 0
jun11 = 0
jun12 = 0

# july
jul1 = 0
jul2 = 0
jul3 = 0
jul4 = 0
jul5 = 0
jul6 = 0
jul7 = 0
jul8 = 0
jul9 = 0
jul10 = 0
jul11 = 0
jul12 = 0

# august
aug1 = 0
aug2 = 0
aug3 = 0
aug4 = 0
aug5 = 0
aug6 = 0
aug7 = 0
aug8 = 0
aug9 = 0
aug10 = 0
aug11 = 0
aug12 = 0

# september
sep1 = 0
sep2 = 0
sep3 = 0
sep4 = 0
sep5 = 0
sep6 = 0
sep7 = 0
sep8 = 0
sep9 = 0
sep10 = 0
sep11 = 0
sep12 = 0

# october
oct1 = 0
oct2 = 0
oct3 = 0
oct4 = 0
oct5 = 0
oct6 = 0
oct7 = 0
oct8 = 0
oct9 = 0
oct10 = 0
oct11 = 0
oct12 = 0

# november
nov1 = 0
nov2 = 0
nov3 = 0
nov4 = 0
nov5 = 0
nov6 = 0
nov7 = 0
nov8 = 0
nov9 = 0
nov10 = 0
nov11 = 0
nov12 = 0

# december
dec1 = 0
dec2 = 0
dec3 = 0
dec4 = 0
dec5 = 0
dec6 = 0
dec7 = 0
dec8 = 0
dec9 = 0
dec10 = 0
dec11 = 0
dec12 = 0

for (i in 1:rows) {
  if (meses[i] == 1) {
    if (i %in% a1) {
      jan1 = jan1 + 1}
    else if (i %in% a2) {
      jan2 = jan2 + 1}
    else if (i %in% a3) {
      jan3 = jan3 + 1}
    else if (i %in% a4) {
      jan4 = jan4 + 1}
    else if (i %in% a5) {
      jan5 = jan5 + 1}
    else if (i %in% a6) {
      jan6 = jan6 + 1}
    else if (i %in% a7) {
      jan7 = jan7 + 1}
    else if (i %in% a8) {
      jan8 = jan8 + 1}
    else if (i %in% a9) {
      jan9 = jan9 + 1}
    else if (i %in% a10) {
      jan10 = jan10 + 1}
    else if (i %in% a11) {
      jan11 = jan11 + 1}
    else if (i %in% a12) {
      jan12 = jan12 + 1}}
  else if (meses[i] == 2) {
    if (i %in% a1) {
      feb1 = feb1 + 1}
    else if (i %in% a2) {
      feb2 = feb2 + 1}
    else if (i %in% a3) {
      feb3 = feb3 + 1}
    else if (i %in% a4) {
      feb4 = feb4 + 1}
    else if (i %in% a5) {
      feb5 = feb5 + 1}
    else if (i %in% a6) {
      feb6 = feb6 + 1}
    else if (i %in% a7) {
      feb7 = feb7 + 1}
    else if (i %in% a8) {
      feb8 = feb8 + 1}
    else if (i %in% a9) {
      feb9 = feb9 + 1}
    else if (i %in% a10) {
      feb10 = feb10 + 1}
    else if (i %in% a11) {
      feb11 = feb11 + 1}
    else if (i %in% a12) {
      feb12 = feb12 + 1}}
  else if (meses[i] == 3) {
    if (i %in% a1) {
      mar1 = mar1 + 1}
    else if (i %in% a2) {
      mar2 = mar2 + 1}
    else if (i %in% a3) {
      mar3 = mar3 + 1}
    else if (i %in% a4) {
      mar4 = mar4 + 1}
    else if (i %in% a5) {
      mar5 = mar5 + 1}
    else if (i %in% a6) {
      mar6 = mar6 + 1}
    else if (i %in% a7) {
      mar7 = mar7 + 1}
    else if (i %in% a8) {
      mar8 = mar8 + 1}
    else if (i %in% a9) {
      mar9 = mar9 + 1}
    else if (i %in% a10) {
      mar10 = mar10 + 1}
    else if (i %in% a11) {
      mar11 = mar11 + 1}
    else if (i %in% a12) {
      mar12 = mar12 + 1}}
  else if (meses[i] == 4) {
    if (i %in% a1) {
      apr1 = apr1 + 1}
    else if (i %in% a2) {
      apr2 = apr2 + 1}
    else if (i %in% a3) {
      apr3 = apr3 + 1}
    else if (i %in% a4) {
      apr4 = apr4 + 1}
    else if (i %in% a5) {
      apr5 = apr5 + 1}
    else if (i %in% a6) {
      apr6 = apr6 + 1}
    else if (i %in% a7) {
      apr7 = apr7 + 1}
    else if (i %in% a8) {
      apr8 = apr8 + 1}
    else if (i %in% a9) {
      apr9 = apr9 + 1}
    else if (i %in% a10) {
      apr10 = apr10 + 1}
    else if (i %in% a11) {
      apr11 = apr11 + 1}
    else if (i %in% a12) {
      apr12 = apr12 + 1}}
  else if (meses[i] == 5) {
    if (i %in% a1) {
      may1 = may1 + 1}
    else if (i %in% a2) {
      may2 = may2 + 1}
    else if (i %in% a3) {
      may3 = may3 + 1}
    else if (i %in% a4) {
      may4 = may4 + 1}
    else if (i %in% a5) {
      may5 = may5 + 1}
    else if (i %in% a6) {
      may6 = may6 + 1}
    else if (i %in% a7) {
      may7 = may7 + 1}
    else if (i %in% a8) {
      may8 = may8 + 1}
    else if (i %in% a9) {
      may9 = may9 + 1}
    else if (i %in% a10) {
      may10 = may10 + 1}
    else if (i %in% a11) {
      may11 = may11 + 1}
    else if (i %in% a12) {
      may12 = may12 + 1}}
  else if (meses[i] == 6) {
    if (i %in% a1) {
      jun1 = jun1 + 1}
    else if (i %in% a2) {
      jun2 = jun2 + 1}
    else if (i %in% a3) {
      jun3 = jun3 + 1}
    else if (i %in% a4) {
      jun4 = jun4 + 1}
    else if (i %in% a5) {
      jun5 = jun5 + 1}
    else if (i %in% a6) {
      jun6 = jun6 + 1}
    else if (i %in% a7) {
      jun7 = jun7 + 1}
    else if (i %in% a8) {
      jun8 = jun8 + 1}
    else if (i %in% a9) {
      jun9 = jun9 + 1}
    else if (i %in% a10) {
      jun10 = jun10 + 1}
    else if (i %in% a11) {
      jun11 = jun11 + 1}
    else if (i %in% a12) {
      jun12 = jun12 + 1}}
  else if (meses[i] == 7) {
    if (i %in% a1) {
      jul1 = jul1 + 1}
    else if (i %in% a2) {
      jul2 = jul2 + 1}
    else if (i %in% a3) {
      jul3 = jul3 + 1}
    else if (i %in% a4) {
      jul4 = jul4 + 1}
    else if (i %in% a5) {
      jul5 = jul5 + 1}
    else if (i %in% a6) {
      jul6 = jul6 + 1}
    else if (i %in% a7) {
      jul7 = jul7 + 1}
    else if (i %in% a8) {
      jul8 = jul8 + 1}
    else if (i %in% a9) {
      jul9 = jul9 + 1}
    else if (i %in% a10) {
      jul10 = jul10 + 1}
    else if (i %in% a11) {
      jul11 = jul11 + 1}
    else if (i %in% a12) {
      jul12 = jul12 + 1}}
  else if (meses[i] == 8) {
    if (i %in% a1) {
      aug1 = aug1 + 1}
    else if (i %in% a2) {
      aug2 = aug2 + 1}
    else if (i %in% a3) {
      aug3 = aug3 + 1}
    else if (i %in% a4) {
      aug4 = aug4 + 1}
    else if (i %in% a5) {
      aug5 = aug5 + 1}
    else if (i %in% a6) {
      aug6 = aug6 + 1}
    else if (i %in% a7) {
      aug7 = aug7 + 1}
    else if (i %in% a8) {
      aug8 = aug8 + 1}
    else if (i %in% a9) {
      aug9 = aug9 + 1}
    else if (i %in% a10) {
      aug10 = aug10 + 1}
    else if (i %in% a11) {
      aug11 = aug11 + 1}
    else if (i %in% a12) {
      aug12 = aug12 + 1}}
  else if (meses[i] == 9) {
    if (i %in% a1) {
      sep1 = sep1 + 1}
    else if (i %in% a2) {
      sep2 = sep2 + 1}
    else if (i %in% a3) {
      sep3 = sep3 + 1}
    else if (i %in% a4) {
      sep4 = sep4 + 1}
    else if (i %in% a5) {
      sep5 = sep5 + 1}
    else if (i %in% a6) {
      sep6 = sep6 + 1}
    else if (i %in% a7) {
      sep7 = sep7 + 1}
    else if (i %in% a8) {
      sep8 = sep8 + 1}
    else if (i %in% a9) {
      sep9 = sep9 + 1}
    else if (i %in% a10) {
      sep10 = sep10 + 1}
    else if (i %in% a11) {
      sep11 = sep11 + 1}
    else if (i %in% a12) {
      sep12 = sep12 + 1}}
  else if (meses[i] == 10) {
    if (i %in% a1) {
      oct1 = oct1 + 1}
    else if (i %in% a2) {
      oct2 = oct2 + 1}
    else if (i %in% a3) {
      oct3 = oct3 + 1}
    else if (i %in% a4) {
      oct4 = oct4 + 1}
    else if (i %in% a5) {
      oct5 = oct5 + 1}
    else if (i %in% a6) {
      oct6 = oct6 + 1}
    else if (i %in% a7) {
      oct7 = oct7 + 1}
    else if (i %in% a8) {
      oct8 = oct8 + 1}
    else if (i %in% a9) {
      oct9 = oct9 + 1}
    else if (i %in% a10) {
      oct10 = oct10 + 1}
    else if (i %in% a11) {
      oct11 = oct11 + 1}
    else if (i %in% a12) {
      oct12 = oct12 + 1}}
  else if (meses[i] == 11) {
    if (i %in% a1) {
      nov1 = nov1 + 1}
    else if (i %in% a2) {
      nov2 = nov2 + 1}
    else if (i %in% a3) {
      nov3 = nov3 + 1}
    else if (i %in% a4) {
      nov4 = nov4 + 1}
    else if (i %in% a5) {
      nov5 = nov5 + 1}
    else if (i %in% a6) {
      nov6 = nov6 + 1}
    else if (i %in% a7) {
      nov7 = nov7 + 1}
    else if (i %in% a8) {
      nov8 = nov8 + 1}
    else if (i %in% a9) {
      nov9 = nov9 + 1}
    else if (i %in% a10) {
      nov10 = nov10 + 1}
    else if (i %in% a11) {
      nov11 = nov11 + 1}
    else if (i %in% a12) {
      nov12 = nov12 + 1}}
  else if (meses[i] == 12) {
    if (i %in% a1) {
      dec1 = dec1 + 1}
    else if (i %in% a2) {
      dec2 = dec2 + 1}
    else if (i %in% a3) {
      dec3 = dec3 + 1}
    else if (i %in% a4) {
      dec4 = dec4 + 1}
    else if (i %in% a5) {
      dec5 = dec5 + 1}
    else if (i %in% a6) {
      dec6 = dec6 + 1}
    else if (i %in% a7) {
      dec7 = dec7 + 1}
    else if (i %in% a8) {
      dec8 = dec8 + 1}
    else if (i %in% a9) {
      dec9 = dec9 + 1}
    else if (i %in% a10) {
      dec10 = dec10 + 1}
    else if (i %in% a11) {
      dec11 = dec11 + 1}
    else if (i %in% a12) {
      dec12 = dec12 + 1}}
}

January = c(jan1, jan2, jan3, jan4, jan5, jan6, jan7, jan8, jan9, jan10, jan11, jan12)
February = c(feb1, feb2, feb3, feb4, feb5, feb6, feb7, feb8, feb9, feb10, feb11, feb12)
March = c(mar1, mar2, mar3, mar4, mar5, mar6, mar7, mar8, mar9, mar10, mar11, mar12)
April = c(apr1, apr2, apr3, apr4, apr5, apr6, apr7, apr8, apr9, apr10, apr11, apr12)
May = c(may1, may2, may3, may4, may5, may6, may7, may8, may9, may10, may11, may12)
June = c(jun1, jun2, jun3, jun4, jun5, jun6, jun7, jun8, jun9, jun10, jun11, jun12)
July = c(jul1, jul2, jul3, jul4, jul5, jul6, jul7, jul8, jul9, jul10, jul11, jul12)
August = c(aug1, aug2, aug3, aug4, aug5, aug6, aug7, aug8, aug9, aug10, aug11, aug12)
September = c(sep1, sep2, sep3, sep4, sep5, sep6, sep7, sep8, sep9, sep10, sep11, sep12)
October = c(oct1, oct2, oct3, oct4, oct5, oct6, oct7, oct8, oct9, oct10, oct11, oct12)
November = c(nov1, nov2, nov3, nov4, nov5, nov6, nov7, nov8, nov9, nov10, nov11, nov12)
December = c(dec1, dec2, dec3, dec4, dec5, dec6, dec7, dec8, dec9, dec10, dec11, dec12)

df = data.frame(January, February, March, April, May, June, July, 
                August, September, October, November, December)

sum(df)

df
write.csv(df,"byMonth_cidreira_hmm_2K.csv", row.names = FALSE)
write.csv(aa1,"profile1_cidreira_hmm_2K.csv", row.names = FALSE)
write.csv(aa2,"profile2_cidreira_hmm_2K.csv", row.names = FALSE)
write.csv(aa3,"profile3_cidreira_hmm_12K.csv", row.names = FALSE)
write.csv(aa4,"profile4_cidreira_hmm_12K.csv", row.names = FALSE)
write.csv(aa5,"profile5_cidreira_hmm_12K.csv", row.names = FALSE)
write.csv(aa6,"profile6_cidreira_hmm_12K.csv", row.names = FALSE)
write.csv(aa7,"profile7_cidreira_hmm_12K.csv", row.names = FALSE)
write.csv(aa8,"profile8_cidreira_hmm_12K.csv", row.names = FALSE)
write.csv(aa9,"profile9_cidreira_hmm_12K.csv", row.names = FALSE)
write.csv(aa10,"profile10_cidreira_hmm_12K.csv", row.names = FALSE)
write.csv(aa11,"profile11_cidreira_hmm_12K.csv", row.names = FALSE)
write.csv(aa12,"profile12_cidreira_hmm_12K.csv", row.names = FALSE)
write.csv(clust$clustering, "cluster_cidreira_hmm_2.csv", row.names = FALSE)

