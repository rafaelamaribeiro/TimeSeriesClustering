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
install.packages("readxl")
library(readxl)

#############Opening the data###############

#setwd("C:\\Users\\Rafaela\\Downloads\\2022\\TFC\\codigos eol\\janauba")
#setwd("/home/rafaela/Documents/PUC/2022.1/TFC/codigos")
getwd()

janauba <- read.csv2("janauba_CE_UCT.csv")
janauba <- janauba[1:359160,]
head(janauba)
tail(janauba)

#janauba$dates_new <- as.POSIXct(janauba$dates_UCT, format = "%Y/%m/%d %H:%M")
#janauba$dates_new <- as.Date(janauba$dates_UCT, format = "%Y/%m/%d %H:%M")
#janauba$dates_new2 <- as.Date(janauba$dates_UCT, format = "%Y/%m/%d")

mes <- read.csv("months.csv")
#day <- read_excel("days.xlsx")
#day <- as.vector(day)
fernando <- read.csv2("hweol_FERNANDO_RN.csv")

months <- month(mes$dates)
length(months)

janauba$dates_new <- fernando$dates

maximo = max(janauba$SWGDN)
minimo = min(janauba$SWGDN)
lim <- ceiling(maximo)

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
  xlab("Months") + ylab("Global Horizontal Irradaytion (Wh/m²)") + 
  theme(legend.position = "none") + theme_classic()
PlotjanaubaBoxplot

for (i in 1:length(janauba$SWGDN))
{janauba$wind_speed_norm[i] = (janauba$SWGDN[i]-minimo)/(maximo-minimo)}

#Plotjanauba <- ggplot(janauba) + 
#  geom_line(aes(x=dates_new, y=SWGDN)) +
#  xlab("Hours") + ylab("Irradaynce Normal (kM/m2)") + 
#  theme(legend.position = "none") + theme_classic() + ylim(c(0,lim))
#Plotjanauba

min(janauba$wind_speed_norm)

############Data Preparation#######################

vec_norm = as.vector(janauba$wind_speed_norm)
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

set.seed(321)
means.res = kmeans(Psi_norm, 3, nstart = 25)

a1 = matrix(0, means.res$size[1])
aa1 = matrix(0, means.res$size[1], 24)
ii = 1
for (j in 1:rows) {
  if (means.res$cluster[j] == 1) {
    a1[ii] = j
    aa1[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a2 = matrix(0, means.res$size[2])
aa2 = matrix(0, means.res$size[2], 24)
ii = 1
for (j in 1:rows) {
  if (means.res$cluster[j] == 2) {
    a2[ii] = j
    aa2[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a3 = matrix(0, means.res$size[3])
aa3 = matrix(0, means.res$size[3], 24)
ii = 1
for (j in 1:rows) {
  if (means.res$cluster[j] == 3) {
    a3[ii] = j
    aa3[ii,] = Psi[j,]
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

means.res$cluster[1]

for (i in 1:rows) {
  if (months[i] == 1) {
    if (means.res$cluster[i] == 1) {
      jan1 = jan1 + 1
    }
    else if (means.res$cluster[i] == 2) {
      jan2 = jan2 + 1
    }
    else if (means.res$cluster[i] == 3) {
      jan3 = jan3 + 1
    }}
  else if (months[i] == 2) {
    if (means.res$cluster[i] == 1) {
      feb1 = feb1 + 1
    }
    else if (means.res$cluster[i] == 2) {
      feb2 = feb2 + 1
    }
    else if (means.res$cluster[i] == 3) {
      feb3 = feb3 + 1
    }}
  else if (months[i] == 3) {
    if (means.res$cluster[i] == 1) {
      mar1 = mar1 + 1
    }
    else if (means.res$cluster[i] == 2) {
      mar2 = mar2 + 1
    }
    else if (means.res$cluster[i] == 3) {
      mar3 = mar3 + 1
    }}
  else if (months[i] == 4) {
    if (means.res$cluster[i] == 1) {
      apr1 = apr1 + 1
    }
    else if (means.res$cluster[i] == 2) {
      apr2 = apr2 + 1
    }
    else if (means.res$cluster[i] == 3) {
      apr3 = apr3 + 1
    }}
  else if (months[i] == 5) {
    if (means.res$cluster[i] == 1) {
      may1 = may1 + 1
    }
    else if (means.res$cluster[i] == 2) {
      may2 = may2 + 1
    }
    else if (means.res$cluster[i] == 3) {
      may3 = may3 + 1
    }}
  else if (months[i] == 6) {
    if (means.res$cluster[i] == 1) {
      jun1 = jun1 + 1
    }
    else if (means.res$cluster[i] == 2) {
      jun2 = jun2 + 1
    }
    else if (means.res$cluster[i] == 3) {
      jun3 = jun3 + 1
    }}
  else if (months[i] == 7) {
    if (means.res$cluster[i] == 1) {
      jul1 = jul1 + 1
    }
    else if (means.res$cluster[i] == 2) {
      jul2 = jul2 + 1
    }
    else if (means.res$cluster[i] == 3) {
      jul3 = jul3 + 1
    }}
  else if (months[i] == 8) {
    if (means.res$cluster[i] == 1) {
      aug1 = aug1 + 1
    }
    else if (means.res$cluster[i] == 2) {
      aug2 = aug2 + 1
    }
    else if (means.res$cluster[i] == 3) {
      aug3 = aug3 + 1
    }}
  else if (months[i] == 9) {
    if (means.res$cluster[i] == 1) {
      sep1 = sep1 + 1
    }
    else if (means.res$cluster[i] == 2) {
      sep2 = sep2 + 1
    }
    else if (means.res$cluster[i] == 3) {
      sep3 = sep3 + 1
    }}
  else if (months[i] == 10) {
    if (means.res$cluster[i] == 1) {
      oct1 = oct1 + 1
    }
    else if (means.res$cluster[i] == 2) {
      oct2 = oct2 + 1
    }
    else if (means.res$cluster[i] == 3) {
      oct3 = oct3 + 1
    }}
  else if (months[i] == 11) {
    if (means.res$cluster[i] == 1) {
      nov1 = nov1 + 1
    }
    else if (means.res$cluster[i] == 2) {
      nov2 = nov2 + 1
    }
    else if (means.res$cluster[i] == 3) {
      nov3 = nov3 + 1
    }}
  else if (months[i] == 12) {
    if (means.res$cluster[i] == 1) {
      dec1 = dec1 + 1
    }
    else if (means.res$cluster[i] == 2) {
      dec2 = dec2 + 1
    }
    else if (means.res$cluster[i] == 3) {
      dec3 = dec3 + 1
    }}
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
write.csv(df,"porMes_janauba_k_means_3K.csv", row.names = FALSE)
write.csv(aa1,"profile1_janauba_k_means_3K.csv", row.names = FALSE)
write.csv(aa2,"profile2_janauba_k_means_3K.csv", row.names = FALSE)
write.csv(aa3,"profile3_janauba_k_means_3K.csv", row.names = FALSE)
write.csv(means.res$cluster, "cluster_janauba_k_means_3K.csv", row.names = FALSE)
