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

floresta <- read.csv2("FLORESTA_RN_UCT.csv")
floresta <- floresta[1:359160,]
head(floresta)
tail(floresta)

#floresta$dates_new <- as.POSIXct(floresta$dates_UCT, format = "%Y/%m/%d %H:%M")
#floresta$dates_new <- as.Date(floresta$dates_UCT, format = "%Y/%m/%d %H:%M")

mes <- read.csv("months.csv")
#day <- read_excel("days.xlsx")
#day <- as.vector(day)
fernando <- read.csv2("hweol_FERNANDO_RN.csv")

floresta$dates_new <- fernando$dates

#colors <- c("Global Horizontal Irradaytion" = "RoyalBlue", 
#            "Direct Normal Irradaytion" = "DeepPink")

#ggplot(floresta[1:24,], aes(x = dates_new)) +
#  geom_line(aes(y = SWGDN, color = "Global Horizontal Irradaytion"), size = 1) +
#  geom_line(aes(y = SWTDN, color = "Direct Normal Irradaytion"), size = 1) +
#  labs(x = "Hours",
#       y = "Irradaytions (Wh/m?)",
#       color = "Legend") +
#  scale_color_manual(values = colors) + theme_classic()

maximo = max(floresta$SWGDN)
minimo = min(floresta$SWGDN)
lim <- ceiling(maximo)
sort(floresta$SWGDN, decreasing = TRUE)

floresta$months <- month(floresta$dates_new)
floresta$months1 <- as.character(floresta$months)

head(floresta)

library(dplyr)
floresta_pos <- filter(floresta, SWGDN>0.1)
length(floresta_pos$SWGDN)

PlotapodiBoxplot <- ggplot(floresta_pos, 
                           aes(y = SWGDN, x = reorder(months1, months))) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2) +
  xlab("Months") + ylab("Global Horizontal Irradaytion (Wh/m²)") + 
  theme(legend.position = "none") + theme_classic()
PlotapodiBoxplot

for (i in 1:length(floresta$SWGDN))
  {floresta$wind_speed_norm[i] = (floresta$SWGDN[i]-minimo)/(maximo-minimo)}

#Plotfloresta <- ggplot(floresta) + 
#  geom_line(aes(x=dates_new, y=SWGDN)) +
#  xlab("Hours") + ylab("Irradaynce Normal (kM/m2)") + 
#  theme(legend.position = "none") + theme_classic() + ylim(c(0,lim))
#Plotfloresta

min(floresta$wind_speed_norm)

############Data Preparation#######################

vec_norm = as.vector(floresta$wind_speed_norm)
vec = as.vector(floresta$SWGDN)
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

###############10 clusters###############

set.seed(321)
pam.res = pam(d, 10, diss = TRUE)

#dates <- unique(apodi$dates_new2)
#tail(dates)
#mes <- as.POSIXct(mes$dates, format = "%Y/%m/%d")
months <- month(mes$dates)
length(months)

a1 = matrix(0, pam.res$clusinfo[1,1])
aa1 = matrix(0, pam.res$clusinfo[1,1], 24)
ii = 1
for (j in 1:rows) {
  if (pam.res$clustering[j] == 1) {
    a1[ii] = j
    aa1[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a2 = matrix(0, pam.res$clusinfo[2,1])
aa2 = matrix(0, pam.res$clusinfo[2,1], 24)
ii = 1
for (j in 1:rows) {
  if (pam.res$clustering[j] == 2) {
    a2[ii] = j
    aa2[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a3 = matrix(0, pam.res$clusinfo[3,1])
aa3 = matrix(0, pam.res$clusinfo[3,1], 24)
ii = 1
for (j in 1:rows) {
  if (pam.res$clustering[j] == 3) {
    a3[ii] = j
    aa3[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a4 = matrix(0, pam.res$clusinfo[4,1])
aa4 = matrix(0, pam.res$clusinfo[4,1], 24)
ii = 1
for (j in 1:rows) {
  if (pam.res$clustering[j] == 4) {
    a4[ii] = j
    aa4[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a5 = matrix(0, pam.res$clusinfo[5,1])
aa5 = matrix(0, pam.res$clusinfo[5,1], 24)
ii = 1
for (j in 1:rows) {
  if (pam.res$clustering[j] == 5) {
    a5[ii] = j
    aa5[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a6 = matrix(0, pam.res$clusinfo[6,1])
aa6 = matrix(0, pam.res$clusinfo[6,1], 24)
ii = 1
for (j in 1:rows) {
  if (pam.res$clustering[j] == 6) {
    a6[ii] = j
    aa6[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a7 = matrix(0, pam.res$clusinfo[7,1])
aa7 = matrix(0, pam.res$clusinfo[7,1], 24)
ii = 1
for (j in 1:rows) {
  if (pam.res$clustering[j] == 7) {
    a7[ii] = j
    aa7[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a8 = matrix(0, pam.res$clusinfo[8,1])
aa8 = matrix(0, pam.res$clusinfo[8,1], 24)
ii = 1
for (j in 1:rows) {
  if (pam.res$clustering[j] == 8) {
    a8[ii] = j
    aa8[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a9 = matrix(0, pam.res$clusinfo[9,1])
aa9 = matrix(0, pam.res$clusinfo[9,1], 24)
ii = 1
for (j in 1:rows) {
  if (pam.res$clustering[j] == 9) {
    a9[ii] = j
    aa9[ii,] = Psi[j,]
    ii = ii + 1
  }
}

a10 = matrix(0, pam.res$clusinfo[10,1])
aa10 = matrix(0, pam.res$clusinfo[10,1], 24)
ii = 1
for (j in 1:rows) {
  if (pam.res$clustering[j] == 10) {
    a10[ii] = j
    aa10[ii,] = Psi[j,]
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

for (i in 1:rows) {
    if (months[i] == 1) {
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
      jan10 = jan10 + 1}}
 else if (months[i] == 2) {
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
      feb10 = feb10 + 1}}
 else if (months[i] == 3) {
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
     mar10 = mar10 + 1}}
 else if (months[i] == 4) {
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
     apr10 = apr10 + 1}}
 else if (months[i] == 5) {
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
     may10 = may10 + 1}}
 else if (months[i] == 6) {
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
     jun10 = jun10 + 1}}
 else if (months[i] == 7) {
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
     jul10 = jul10 + 1}}
 else if (months[i] == 8) {
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
     aug10 = aug10 + 1}}
 else if (months[i] == 9) {
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
     sep10 = sep10 + 1}}
 else if (months[i] == 10) {
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
     oct10 = oct10 + 1}}
 else if (months[i] == 11) {
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
     nov10 = nov10 + 1}}
 else if (months[i] == 12) {
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
     dec10 = dec10 + 1}}
}

January = c(jan1, jan2, jan3, jan4, jan5, jan6, jan7, jan8, jan9, jan10)
February = c(feb1, feb2, feb3, feb4, feb5, feb6, feb7, feb8, feb9, feb10)
March = c(mar1, mar2, mar3, mar4, mar5, mar6, mar7, mar8, mar9, mar10)
April = c(apr1, apr2, apr3, apr4, apr5, apr6, apr7, apr8, apr9, apr10)
May = c(may1, may2, may3, may4, may5, may6, may7, may8, may9, may10)
June = c(jun1, jun2, jun3, jun4, jun5, jun6, jun7, jun8, jun9, jun10)
July = c(jul1, jul2, jul3, jul4, jul5, jul6, jul7, jul8, jul9, jul10)
August = c(aug1, aug2, aug3, aug4, aug5, aug6, aug7, aug8, aug9, aug10)
September = c(sep1, sep2, sep3, sep4, sep5, sep6, sep7, sep8, sep9, sep10)
October = c(oct1, oct2, oct3, oct4, oct5, oct6, oct7, oct8, oct9, oct10)
November = c(nov1, nov2, nov3, nov4, nov5, nov6, nov7, nov8, nov9, nov10)
December = c(dec1, dec2, dec3, dec4, dec5, dec6, dec7, dec8, dec9, dec10)

df = data.frame(January, February, March, April, May, June, July, 
                August, September, October, November, December)

sum(df)

df
write.csv(df,"porMes_floresta.csv", row.names = FALSE)
write.csv(aa1,"profile1_floresta.csv", row.names = FALSE)
write.csv(aa2,"profile2_floresta.csv", row.names = FALSE)
write.csv(aa3,"profile3_floresta.csv", row.names = FALSE)
write.csv(aa4,"profile4_floresta.csv", row.names = FALSE)
write.csv(aa5,"profile5_floresta.csv", row.names = FALSE)
write.csv(aa6,"profile6_floresta.csv", row.names = FALSE)
write.csv(aa7,"profile7_floresta.csv", row.names = FALSE)
write.csv(aa8,"profile8_floresta.csv", row.names = FALSE)
write.csv(aa9,"profile9_floresta.csv", row.names = FALSE)
write.csv(aa10,"profile10_floresta.csv", row.names = FALSE)
write.csv(pam.res$clustering, "cluster.csv", row.names = FALSE)
