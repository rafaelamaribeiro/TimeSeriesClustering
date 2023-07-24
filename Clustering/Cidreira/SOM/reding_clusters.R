rm(list=ls())

library(lubridate)

#setwd("/home/rafaela/Documents/PUC/2022.1/TFC/codigos/cidreira")
setwd("C:\\Users\\Rafaela\\Downloads\\2022\\TFC\\codigos eol\\Cidreira")

que_bola <- read.csv("clust_som.csv")

head(que_bola)
max(que_bola$x)
min(que_bola$x)

que_cluster <- read.csv("12k_som.csv")

unique(que_cluster$x)

rows <- length(que_bola$X)
clust1 <- c()
clust2 <- c()
clust3 <- c()
clust4 <- c()
clust5 <- c() 
clust6 <- c() 
clust7 <- c() 
clust8 <- c() 
clust9 <- c() 
clust10 <- c() 
clust11 <- c() 
clust12 <- c() 
m <- 1
n <- 1
o <- 1
p <- 1
q <- 1
s <- 1
t <- 1
u <- 1
v <- 1
w <- 1
x <- 1
y <- 1
for (r in 1:rows) {
  bola = que_bola$x[r]
  if (que_cluster$x[bola] == 1) {
    clust1[m] <- r;
    m <- m + 1
  }
  else if (que_cluster$x[bola] == 2) {
    clust2[n] <- r;
    n <- n + 1
  }
  else if (que_cluster$x[bola] == 3) {
    clust3[o] <- r;
    o <- o + 1
  }
  else if (que_cluster$x[bola] == 4) {
    clust4[p] <- r;
    p <- p + 1
  }
  else if (que_cluster$x[bola] == 5) {
    clust5[q] <- r;
    q <- q + 1
  }
  else if (que_cluster$x[bola] == 6) {
    clust6[s] <- r;
    s <- s + 1
  }
  else if (que_cluster$x[bola] == 7) {
    clust7[t] <- r;
    t <- t + 1
  }
  else if (que_cluster$x[bola] == 8) {
    clust8[u] <- r;
    u <- u + 1
  }
  else if (que_cluster$x[bola] == 9) {
    clust9[v] <- r;
    v <- v + 1
  }
  else if (que_cluster$x[bola] == 10) {
    clust10[w] <- r;
    w <- w + 1
  }
  else if (que_cluster$x[bola] == 11) {
    clust11[x] <- r;
    x <- x + 1
  }
  else if (que_cluster$x[bola] == 12) {
    clust12[y] <- r;
    y <- y + 1
  }
}

write.csv(clust1, "cidreira_clust1_som_12K_16n.csv")
write.csv(clust2, "cidreira_clust2_som_12K_16n.csv")
write.csv(clust3, "cidreira_clust3_som_12K_16n.csv")
write.csv(clust4, "cidreira_clust4_som_12K_16n.csv")
write.csv(clust5, "cidreira_clust5_som_12K_16n.csv")
write.csv(clust6, "cidreira_clust6_som_12K_16n.csv")
write.csv(clust7, "cidreira_clust7_som_12K_16n.csv")
write.csv(clust8, "cidreira_clust8_som_12K_16n.csv")
write.csv(clust9, "cidreira_clust9_som_12K_16n.csv")
write.csv(clust10, "cidreira_clust10_som_12K_16n.csv")
write.csv(clust11, "cidreira_clust11_som_12K_16n.csv")
write.csv(clust12, "cidreira_clust12_som_12K_16n.csv")

length(clust1)
length(clust2)
length(clust3)

cidreira <- read.csv2("hweol_CIDREIRA_RS.csv")
#cidreira <- cidreira[1:359160,]
head(cidreira)
tail(cidreira)

cidreira$dates_new <- as.Date(cidreira$dates, format = "%Y-%m-%d %H:%M")

dates <- unique(cidreira$dates_new)
tail(dates)
meses <- month(dates)
length(meses)

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
    if (i %in% clust1) {
      jan1 = jan1 + 1}
    else if (i %in% clust2) {
      jan2 = jan2 + 1}
    else if (i %in% clust3) {
      jan3 = jan3 + 1}
    else if (i %in% clust4) {
      jan4 = jan4 + 1}
    else if (i %in% clust5) {
      jan5 = jan5 + 1}
    else if (i %in% clust6) {
      jan6 = jan6 + 1}
    else if (i %in% clust7) {
      jan7 = jan7 + 1}
    else if (i %in% clust8) {
      jan8 = jan8 + 1}
    else if (i %in% clust9) {
      jan9 = jan9 + 1}
    else if (i %in% clust10) {
      jan10 = jan10 + 1}
    else if (i %in% clust11) {
      jan11 = jan11 + 1}
    else if (i %in% clust12) {
      jan12 = jan12 + 1}}
  else if (meses[i] == 2) {
    if (i %in% clust1) {
      feb1 = feb1 + 1}
    else if (i %in% clust2) {
      feb2 = feb2 + 1}
    else if (i %in% clust3) {
      feb3 = feb3 + 1}
    else if (i %in% clust4) {
      feb4 = feb4 + 1}
    else if (i %in% clust5) {
      feb5 = feb5 + 1}
    else if (i %in% clust6) {
      feb6 = feb6 + 1}
    else if (i %in% clust7) {
      feb7 = feb7 + 1}
    else if (i %in% clust8) {
      feb8 = feb8 + 1}
    else if (i %in% clust9) {
      feb9 = feb9 + 1}
    else if (i %in% clust10) {
      feb10 = feb10 + 1}
    else if (i %in% clust11) {
      feb11 = feb11 + 1}
    else if (i %in% clust12) {
      feb12 = feb12 + 1}}
  else if (meses[i] == 3) {
    if (i %in% clust1) {
      mar1 = mar1 + 1}
    else if (i %in% clust2) {
      mar2 = mar2 + 1}
    else if (i %in% clust3) {
      mar3 = mar3 + 1}
    else if (i %in% clust4) {
      mar4 = mar4 + 1}
    else if (i %in% clust5) {
      mar5 = mar5 + 1}
    else if (i %in% clust6) {
      mar6 = mar6 + 1}
    else if (i %in% clust7) {
      mar7 = mar7 + 1}
    else if (i %in% clust8) {
      mar8 = mar8 + 1}
    else if (i %in% clust9) {
      mar9 = mar9 + 1}
    else if (i %in% clust10) {
      mar10 = mar10 + 1}
    else if (i %in% clust11) {
      mar11 = mar11 + 1}
    else if (i %in% clust12) {
      mar12 = mar12 + 1}}
  else if (meses[i] == 4) {
    if (i %in% clust1) {
      apr1 = apr1 + 1}
    else if (i %in% clust2) {
      apr2 = apr2 + 1}
    else if (i %in% clust3) {
      apr3 = apr3 + 1}
    else if (i %in% clust4) {
      apr4 = apr4 + 1}
    else if (i %in% clust5) {
      apr5 = apr5 + 1}
    else if (i %in% clust6) {
      apr6 = apr6 + 1}
    else if (i %in% clust7) {
      apr7 = apr7 + 1}
    else if (i %in% clust8) {
      apr8 = apr8 + 1}
    else if (i %in% clust9) {
      apr9 = apr9 + 1}
    else if (i %in% clust10) {
      apr10 = apr10 + 1}
    else if (i %in% clust11) {
      apr11 = apr11 + 1}
    else if (i %in% clust12) {
      apr12 = apr12 + 1}}
  else if (meses[i] == 5) {
    if (i %in% clust1) {
      may1 = may1 + 1}
    else if (i %in% clust2) {
      may2 = may2 + 1}
    else if (i %in% clust3) {
      may3 = may3 + 1}
    else if (i %in% clust4) {
      may4 = may4 + 1}
    else if (i %in% clust5) {
      may5 = may5 + 1}
    else if (i %in% clust6) {
      may6 = may6 + 1}
    else if (i %in% clust7) {
      may7 = may7 + 1}
    else if (i %in% clust8) {
      may8 = may8 + 1}
    else if (i %in% clust9) {
      may9 = may9 + 1}
    else if (i %in% clust10) {
      may10 = may10 + 1}
    else if (i %in% clust11) {
      may11 = may11 + 1}
    else if (i %in% clust12) {
      may12 = may12 + 1}}
  else if (meses[i] == 6) {
    if (i %in% clust1) {
      jun1 = jun1 + 1}
    else if (i %in% clust2) {
      jun2 = jun2 + 1}
    else if (i %in% clust3) {
      jun3 = jun3 + 1}
    else if (i %in% clust4) {
      jun4 = jun4 + 1}
    else if (i %in% clust5) {
      jun5 = jun5 + 1}
    else if (i %in% clust6) {
      jun6 = jun6 + 1}
    else if (i %in% clust7) {
      jun7 = jun7 + 1}
    else if (i %in% clust8) {
      jun8 = jun8 + 1}
    else if (i %in% clust9) {
      jun9 = jun9 + 1}
    else if (i %in% clust10) {
      jun10 = jun10 + 1}
    else if (i %in% clust11) {
      jun11 = jun11 + 1}
    else if (i %in% clust12) {
      jun12 = jun12 + 1}}
  else if (meses[i] == 7) {
    if (i %in% clust1) {
      jul1 = jul1 + 1}
    else if (i %in% clust2) {
      jul2 = jul2 + 1}
    else if (i %in% clust3) {
      jul3 = jul3 + 1}
    else if (i %in% clust4) {
      jul4 = jul4 + 1}
    else if (i %in% clust5) {
      jul5 = jul5 + 1}
    else if (i %in% clust6) {
      jul6 = jul6 + 1}
    else if (i %in% clust7) {
      jul7 = jul7 + 1}
    else if (i %in% clust8) {
      jul8 = jul8 + 1}
    else if (i %in% clust9) {
      jul9 = jul9 + 1}
    else if (i %in% clust10) {
      jul10 = jul10 + 1}
    else if (i %in% clust11) {
      jul11 = jul11 + 1}
    else if (i %in% clust12) {
      jul12 = jul12 + 1}}
  else if (meses[i] == 8) {
    if (i %in% clust1) {
      aug1 = aug1 + 1}
    else if (i %in% clust2) {
      aug2 = aug2 + 1}
    else if (i %in% clust3) {
      aug3 = aug3 + 1}
    else if (i %in% clust4) {
      aug4 = aug4 + 1}
    else if (i %in% clust5) {
      aug5 = aug5 + 1}
    else if (i %in% clust6) {
      aug6 = aug6 + 1}
    else if (i %in% clust7) {
      aug7 = aug7 + 1}
    else if (i %in% clust8) {
      aug8 = aug8 + 1}
    else if (i %in% clust9) {
      aug9 = aug9 + 1}
    else if (i %in% clust10) {
      aug10 = aug10 + 1}
    else if (i %in% clust11) {
      aug11 = aug11 + 1}
    else if (i %in% clust12) {
      aug12 = aug12 + 1}}
  else if (meses[i] == 9) {
    if (i %in% clust1) {
      sep1 = sep1 + 1}
    else if (i %in% clust2) {
      sep2 = sep2 + 1}
    else if (i %in% clust3) {
      sep3 = sep3 + 1}
    else if (i %in% clust4) {
      sep4 = sep4 + 1}
    else if (i %in% clust5) {
      sep5 = sep5 + 1}
    else if (i %in% clust6) {
      sep6 = sep6 + 1}
    else if (i %in% clust7) {
      sep7 = sep7 + 1}
    else if (i %in% clust8) {
      sep8 = sep8 + 1}
    else if (i %in% clust9) {
      sep9 = sep9 + 1}
    else if (i %in% clust10) {
      sep10 = sep10 + 1}
    else if (i %in% clust11) {
      sep11 = sep11 + 1}
    else if (i %in% clust12) {
      sep12 = sep12 + 1}}
  else if (meses[i] == 10) {
    if (i %in% clust1) {
      oct1 = oct1 + 1}
    else if (i %in% clust2) {
      oct2 = oct2 + 1}
    else if (i %in% clust3) {
      oct3 = oct3 + 1}
    else if (i %in% clust4) {
      oct4 = oct4 + 1}
    else if (i %in% clust5) {
      oct5 = oct5 + 1}
    else if (i %in% clust6) {
      oct6 = oct6 + 1}
    else if (i %in% clust7) {
      oct7 = oct7 + 1}
    else if (i %in% clust8) {
      oct8 = oct8 + 1}
    else if (i %in% clust9) {
      oct9 = oct9 + 1}
    else if (i %in% clust10) {
      oct10 = oct10 + 1}
    else if (i %in% clust11) {
      oct11 = oct11 + 1}
    else if (i %in% clust12) {
      oct12 = oct12 + 1}}
  else if (meses[i] == 11) {
    if (i %in% clust1) {
      nov1 = nov1 + 1}
    else if (i %in% clust2) {
      nov2 = nov2 + 1}
    else if (i %in% clust3) {
      nov3 = nov3 + 1}
    else if (i %in% clust4) {
      nov4 = nov4 + 1}
    else if (i %in% clust5) {
      nov5 = nov5 + 1}
    else if (i %in% clust6) {
      nov6 = nov6 + 1}
    else if (i %in% clust7) {
      nov7 = nov7 + 1}
    else if (i %in% clust8) {
      nov8 = nov8 + 1}
    else if (i %in% clust9) {
      nov9 = nov9 + 1}
    else if (i %in% clust10) {
      nov10 = nov10 + 1}
    else if (i %in% clust11) {
      nov11 = nov11 + 1}
    else if (i %in% clust12) {
      nov12 = nov12 + 1}}
  else if (meses[i] == 12) {
    if (i %in% clust1) {
      dec1 = dec1 + 1}
    else if (i %in% clust2) {
      dec2 = dec2 + 1}
    else if (i %in% clust3) {
      dec3 = dec3 + 1}
    else if (i %in% clust4) {
      dec4 = dec4 + 1}
    else if (i %in% clust5) {
      dec5 = dec5 + 1}
    else if (i %in% clust6) {
      dec6 = dec6 + 1}
    else if (i %in% clust7) {
      dec7 = dec7 + 1}
    else if (i %in% clust8) {
      dec8 = dec8 + 1}
    else if (i %in% clust9) {
      dec9 = dec9 + 1}
    else if (i %in% clust10) {
      dec10 = dec10 + 1}
    else if (i %in% clust11) {
      dec11 = dec11 + 1}
    else if (i %in% clust12) {
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

df

write.csv(df,"porMes_cidreira_som_12K_16n.csv", row.names = FALSE)

