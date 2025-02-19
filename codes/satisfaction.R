library(rio)
dat1=import("data/fasting2.xlsx", sheet="satisfaction")
dat1$SD.C=as.numeric(dat1$SD.C)
#forest(p1)

# SMD
library(meta)
m1=metacont(mean.e = mean.C, sd.e=SD.C, n.e=N.C, mean.c = mean.E, sd.c=SD.E, n.c=N.E, 
           studlab=study, 
            data=dat1, MH.exact = T, prediction=T  ,
           sm="SMD")
forest(m1, label.left="favors fasting", label.right = "favors control", 
       main="Composite outcome", header=T, xlim=c(-2,2), 
       label.e="Fasting", label.c="Non-fasting")

# sensitivity analysis excluding study with imputed SD
dat2=dat1[dat1$study!="Atkinson et al",]
m2=metacont(mean.e = mean.C, sd.e=SD.C, n.e=N.C, mean.c = mean.E, sd.c=SD.E, n.c=N.E, 
            studlab=study, 
            data=dat2, MH.exact = T, prediction=T  ,
            sm="SMD")
forest(m2, label.left="favors fasting", label.right = "favors control", 
       main="Composite outcome", header=T, xlim=c(-2,2), 
       label.e="Fasting", label.c="Non-fasting")
