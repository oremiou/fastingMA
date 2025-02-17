library(rio)
dat1=import("data/fasting.xlsx", sheet="satisfaction")

#forest(p1)

# MH without cc
library(meta)
m1=metacont(mean.e = mean.E, sd.e=SD.E, n.e=N.E, mean.c = mean.C, sd.c=SD.C, n.c=N.C, 
           studlab=study, 
            data=dat1, MH.exact = T, prediction=F  ,
           sm="SMD")
forest(m1, label.left="favors fasting", label.right = "favors control", 
       main="Composite outcome", header=T, xlim=c(-2,2))

