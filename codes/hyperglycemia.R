library(rio)
dat1=import("data/fasting.xlsx", sheet="HYPERGLYCEMIA")

# estimate baseline event rate
p1=metaprop(event=events.C, n=total.C, data=dat1, method="glm")
p2=metaprop(event=events.E, n=total.E, data=dat1, method="glm")
#forest(p1)

# MH without cc
library(meta)
m1=metabin(event.e = events.E, n.e=total.E, event.c = events.C, n.c=total.C, 
           studlab=study, 
           method = "MH", data=dat1, MH.exact = T, prediction=F  ,
           sm="OR")
forest(m1, label.left="favors fasting", label.right = "favors control", 
       main="Composite outcome",prediction=T,label.e="Fasting", header=T, xlim=c(0.2,5))


# Mantel-Haenszel odds-ratio with 0.5 continuity correction
m2=metabin(event.e = events.E, n.e=total.E, event.c = events.C, n.c=total.C, 
           studlab=study, MH.exact = F,
           method = "inverse", data=dat1 ,prediction=F  ,
           sm="OR", incr=0.5)
forest(m2, label.left="favors fasting", label.right = "favors control", 
       main="Composite outcome",label.e="Fasting", header=T, xlim=c(0.2,5))


# Mantel-Haenszel risk-difference with no continuity correction
MH.RD<- metabin(event.e = events.E, n.e=total.E, event.c = events.C, n.c=total.C, 
                studlab=study, MH.exact = T, data=dat1 ,prediction=F,
                sm="RD")
forest(MH.RD, label.left="favors fasting", label.right = "favors control", 
       main="Composite outcome",prediction=T,label.e="Fasting", header=T, xlim=c(-0.1,.1))


#### Bayesian RE - MA ----------------------------------
library(R2jags)
model="
model {
for (i in 1:NS){ ### i cycles through the number of studies
events.control[i]~dbin(p.control[i],total.control[i])
events.treat[i]~dbin(p.treat[i],total.treat[i])
logit(p.control[i]) <- u[i]
logit(p.treat[i]) <- u[i]+lor[i]
u[i] ~dnorm(0,0.001)
lor[i]~dnorm(mu,prec.tau)
}
prec.tau<-1/tauOR.s
tauOR.s~dlnorm(-1.43,inv.sd2) ##prior based on empirical priors
inv.sd2<-1/(1.45*1.45) ##prior based on empirical priors
mu~dnorm(0,0.001)
or<-exp(mu)
}
"
model1.spec<-textConnection(model) 
data <- list(NS=length(dat1$study), events.treat=dat1$events.E, events.control=dat1$events.C, 
             total.treat=dat1$total.E, total.control=dat1$total.C)
jags.m=0
n.chains =4
jags.m <- jags.model(model1.spec, data = data, n.chains =n.chains, n.adapt = 500)

params <- c("or", "tauOR.s") 
closeAllConnections()
sampsRE<- coda.samples(jags.m, params, n.iter =10000)

library(MCMCvis)
MCMCsummary(sampsRE)
MCMCtrace(sampsRE,pdf = FALSE, params = "or") 


all.sampsRE=sampsRE[[1]]
for(i in 2:n.chains){ all.sampsRE=rbind(all.sampsRE, sampsRE[[i]])}
mean(all.sampsRE[,1]<1)

library(ggplot2)
all.sampsHRRE=data.frame(HR=all.sampsRE[,1])


