library(rio)
dat1=import("data/fasting2.xlsx", sheet="COMPOSITE")
dat1=dat1[dat1$study %in% c("Mishra et al",  "Woods et al", 
                            "Boukantar et al", "Tamborrino et al", "Ferreira et al", "Mitchell et al"),]

# estimate baseline event rate
p1=metaprop(event=events.C, n=total.C, data=dat1, method="glm")
p2=metaprop(event=events.E, n=total.E, data=dat1, method="glm")
forest(p1)

# MH without cc
library(meta)
m1=metabin(event.c = events.E, n.c=total.E, event.e = events.C, n.e=total.C, 
           studlab=study, 
           method = "MH", data=dat1, MH.exact = T, prediction=T  ,
           sm="OR")
forest(m1, label.left="favors fasting", label.right = "favors non-fasting", 
       main="Composite outcome", header=T, xlim=c(0.1,5), 
       label.e="Fasting",label.c="Non-fasting", prediction=T)


# Mantel-Haenszel odds-ratio with 0.5 continuity correction
m2=metabin(event.c = events.E, n.c=total.E, event.e = events.C, n.e=total.C,
           studlab=study, MH.exact = F,
           method = "inverse", data=dat1 ,prediction=T  ,
           sm="OR", incr=0.5)
forest(m2,  label.left="favors fasting", label.right = "favors non-fasting", 
       main="Composite outcome", header=T, xlim=c(0.1,5), 
       label.e="Fasting",label.c="Non-fasting", prediction=T)


# Mantel-Haenszel risk-difference with no continuity correction
MH.RD<- metabin(event.c = events.E, n.c=total.E, event.e = events.C, n.e=total.C,
                studlab=study, MH.exact = T, data=dat1 ,prediction=T  ,
                sm="RD")
forest(MH.RD,  label.left="favors fasting", label.right = "favors non-fasting", 
       main="Composite outcome", header=T, xlim=c(-0.1,0.1), 
       label.e="Fasting",label.c="Non-fasting", prediction=T)


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
data <- list(NS=length(dat1$study), events.treat=dat1$events.C, events.control=dat1$events.E, 
             total.treat=dat1$total.C, total.control=dat1$total.E)
jags.m=0
n.chains =6
jags.m <- jags.model(model1.spec, data = data, n.chains =n.chains, n.adapt = 5000)

params <- c("or", "tauOR.s") 
closeAllConnections()
sampsRE<- coda.samples(jags.m, params, n.iter =100000)

library(MCMCvis)
MCMCsummary(sampsRE)
MCMCtrace(sampsRE,pdf = FALSE, params = "or") 


