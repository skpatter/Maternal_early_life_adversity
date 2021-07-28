##################
##################
##   Mortality  ##
##################
##################

d <- read.csv("Mortality 2yrs.csv", header = TRUE, na.strings = "")

### parity
d$multiparous <- ifelse(d$Parity=="M" , 1 , 0 )
d$primiparous <- ifelse(d$Parity=="P" , 1 , 0 )
d <- subset(d, multiparous==1)

d$n_grpsize <- (d$group.size.at.birth - min(d$group.size.at.birth))/(max(d$group.size.at.birth)-min(d$group.size.at.birth))
d$n_drought <- (d$biomassRev - min(d$biomassRev))/(max(d$biomassRev)-min(d$biomassRev))
d$n_ibi <-(d$ibidays - min(d$ibidays))/(max(d$ibidays)-min(d$ibidays))
d$n_momlossadj <- (d$age.at.mom.death - min(d$age.at.mom.death))/(max(d$age.at.mom.death)-min(d$age.at.mom.death))

d$s_ageopuntia <- (d$age.at.opuntia.benefit - mean(d$age.at.opuntia.benefit))/sd(d$age.at.opuntia.benefit)
d$s_rank <- (d$Rank_birth - mean(d$Rank_birth))/sd(d$Rank_birth)

d$adversity.cumul <- d$n_grpsize + d$n_drought + d$n_ibi + d$n_momlossadj + d$Primiparous
d$s_adversity.cumul <- (d$adversity.cumul - mean(d$adversity.cumul))/sd(d$adversity.cumul)
d$s_adversity.binary <- (d$BinaryIndexAll- mean(d$BinaryIndexAll))/sd(d$BinaryIndexAll)

d$s_grpsizeBirth <- (d$Grpsize_birth - mean(d$Grpsize_birth))/sd(d$Grpsize_birth)
d$s_infibi <- (d$IBI_full - mean(d$IBI_full))/sd(d$IBI_full)

d$MomAge <- (d$infdob - d$DOBmom)/365
d$s_MomAge <- (d$MomAge- mean(d$MomAge))/sd(d$MomAge)

d$Mom <- as.factor(d$Mom)
d$Mom <- as.integer(as.factor(d$Mom))

#Mortality model 1 (rank and ELA)
mortality2yrs1 <- map2stan( 
  alist(
    
    died ~ dbinom( 1 , p ) ,
    
    logit(p) <- ap + ap_Focal[Focal] + bp_ageop*ageop + bp_momage*momage +
      bp_grpsize*grpsize + bp_rank*relrank +
      bp_adv*advNorm  ,
    
    
    c(ap_Focal)[Focal] ~ dnorm(0,sigma_focal) ,
    c(ap,bp_ageop,bp_adv,
      bp_rank,bp_grpsize,bp_momage) ~ dnorm(0,1) ,
    c(sigma_focal)  ~ dcauchy(0,2)  
    
  ) , 
  data=list(
    died=d$Mortality,
    Focal=d$Mom,
    grpsize=d$s_grpsizeBirth,
    momage = d$s_MomAge,
    ageop = d$s_ageopuntia,
    advNorm = d$s_adversity.cumul ,
    relrank = d$s_rank
  ),
  cores=3 , chains=3 , warmup=1500, iter=8000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

#Mortality model 2 (rank X ela)
mortality2yrs2 <- map2stan( 
  alist(
    
    died ~ dbinom( 1 , p ) ,
    
    logit(p) <- ap + ap_Focal[Focal] + bp_ageop*ageop + bp_momage*momage +
      bp_grpsize*grpsize + bp_rank*relrank +
      bp_adv*advNorm  + bp_adv_rank*advNorm*relrank ,
    
    
    c(ap_Focal)[Focal] ~ dnorm(0,sigma_focal) ,
    c(ap,bp_ageop,bp_adv,
      bp_rank,bp_grpsize,bp_adv_rank,bp_momage) ~ dnorm(0,1) ,
    c(sigma_focal)  ~ dcauchy(0,2)  
    
  ) , 
  data=list(
    died=d$Mortality,
    Focal=d$Mom,
    grpsize=d$s_grpsizeBirth,
    momage = d$s_MomAge, 
    ageop = d$s_ageopuntia,
    advNorm = d$s_adversity.cumul ,
    relrank = d$s_rank
  ),
  cores=3 , chains=3 , warmup=1500, iter=8000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

# for mortality separate ela componenets models
d$n_grpsize <- (d$group.size.at.birth - min(d$group.size.at.birth))/(max(d$group.size.at.birth)-min(d$group.size.at.birth))
d$n_drought <- (d$biomass.year - min(d$biomass.year))/(max(d$biomass.year)-min(d$biomass.year))
d$n_ibi <-(d$ibidays - min(d$ibidays))/(max(d$ibidays)-min(d$ibidays))
d$n_momlossadj <- (d$age.at.mom.s.death.sep - min(d$age.at.mom.s.death.sep))/(max(d$age.at.mom.s.death.sep)-min(d$age.at.mom.s.death.sep))

d$s_grpsize <- (d$n_grpsize - mean(d$n_grpsize))/sd(d$n_grpsize)
d$s_drought <- (d$n_drought - mean(d$n_drought))/sd(d$n_drought)
d$s_ibi <- (d$n_ibi - mean(d$n_ibi))/sd(d$n_ibi)
d$s_momlossadj <- (d$n_momlossadj - mean(d$n_momlossadj))/sd(d$n_momlossadj)

dd <- d
d$s_grpsizeBirth <- (d$Grpsize_birth - mean(d$Grpsize_birth))/sd(d$Grpsize_birth)

d$Mom <- as.factor(d$Mom)
d$Mom <- as.integer(as.factor(d$Mom))

#mortality separate ela model 1
mortality2yrs1_separate <- map2stan( 
  alist(
    
    died ~ dbinom( 1 , p ) ,
    
    logit(p) <- ap + ap_Focal[Focal] + bp_ageop*ageop + bp_momage*momage + 
      bp_grpsize*grpsize + bp_rank*relrank +
      bp_grpsizeELA*grpsizeELA + bp_drought*drought + 
      bp_ibiELA*ibiELA + bp_momloss*momloss + bp_firstborn*first  ,
    
    
    c(ap_Focal)[Focal] ~ dnorm(0,sigma_focal) ,
    c(ap,bp_ageop,bp_grpsizeELA,bp_drought,
      bp_ibiELA,bp_momloss,
      bp_rank,bp_grpsize,bp_firstborn,bp_momage) ~ dnorm(0,1) ,
    c(sigma_focal)  ~ dcauchy(0,2)  
    
  ) , 
  data=list(
    died=d$Mortality,
    Focal=d$Mom,
    grpsize=d$s_grpsizeBirth,
    ageop = d$s_ageopuntia,
    grpsizeELA = d$s_grpsize,
    drought = d$s_drought,
    ibiELA = d$s_ibi,
    momloss = d$s_momlossadj,
    relrank = d$s_rank,
    first = d$Primiparous,
    momage = d$s_MomAge
  ),
  cores=3 , chains=3 , warmup=1500, iter=10000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

#mortality separate ela model 2
mortality2yrs2_separate <- map2stan( 
  alist(
    
    died ~ dbinom( 1 , p ) ,
    
    logit(p) <- ap + ap_Focal[Focal] + bp_ageop*ageop + bp_momage*momage +
      bp_grpsize*grpsize + bp_rank*relrank +
      bp_grpsizeELA*grpsizeELA + bp_drought*drought + 
      bp_ibiELA*ibiELA + bp_momloss*momloss  + bp_firstborn*first +
      (bp_rank_grpsizeELA*grpsizeELA + bp_rank_drought*drought + 
         bp_rank_ibiELA*ibiELA + bp_rank_momloss*momloss + bp_rank_firstborn*first)*relrank ,
    
    
    c(ap_Focal)[Focal] ~ dnorm(0,sigma_focal) ,
    c(ap,bp_ageop,bp_grpsizeELA,bp_drought,
      bp_ibiELA,bp_momloss,bp_firstborn,
      bp_rank,bp_grpsize,bp_rank_grpsizeELA,bp_rank_drought,
      bp_rank_ibiELA,bp_rank_momloss,bp_rank_firstborn,bp_momage) ~ dnorm(0,1) ,
    c(sigma_focal)  ~ dcauchy(0,2)  
    
  ) , 
  data=list(
    died=d$Mortality,
    Focal=d$Mom,
    grpsize=d$s_grpsizeBirth,
    ageop = d$s_ageopuntia,
    grpsizeELA = d$s_grpsize,
    drought = d$s_drought,
    ibiELA = d$s_ibi,
    momloss = d$s_momlossadj,
    relrank = d$s_rank,
    first = d$Primiparous,
    momage = d$s_MomAge
  ),
  cores=3 , chains=3 , warmup=1500, iter=8000, WAIC=TRUE, types=list(adapt.delta=0.99)
)


#frequentist models

library(coxme)

full_model1<-coxme(Surv(days_survived,died)~s_cumadvprimip+s_ageopuntia2 + s_rank +
                    s_MomAge + s_grpsizeBirth+(1|Mom),data=dd) 
summary(full_model1)

full_model1_interaction<-coxme(Surv(days_survived,died)~s_cumadvprimip+s_ageopuntia2 + s_rank +
                     s_MomAge + s_grpsizeBirth + s_cumadvprimip*s_rank + (1|Mom),data=dd) 
summary(full_model1_interaction)



##################
###            ###
###   Plots    ###
###            ###
##################



## Mortality plots
a_foc_z <- matrix(0,1000,length(unique(d$Mom)))
age.seq=seq(min(d$s_adversity.cumul),max(d$s_adversity.cumul),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  grpsize=rep(mean(d$s_grpsizeBirth),length(age.seq)),
  momage=rep(mean(d$s_MomAge),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  relrank=rep(mean(d$s_rank),length(age.seq)),
  advNorm=age.seq
)

LM1a <- link(mortality2yrs1, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2a <- link(mortality2yrs2, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)

w <- compare(mortality2yrs1,mortality2yrs2 , sort=FALSE)@output$weight  ##extract weights from table
w <- w/sum(w)  ##should add up to 1
idw <- round( w * 1000 )

#generate predictions from each model that will be averaged
PM1a <- LM1a
PM2a <- LM2a

pred_1a <- PM1a
pred_1a[1:idw[2],] <- PM2a[1:idw[2],] #overwrite first pred with weighted predicitons from pm2a

pred1a <- pred_1a[sample(nrow(pred_1a)),] #randomly mixes
pred1a.median <- apply(pred1a, 2, median )
pred1a.HPDI <-apply(pred1a, 2, HPDI )

par(mfrow = c(1, 1), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))

xxa=age.seq
yya=as.vector(t(pred1a[1:1000,]))
smoothScatter(rep(xxa,1000),yya,xlim=c(min(d$s_adversity.cumul),max(d$s_adversity.cumul)), colramp=colorRampPalette(c("white","#33CCFF")),nbin=200,transformation = function(x) x^.5,ylab='',cex=1.2,xaxt='n',ylim=c(0,1),nrpoints=0)
points(Mortality ~ s_adversity.cumul, data=dd , col=alpha("#33CCFF",0.6),pch=16, cex=0.9)

lines( age.seq , pred1a.median ,lwd=1,col="black")
lines(age.seq, pred1a.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1a.HPDI[2,],lty=2,lwd=1)

lab.a=seq(1,2.5,.5)
lab.sa=(lab.a-mean(d$adversity.cumul))/sd(d$adversity.cumul)
axis(1,labels=NA, at=lab.sa, tck=-0.01)
mtext(lab.a,at=lab.sa,side=1, line=.5)

mtext( "Cumulative early life adversity score" , side=1 , line=2, outer=FALSE , cex=1.5)
mtext( "Probability of mortality" , side=2 , line=2, outer=FALSE , cex=1.5)


#Mortality interaction with rank
a_foc_z <- matrix(0,1000,length(unique(d$Mom)))
age.seq=seq(min(d$s_adversity.cumul),max(d$s_adversity.cumul),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  grpsize=rep(mean(d$s_grpsizeBirth),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  relrank=rep(-1,length(age.seq)),
  momage=rep(mean(d$s_MomAge),length(age.seq)),
  advNorm=age.seq
)

LM1a <- link(mortality2yrs1, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2a <- link(mortality2yrs2, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)

w <- compare(mortality2yrs1,mortality2yrs2 , sort=FALSE)@output$weight
w <- w/sum(w)
idw <- round( w * 1000 )

PM1a <- LM1a
PM2a <- LM2a


pred_1a <- PM1a
pred_1a[1:idw[2],] <- PM2a[1:idw[2],]

pred1a <- pred_1a[sample(nrow(pred_1a)),]
pred1a.median <- apply(pred1a, 2, median )
pred1a.HPDI <-apply(pred1a, 2, HPDI )

d.pred2 <- list(
  Focal=rep(1,length(age.seq)),
  grpsize=rep(mean(d$s_grpsizeBirth),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  relrank=rep(1,length(age.seq)),
  momage=rep(mean(d$s_MomAge),length(age.seq)),
  advNorm=age.seq
)

LM1b <- link(mortality2yrs1, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2b <- link(mortality2yrs2, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)

PM1b <- LM1b
PM2b <- LM2b


pred_1b <- PM1b
pred_1b[1:idw[2],] <- PM2b[1:idw[2],]

pred1b <- pred_1b[sample(nrow(pred_1b)),] 
pred1b.median <- apply(pred1b, 2, median )
pred1b.HPDI <-apply(pred1b, 2, HPDI )

par(mfrow = c(1, 2), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))

plot( Mortality[d$s_rank==-1] ~ s_adversity.cumul[d$s_rank==-1] , data=dd , col=alpha("red",0.5), pch=16 ,ylim=c(0,1),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_adversity.cumul),max(d$s_adversity.cumul)))
pred1a.lines=pred1a[sample(nrow(pred1a),100,replace=F),]
for (i in 1:100){
  lines( age.seq , pred1a.lines[i,] ,lwd=2,col=alpha("red",0.1))
}
lines( age.seq , pred1a.median ,lwd=1,col="black")
lines(age.seq, pred1a.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1a.HPDI[2,],lty=2,lwd=1)
axis(2,cex.axis=.9)
lab.a=seq(1,2.5,.5)
lab.sa=(lab.a-mean(d$adversity.cumul))/sd(d$adversity.cumul)
axis(1,labels=NA, at=lab.sa, tck=-0.01)
mtext(lab.a,at=lab.sa,side=1, line=.5)
mtext( "Cumulative early life adversity" , side=1 , line=1.5, outer=FALSE , cex=1.3)
mtext( "Low rank" , side=3 , line=-1, outer=FALSE , cex=1.5)

plot( Mortality[d$s_rank==1] ~ s_adversity.cumul[d$s_rank==1] , data=dd , col=alpha("#33CCFF",0.5), pch=16 ,ylim=c(0,1),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_adversity.cumul),max(d$s_adversity.cumul)))
pred1b.lines=pred1b[sample(nrow(pred1b),100,replace=F),]

for (i in 1:100){
  lines( age.seq , pred1b.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.1))
}

lines( age.seq , pred1b.median ,lwd=1,col="black")
lines(age.seq, pred1b.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1b.HPDI[2,],lty=2,lwd=1)

lab.a=seq(1,2.5,.5)
lab.sa=(lab.a-mean(d$adversity.cumul))/sd(d$adversity.cumul)
axis(1,labels=NA, at=lab.sa, tck=-0.01)
mtext(lab.a,at=lab.sa,side=1, line=.5)
mtext( "Cumulative early life adversity" , side=1 , line=1.5, outer=FALSE , cex=1.3)
mtext( "High rank" , side=3 , line=-1, outer=FALSE , cex=1.5)
mtext("Probability of mortality", side=2, line=2,cex=1.5, outer=TRUE)



