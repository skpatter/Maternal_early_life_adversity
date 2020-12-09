##################
##################
##   Mortality  ##
##################
##################

d <- read.csv("Mortality 2yrs.csv", header = TRUE, na.strings = "")
o <- read.csv("advscores_short.csv", header = TRUE, na.strings = "")

merge <- merge(d,o)
d <- merge

### parity
d$multiparous <- ifelse(d$Parity=="M" , 1 , 0 )
d$primiparous <- ifelse(d$Parity=="P" , 1 , 0 )
d <- subset(d, multiparous==1)

d$n_grpsize <- (d$group.size.at.birth - min(d$group.size.at.birth))/(max(d$group.size.at.birth)-min(d$group.size.at.birth))
d$n_drought <- (d$biomassRev - min(d$biomassRev))/(max(d$biomassRev)-min(d$biomassRev))
d$n_ibi <-(d$ibidays - min(d$ibidays))/(max(d$ibidays)-min(d$ibidays))
d$n_momlossadj <- (d$age.at.mom.death - min(d$age.at.mom.death))/(max(d$age.at.mom.death)-min(d$age.at.mom.death))

d$s_ageopuntia2 <- (d$age.at.opuntia.benefit - mean(d$age.at.opuntia.benefit))/sd(d$age.at.opuntia.benefit)
d$s_rank <- (d$Rank_birth - mean(d$Rank_birth))/sd(d$Rank_birth)

dd <- d
dd$adversity.cumul <- dd$n_grpsize + dd$n_drought + dd$n_ibi + dd$n_momlossadj + dd$Primiparous #2 is momlossadj and original is momloss
dd$s_adversity.cumul <- (dd$adversity.cumul - mean(dd$adversity.cumul))/sd(dd$adversity.cumul)
dd$s_adversity.binary <- (dd$BinaryIndexAll- mean(dd$BinaryIndexAll))/sd(dd$BinaryIndexAll)

dd$s_grpsizeBirth <- (dd$Grpsize_birth - mean(dd$Grpsize_birth))/sd(dd$Grpsize_birth)
dd$s_infibi <- (dd$IBI_full - mean(dd$IBI_full))/sd(dd$IBI_full)

dd$MomAge <- (dd$infdob - dd$DOBmom)/365
dd$s_MomAge <- (dd$MomAge- mean(dd$MomAge))/sd(dd$MomAge)

dd$Mom <- as.factor(dd$Mom)
dd$Mom <- as.integer(as.factor(dd$Mom))

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
    died=dd$Mortality,
    Focal=dd$Mom,
    grpsize=dd$s_grpsizeBirth,
    momage = dd$s_MomAge,
    ageop = dd$s_ageopuntia2,
    advNorm = dd$s_adversity.cumul ,
    relrank = dd$s_rank
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
    died=dd$Mortality,
    Focal=dd$Mom,
    grpsize=dd$s_grpsizeBirth,
    momage = dd$s_MomAge, 
    ageop = dd$s_ageopuntia2,
    advNorm = dd$s_adversity.cumul ,
    relrank = dd$s_rank
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
dd$s_grpsizeBirth <- (dd$Grpsize_birth - mean(dd$Grpsize_birth))/sd(dd$Grpsize_birth)

dd$Mom <- as.factor(dd$Mom)
dd$Mom <- as.integer(as.factor(dd$Mom))

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
    died=dd$Mortality,
    Focal=dd$Mom,
    grpsize=dd$s_grpsizeBirth,
    ageop = dd$s_ageopuntia2,
    grpsizeELA = dd$s_grpsize,
    drought = dd$s_drought,
    ibiELA = dd$s_ibi,
    momloss = dd$s_momlossadj,
    relrank = dd$s_rank,
    first = dd$Primiparous,
    momage = dd$s_MomAge
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
    died=dd$Mortality,
    Focal=dd$Mom,
    grpsize=dd$s_grpsizeBirth,
    ageop = dd$s_ageopuntia,
    grpsizeELA = dd$s_grpsize,
    drought = dd$s_drought,
    ibiELA = dd$s_ibi,
    momloss = dd$s_momlossadj,
    relrank = dd$s_rank,
    first = dd$Primiparous,
    momage = dd$s_MomAge
  ),
  cores=3 , chains=3 , warmup=1500, iter=8000, WAIC=TRUE, types=list(adapt.delta=0.99)
)


##################
###            ###
###   Plots    ###
###            ###
##################



## Mortality plots
a_foc_z <- matrix(0,1000,length(unique(dd$Mom)))
age.seq=seq(min(dd$s_adversity.cumul),max(dd$s_adversity.cumul),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  grpsize=rep(mean(dd$s_grpsizeBirth),length(age.seq)),
  momage=rep(mean(dd$s_MomAge),length(age.seq)),
  ageop=rep(mean(dd$s_ageopuntia),length(age.seq)),
  relrank=rep(mean(dd$s_rank),length(age.seq)),
  advNorm=age.seq
)

LM1a <- link(mortality2yrs1, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2a <- link(mortality2yrs2, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)

w <- compare(mortality2yrs1,mortality2yrs2 , sort=FALSE)@output$weight  ##extract weights from table, alphabetically ordered, important to keep models in alph ordered so proper weights are attributed
w <- w/sum(w)  ##make sure everything adds up to 1
idw <- round( w * 1000 ) ###round weights to nearest integer so samples we extract sum up to 1000

#generate predictions from each model that will be averaged
PM1a <- LM1a #for bndur_M32
PM2a <- LM2a #for bndur_M3adv

#make a vector or predictions sampled according to model weight, should be ordereed alphabetically/numerically
pred_1a <- PM1a #makes pred nurse all predictions from bndur_M32
pred_1a[1:idw[2],] <- PM2a[1:idw[2],] #overwrite predfather with approriate number of weighted predicitons from PF3a

pred1a <- pred_1a[sample(nrow(pred_1a)),] #randomly mixes vector-- important to make graphing individual model predictions random
pred1a.median <- apply(pred1a, 2, median )
pred1a.HPDI <-apply(pred1a, 2, HPDI )

par(mfrow = c(1, 1), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))

xxa=age.seq
yya=as.vector(t(pred1a[1:1000,]))
smoothScatter(rep(xxa,1000),yya,xlim=c(min(dd$s_adversity.cumul),max(dd$s_adversity.cumul)), colramp=colorRampPalette(c("white","#33CCFF")),nbin=200,transformation = function(x) x^.5,ylab='',cex=1.2,xaxt='n',ylim=c(0,1),nrpoints=0)
points(Mortality ~ s_adversity.cumul, data=dd , col=alpha("#33CCFF",0.6),pch=16, cex=0.9)

lines( age.seq , pred1a.median ,lwd=1,col="black")
lines(age.seq, pred1a.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1a.HPDI[2,],lty=2,lwd=1)

lab.a=seq(1,2.5,.5)
lab.sa=(lab.a-mean(dd$adversity.cumul))/sd(dd$adversity.cumul)
axis(1,labels=NA, at=lab.sa, tck=-0.01)
mtext(lab.a,at=lab.sa,side=1, line=.5)

mtext( "Cumulative early life adversity score" , side=1 , line=2, outer=FALSE , cex=1.5)
mtext( "Probability of mortality" , side=2 , line=2, outer=FALSE , cex=1.5)


#Mortality interaction with rank
a_foc_z <- matrix(0,1000,length(unique(dd$Mom)))
age.seq=seq(min(dd$s_adversity.cumul),max(dd$s_adversity.cumul),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  grpsize=rep(mean(dd$s_grpsizeBirth),length(age.seq)),
  ageop=rep(mean(dd$s_ageopuntia),length(age.seq)),
  relrank=rep(-1,length(age.seq)),
  momage=rep(mean(dd$s_MomAge),length(age.seq)),
  advNorm=age.seq
)

LM1a <- link(mortality2yrs1, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2a <- link(mortality2yrs2, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)

w <- compare(mortality2yrs1,mortality2yrs2 , sort=FALSE)@output$weight  ##extract weights from table, alphabetically ordered, important to keep models in alph ordered so proper weights are attributed
w <- w/sum(w)  ##make sure everything adds up to 1
idw <- round( w * 1000 ) ###round weights to nearest integer so samples we extract sum up to 1000

#generate predictions from each model that will be averaged
PM1a <- LM1a #for bndur_M32
PM2a <- LM2a #for bndur_M3adv

#make a vector or predictions sampled according to model weight, should be ordereed alphabetically/numerically
pred_1a <- PM1a #makes pred nurse all predictions from bndur_M32
pred_1a[1:idw[2],] <- PM2a[1:idw[2],] #overwrite predfather with approriate number of weighted predicitons from PF3a

pred1a <- pred_1a[sample(nrow(pred_1a)),] #randomly mixes vector-- important to make graphing individual model predictions random
pred1a.median <- apply(pred1a, 2, median )
pred1a.HPDI <-apply(pred1a, 2, HPDI )

d.pred2 <- list(
  Focal=rep(1,length(age.seq)),
  grpsize=rep(mean(dd$s_grpsizeBirth),length(age.seq)),
  ageop=rep(mean(dd$s_ageopuntia),length(age.seq)),
  relrank=rep(1,length(age.seq)),
  momage=rep(mean(dd$s_MomAge),length(age.seq)),
  advNorm=age.seq
)

LM1b <- link(mortality2yrs1, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2b <- link(mortality2yrs2, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)

#generate predictions from each model that will be averaged
PM1b <- LM1b #for bndur_M32
PM2b <- LM2b #for bndur_M3adv

#make a vector or predictions sampled according to model weight, should be ordereed alphabetically/numerically
pred_1b <- PM1b #makes pred nurse all predictions from bndur_M32
pred_1b[1:idw[2],] <- PM2b[1:idw[2],] #overwrite predfather with approriate number of weighted predicitons from PF3a

pred1b <- pred_1b[sample(nrow(pred_1b)),] #randomly mixes vector-- important to make graphing individual model predictions random
pred1b.median <- apply(pred1b, 2, median )
pred1b.HPDI <-apply(pred1b, 2, HPDI )

par(mfrow = c(1, 2), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))

plot( Mortality[dd$s_rank==-1] ~ s_adversity.cumul[dd$s_rank==-1] , data=dd , col=alpha("red",0.5), pch=16 ,ylim=c(0,1),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(dd$s_adversity.cumul),max(dd$s_adversity.cumul)))
pred1a.lines=pred1a[sample(nrow(pred1a),100,replace=F),]
for (i in 1:100){
  lines( age.seq , pred1a.lines[i,] ,lwd=2,col=alpha("red",0.1))
}
lines( age.seq , pred1a.median ,lwd=1,col="black")
lines(age.seq, pred1a.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1a.HPDI[2,],lty=2,lwd=1)
axis(2,cex.axis=.9)
lab.a=seq(1,2.5,.5)
lab.sa=(lab.a-mean(dd$adversity.cumul))/sd(dd$adversity.cumul)
axis(1,labels=NA, at=lab.sa, tck=-0.01)
mtext(lab.a,at=lab.sa,side=1, line=.5)
mtext( "Cumulative early life adversity" , side=1 , line=1.5, outer=FALSE , cex=1.3)
mtext( "Low rank" , side=3 , line=-1, outer=FALSE , cex=1.5)

plot( Mortality[dd$s_rank==1] ~ s_adversity.cumul[dd$s_rank==1] , data=dd , col=alpha("#33CCFF",0.5), pch=16 ,ylim=c(0,1),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(dd$s_adversity.cumul),max(dd$s_adversity.cumul)))
pred1b.lines=pred1b[sample(nrow(pred1b),100,replace=F),]

for (i in 1:100){
  lines( age.seq , pred1b.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.1))
}

lines( age.seq , pred1b.median ,lwd=1,col="black")
lines(age.seq, pred1b.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1b.HPDI[2,],lty=2,lwd=1)

lab.a=seq(1,2.5,.5)
lab.sa=(lab.a-mean(dd$adversity.cumul))/sd(dd$adversity.cumul)
axis(1,labels=NA, at=lab.sa, tck=-0.01)
mtext(lab.a,at=lab.sa,side=1, line=.5)
mtext( "Cumulative early life adversity" , side=1 , line=1.5, outer=FALSE , cex=1.3)
mtext( "High rank" , side=3 , line=-1, outer=FALSE , cex=1.5)
mtext("Probability of mortality", side=2, line=2,cex=1.5, outer=TRUE)



