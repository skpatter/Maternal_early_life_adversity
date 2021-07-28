##################
###            ###
### GCM models ###
###            ###
##################

d <- read.csv("hormones ela.csv", header = TRUE, na.strings = "")

### parity
d$multiparous <- ifelse(d$parity=="M" , 1 , 0 )
d$primiparous <- ifelse(d$parity=="P" , 1 , 0 )
d <- subset(d, multiparous==1)
### groups
d$phg <- ifelse(d$troop=="phg" , 1 , 0 )
d$enk <- ifelse(d$troop=="enk" , 1 , 0 )
d$ynt <- ifelse(d$troop=="ynt" , 1 , 0 )
d$nmu <- ifelse(d$troop=="nmu" , 1 , 0 )
### infant sex
d$female <- ifelse(d$sex=="f" , 1 , 0 )
d$male <- ifelse(d$sex=="m" , 1 , 0 )
d$s_infAge <- (d$infageRel - mean(d$infageRel))/sd(d$infageRel)
d$s_infAgeSquared <- (d$s_infAge)^2

#current environment
d <- d[complete.cases(as.numeric(d$avgTotal)), ]
d$avgTotal <- as.numeric(d$avgTotal)
d$s_grass <- (d$avgTotal - mean(d$avgTotal))/sd(d$avgTotal)
d$Human_count <- as.numeric(d$Human_count)
d$s_human <- (d$Human_count - mean(d$Human_count))/sd(d$Human_count)
d$StrangeMales_count <- as.numeric(d$StrangeMales_count)
d$s_strangemales <- (d$StrangeMales_count - mean(d$StrangeMales_count))/sd(d$StrangeMales_count)
d$challenges <- d$Human_count + d$StrangeMales_count
d$s_challenges <- (d$challenges - mean(d$challenges))/sd(d$challenges)
d$s_grpsize_cur <- (d$group_size - mean(d$group_size))/sd(d$group_size)

d <- d[complete.cases(d$relRank), ]
d$relRank <- as.numeric(d$relRank)
d$s_rank <- (d$relRank - mean(d$relRank))/sd(d$relRank)

#early life adversity
d <- d[complete.cases(d$group.size.at.birth,d$biomassRev,d$ageatmomdeathRev,d$relRank), ]
d$n_grpsize <- (d$group.size.at.birth - min(d$group.size.at.birth))/(max(d$group.size.at.birth)-min(d$group.size.at.birth))
d$n_drought <- (d$biomassRev - min(d$biomassRev))/(max(d$biomassRev)-min(d$biomassRev))
d$n_ibi <-(d$ibidays - min(d$ibidays))/(max(d$ibidays)-min(d$ibidays))
d$n_momlossadj <- (d$age.at.mom.death - min(d$age.at.mom.death))/(max(d$age.at.mom.death)-min(d$age.at.mom.death))
# Following is for separate ela components model:
# d$n_grpsize <- (d$group.size.at.birth - min(d$group.size.at.birth))/(max(d$group.size.at.birth)-min(d$group.size.at.birth))
# d$n_drought <- (d$biomass.year - min(d$biomass.year))/(max(d$biomass.year)-min(d$biomass.year))
# d$n_ibi <-(d$ibidays - min(d$ibidays))/(max(d$ibidays)-min(d$ibidays))
# d$n_momlossadj <- (d$age.at.mom.s.death.sep - min(d$age.at.mom.s.death.sep))/(max(d$age.at.mom.s.death.sep)-min(d$age.at.mom.s.death.sep))
d$s_grpsize <- (d$n_grpsize - mean(d$n_grpsize))/sd(d$n_grpsize)
d$s_drought <- (d$n_drought - mean(d$n_drought))/sd(d$n_drought)
d$s_ibi <- (d$n_ibi  - mean(d$n_ibi))/sd(d$n_ibi)
d$s_momlossadj <- (d$n_momlossadj  - mean(d$n_momlossadj))/sd(d$n_momlossadj)

d$adversity.cumul <- d$n_grpsize + d$n_drought + d$n_ibi + d$n_momlossadj + d$Primiparous
d$s_adversity.cumul <- (d$adversity.cumul - mean(d$adversity.cumul))/sd(d$adversity.cumul)
d$s_adversity.binary <-(d$BinaryIndexAll- mean(d$BinaryIndexAll))/sd(d$BinaryIndexAll)
d$s_ageopuntia <- (d$age.at.opuntia.benefit - mean(d$age.at.opuntia.benefit))/sd(d$age.at.opuntia.benefit)


d$GC.ng.g <- as.numeric(d$GC.ng.g)
d$s_GC <- (d$GC.ng.g - mean(d$GC.ng.g))/sd(d$GC.ng.g) #check distribution
d$logGC <- log10(d$GC.ng.g)
d$s_logGC <- (d$logGC - mean(d$logGC))/sd(d$logGC) #check distribution

d$Mom <- as.factor(d$Mom)
d$Mom <- as.integer(as.factor(d$Mom))

#GCM model 1
GC_model1 <- map2stan(
  alist(
    
    log_gc ~ dnorm(mu , sigma ),
    
    mu <- ap + ap_Focal[Focal] +  bp_female*female + bp_grpsize*grpsize +
      bp_ageS*infageS + bp_age*infage + bp_rank*relrank + bp_momage*momage +
      bp_challenges*challenges + bp_grass*grass + 
      bp_ageopuntia*ageop + bp_adv*advNorm   ,
    
    c(ap_Focal)[Focal] ~ dnorm(0,sigma_focal) ,
    c(ap,bp_female,bp_grpsize,bp_age,bp_ageS,bp_rank,bp_momage,
      bp_challenges,bp_grass,bp_ageopuntia,bp_adv) ~ dnorm(0,1) ,
    c(sigma_focal)  ~ dcauchy(0,1)  ,
    sigma ~ dcauchy(0,2) 
    
  ),
  data=list(
    log_gc=d$s_GC,
    Focal=d$Mom,
    grpsize=d$s_grpsize_cur,
    ageop = d$s_ageopuntia2,
    advNorm = d$s_adversity.cumul ,
    female = d$female,
    challenges = d$s_challenges,
    grass = d$s_grass,
    relrank = d$s_rank,
    momage = d$s_MomAge,
    infageS = d$s_infAgeSquared,
    infage = d$s_infAge
  ),
  
  cores=3 , chains=3 , warmup=1500, iter=5000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

#GCM model 2
GC_model2 <- map2stan(
  alist(
    
    log_gc ~ dnorm(mu , sigma ),
    
    mu <- ap + ap_Focal[Focal] +  bp_female*female + bp_grpsize*grpsize +
      bp_age*infage + bp_ageS*infageS + bp_rank*relrank + bp_momage*momage +
      bp_challenges*challenges + bp_grass*grass + 
      bp_ageopuntia*ageop + bp_adv*advNorm + bp_rank_adv*relrank*advNorm  ,
    
    c(ap_Focal)[Focal] ~ dnorm(0,sigma_focal) ,
    c(ap,bp_age,bp_ageS,bp_female,bp_grpsize,bp_rank,bp_momage,
      bp_challenges,bp_grass,bp_ageopuntia,bp_adv,bp_rank_adv) ~ dnorm(0,1) ,
    c(sigma_focal)  ~ dcauchy(0,1)  ,
    sigma ~ dcauchy(0,2) 
    
  ),
  data=list(
    log_gc=d$s_GC,
    Focal=d$Mom,
    grpsize=d$s_grpsize_cur,
    ageop = d$s_ageopuntia2,
    advNorm = d$s_adversity.cumul ,
    female = d$female,
    challenges = d$s_challenges,
    grass = d$s_grass,
    relrank = d$s_rank,
    momage = d$s_MomAge,
    infageS = d$s_infAgeSquared,
    infage = d$s_infAge
  ),
  
  cores=3 , chains=3 , warmup=1500, iter=5000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

#GCM separate model 1 (rank and ela)
GC1_separate <- map2stan(
  alist(
    
    log_gc ~ dnorm(mu , sigma ),
    
    mu <- ap + ap_Focal[Focal] + bp_age*infage + bp_ageS*infageS + bp_momage*momage + bp_female*female + 
      bp_ageop*ageop + bp_grpsizeCur*grpsizecur + bp_challenge*challenge +
      bp_grass*grass +  bp_rank*relrank +
      bp_grpsizeELA*grpsizeELA + bp_drought*drought + 
      bp_ibiELA*ibiELA + bp_momloss*momloss + bp_firstbirth*first,   
    
    c(ap_Focal)[Focal] ~ dnorm(0,sigma_focal) ,
    c(ap,bp_age,bp_ageS,bp_momage,bp_ageop,bp_grpsizeCur,bp_challenge,bp_grass,
      bp_female,bp_rank,bp_grpsizeELA,bp_drought,
      bp_ibiELA,bp_momloss,bp_firstbirth) ~ dnorm(0,1) ,
    
    c(sigma_focal)  ~ dcauchy(0,1)  ,
    
    sigma ~ dcauchy(0,2) 
    
  ),
  data=list(
    log_gc=d$s_GC,
    Focal=d$Mom,
    infage = d$s_infAge,
    infageS = d$s_infAgeSquared,
    momage = d$s_MomAge,
    grpsizecur=d$s_grpsize_cur,
    challenge=d$s_challenges,
    grass=d$s_grass,
    relrank=d$s_rank,
    ageop = d$s_ageopuntia,
    grpsizeELA = d$s_grpsize,
    drought = d$s_drought,
    ibiELA = d$s_ibi,
    momloss = d$s_momlossadj,
    female = d$female,
    first = d$Primiparous
  ),
  
  cores=3 , chains=3 , warmup=1500, iter=10000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

#GCM separate model 2 (rank x ela)
GC2_separate <- map2stan(
  alist(
    
    log_gc ~ dnorm(mu , sigma ),
    
    mu <- ap + ap_Focal[Focal] + bp_age*infage + bp_ageS*infageS + bp_momage*momage +
      bp_ageop*ageop + bp_grpsizeCur*grpsizecur + bp_challenge*challenge +
      bp_grass*grass + bp_female*female + bp_rank*relrank +
      bp_grpsizeELA*grpsizeELA + bp_drought*drought + 
      bp_ibiELA*ibiELA + bp_momloss*momloss + bp_firstbirth*first +
      (bp_rank_grpsizeELA*grpsizeELA + bp_rank_drought*drought + 
         bp_rank_ibiELA*ibiELA + bp_rank_momloss*momloss + bp_rank_firstbirth*first)*relrank ,   
    
    c(ap_Focal)[Focal] ~ dnorm(0,sigma_focal) ,
    c(ap,bp_age,bp_ageS,bp_momage,bp_ageop,bp_grpsizeCur,bp_challenge,bp_grass,
      bp_female,bp_rank,bp_grpsizebirth, bp_grpsizeELA,bp_drought,
      bp_ibiELA,bp_momloss,bp_firstbirth,bp_rank_grpsizeELA,bp_rank_drought,
      bp_rank_ibiELA,bp_rank_momloss,bp_rank_firstbirth) ~ dnorm(0,1) ,
    
    c(sigma_focal)  ~ dcauchy(0,1)  ,
    
    sigma ~ dcauchy(0,2) 
    
  ),
  data=list(
    log_gc=d$s_GC,
    Focal=d$Mom,
    infage = d$s_infAge,
    infageS = d$s_infAgeSquared,
    momage = d$s_MomAge,
    grpsizecur=d$s_grpsize_cur,
    challenge=d$s_challenges,
    grass=d$s_grass,
    relrank=d$s_rank,
    ageop = d$s_ageopuntia,
    grpsizeELA = d$s_grpsize,
    drought = d$s_drought,
    ibiELA = d$s_ibi,
    momloss = d$s_momlossadj,
    female = d$female,
    first = d$Primiparous
  ),
  
  cores=3 , chains=3 , warmup=1500, iter=10000, WAIC=TRUE, types=list(adapt.delta=0.99)
)



##################
###            ###
###   Plots    ###
###            ###
##################



#GCM index plots
a_foc_z <- matrix(0,1000,length(unique(d$Mom)))
age.seq=seq(min(d$s_adversity.cumul),max(d$s_adversity.cumul),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  infageS=rep(mean(d$s_infAgeSquared),length(age.seq)),
  infage=rep(mean(d$s_infAge),length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenges=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  relrank=rep(mean(d$s_rank),length(age.seq)),
  momage=rep(mean(d$s_MomAge),length(age.seq)),
  advNorm=age.seq
)


LM1 <- link(GC_model1, n=1000 , data=d.pred1, replace=
              list(am_Focal=a_foc_z), WAIC=TRUE)
LM2 <- link(GC_model2, n=1000 , data=d.pred1, replace=
              list(am_Focal=a_foc_z), WAIC=TRUE)

w <- compare(GC_model1,GC_model2 , sort=FALSE)@output$weight  ##extract weights from table
w <- w/sum(w)  ##adds up to 1?
idw <- round( w * 1000 ) ###round weights so extracted samples add to 1000

PM1 <- LM1
PM2 <- LM2

#make a vector according to model weight
pred_1 <- PM1 #first one
pred_1[1:idw[2],] <- PM2[1:idw[2],] #overwrite first one with eighted predicitons from pm2

pred1 <- pred_1[sample(nrow(pred_1)),] #randomly mixes
pred1.median <- apply(pred1, 2, median )
pred1.HPDI <- apply(pred1, 2, HPDI )

par(mfrow = c(1, 1), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))

xxa=age.seq
yya=as.vector(t(pred1[1:1000,]))
smoothScatter(rep(xxa,1000),yya,xlim=c(min(d$s_adversity.cumul),max(d$s_adversity.cumul)), colramp=colorRampPalette(c("white","#33CCFF")),nbin=200,transformation = function(x) x^.5,ylab='',cex=1.2,xaxt='n',yaxt='n',ylim=c(-1.7,3),nrpoints=0)
points(s_GC ~ s_adversity.cumul, data=ddM , col=alpha("#33CCFF",0.6),pch=16, cex=0.9)

lines( age.seq , pred1.median ,lwd=1,col="black")
lines(age.seq, pred1.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1.HPDI[2,],lty=2,lwd=1)

#axis(2,cex.axis=.9)

lab.a=seq(1,2.6,.5)
lab.sa=(lab.a-mean(d$adversity.cumul))/sd(d$adversity.cumul)

lab.b=seq(300,5300,500)
lab.sb=(lab.b-mean(d$GC.ng.g))/sd(d$GC.ng.g)

axis(1,labels=NA, at=lab.sa, tck=-0.01)
mtext(lab.a,at=lab.sa,side=1, line=.5)

axis(2,labels=NA, at=lab.sb, tck=-0.01)
mtext(lab.b,at=lab.sb,side=2, line=.5)

mtext( "Cumulative early life adversity score" , side=1 , line=2, outer=FALSE , cex=1.5)
mtext("GCM ng/g", side=2, line=2,cex=1.5, outer=TRUE)


#GCM interaction with rank
a_foc_z <- matrix(0,1000,length(unique(d$Mom)))
age.seq=seq(min(d$s_adversity.cumul),max(d$s_adversity.cumul),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  infage=rep(mean(d$s_infAge),length(age.seq)),
  infageS=rep(mean(d$s_infAgeSquared),length(age.seq)),
  momage=rep(mean(d$s_MomAge),length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenges=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  relrank=rep(-1,length(age.seq)),
  advNorm=age.seq
)

LM1a <- link(GC_model1, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2a <- link(GC_model2, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)

w <- compare(GC_model1,GC_model2 , sort=FALSE)@output$weight
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
  infage=rep(mean(d$s_infAge),length(age.seq)),
  infageS=rep(mean(d$s_infAgeSquared),length(age.seq)),
  momage=rep(mean(d$s_MomAge),length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenges=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  relrank=rep(1,length(age.seq)),
  advNorm=age.seq
)

LM1b <- link(GC_model1, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2b <- link(GC_model2, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)


w <- w/sum(w)
idw <- round( w * 1000 )

PM1b <- LM1b
PM2b <- LM2b


pred_1b <- PM1b
pred_1b[1:idw[2],] <- PM2b[1:idw[2],]

pred1b <- pred_1b[sample(nrow(pred_1b)),]
pred1b.median <- apply(pred1b, 2, median )
pred1b.HPDI <-apply(pred1b, 2, HPDI )

par(mfrow = c(1, 2), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))

plot( s_GC[d$s_rank==-1] ~ s_adversity.cumul[d$s_rank==-1] , data=ddM , col=alpha("red",0.5), pch=16 ,ylim=c(-1.5,3),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_adversity.cumul),max(d$s_adversity.cumul)))
pred1a.lines=pred1a[sample(nrow(pred1a),100,replace=F),]
for (i in 1:100){
  lines( age.seq , pred1a.lines[i,] ,lwd=2,col=alpha("red",0.1))
}
lines( age.seq , pred1a.median ,lwd=1,col="black")
lines(age.seq, pred1a.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1a.HPDI[2,],lty=2,lwd=1)
lab.a=seq(1,2.6,.5)
lab.sa=(lab.a-mean(d$adversity.cumul))/sd(d$adversity.cumul)

lab.b=seq(300,5300,500)
lab.sb=(lab.b-mean(d$GC.ng.g))/sd(d$GC.ng.g)

axis(1,labels=NA, at=lab.sa, tck=-0.01)
mtext(lab.a,at=lab.sa,side=1, line=.5)

axis(2,labels=NA, at=lab.sb, tck=-0.01)
mtext(lab.b,at=lab.sb,side=2, line=.5)

mtext(lab.a,at=lab.sa,side=1, line=.5)
mtext( "Cumulative early life adversity" , side=1 , line=1.5, outer=FALSE , cex=1.3)
mtext( "Low rank" , side=3 , line=-1, outer=FALSE , cex=1.5)

plot( s_GC[d$s_rank==1] ~ s_adversity.cumul[d$s_rank==1] , data=ddM , col=alpha("#33CCFF",0.5), pch=16 ,ylim=c(-1.5,3),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_adversity.cumul),max(d$s_adversity.cumul)))
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
mtext("GCM ng/g", side=2, line=2,cex=1.5, outer=TRUE)

## example code for the supplementary plots is in nursing and carrying file

