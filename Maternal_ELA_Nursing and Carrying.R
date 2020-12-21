
########################
########################
## Nursing, Carrying ###
########################
########################


d <- read.csv("nursing carrying ELA.csv", header = TRUE, na.strings = "")

### parity
d$multiparous <- ifelse(d$parity=="M" , 1 , 0 )
d$primiparous <- ifelse(d$parity=="P" , 1 , 0 )
d <- subset(d, multiparous==1)
### infant sex
d$female <- ifelse(d$sex=="f" , 1 , 0 )
d$male <- ifelse(d$sex=="m" , 1 , 0 )
### groups
d$phg <- ifelse(d$troop=="phg" , 1 , 0 )
d$enk <- ifelse(d$troop=="enk" , 1 , 0 )
d$ynt <- ifelse(d$troop=="ynt" , 1 , 0 )
d$nmu <- ifelse(d$troop=="nmu" , 1 , 0 )
#age
d$s_momage <- (d$FocalAge- mean(d$FocalAge))/sd(d$FocalAge)
d$s_age <- (d$age - mean(d$age))/sd(d$age) #infant age
d <- subset(d, age < 365)

#current conditions:
d <- d[complete.cases(d$avgTotal), ]
d$s_grass <- (d$avgTotal - mean(d$avgTotal))/sd(d$avgTotal)
d$Human_count <- as.numeric(d$Human_count)
d$s_human <- (d$Human_count - mean(d$Human_count))/sd(d$Human_count)
d$StrangeMales_count <- as.numeric(d$StrangeMales_count)
d$s_strangemales <- (d$StrangeMales_count - mean(d$StrangeMales_count))/sd(d$StrangeMales_count)
d$challenges <- d$Human_count + d$StrangeMales_count
d$s_challenges <- (d$challenges - mean(d$challenges))/sd(d$challenges)
d$s_grpsize_cur <- (d$group_size - mean(d$group_size))/sd(d$group_size)
d$s_ageopuntia <- (d$age.at.opuntia.benefit - mean(d$age.at.opuntia.benefit))/sd(d$age.at.opuntia.benefit)

#early life adversity
d <- d[complete.cases(d$group.size.at.birth,d$biomassRev,d$ibidays,d$ageatmomdeathRev,d$Primiparous,d$relRank), ]
d$n_grpsize <- (d$group.size.at.birth - min(d$group.size.at.birth))/(max(d$group.size.at.birth)-min(d$group.size.at.birth))
d$n_drought <- (d$biomassRev - min(d$biomassRev))/(max(d$biomassRev)-min(d$biomassRev))
d$n_ibi <-(d$ibidays - min(d$ibidays))/(max(d$ibidays)-min(d$ibidays))
d$n_momlossadj <- (d$age.at.mom.death - min(d$age.at.mom.death))/(max(d$age.at.mom.death)-min(d$age.at.mom.death))
##to model componenets of early life adversity separately:
# d$n_grpsize <- (d$group.size.at.birth - min(d$group.size.at.birth))/(max(d$group.size.at.birth)-min(d$group.size.at.birth))
# d$n_drought <- (d$biomass.year - min(d$biomass.year))/(max(d$biomass.year)-min(d$biomass.year))
# d$n_ibi <-(d$ibidays - min(d$ibidays))/(max(d$ibidays)-min(d$ibidays))
# d$n_momlossadj <- (d$age.at.mom.s.death.sep - min(d$age.at.mom.s.death.sep))/(max(d$age.at.mom.s.death.sep)-min(d$age.at.mom.s.death.sep))
d$s_grpsize <- (d$n_grpsize - mean(d$n_grpsize))/sd(d$n_grpsize)
d$s_drought <- (d$n_drought - mean(d$n_drought))/sd(d$n_drought)
d$s_ibi <- (d$n_ibi - mean(d$n_ibi))/sd(d$n_ibi)
d$s_momlossadj <- (d$n_momlossadj - mean(d$n_momlossadj))/sd(d$n_momlossadj)

#continuous and binary indices
d$adversity.cumul <- d$n_grpsize + d$n_drought + d$n_ibi + d$n_momlossadj + d$Primiparous
d$s_adversity.cumul <- (d$adversity.cumul - mean(d$adversity.cumul))/sd(d$adversity.cumul)
d$s_adversity.binary <- (d$BinaryIndexAll - mean(d$BinaryIndexAll))/sd(d$BinaryIndexAll)

d$Focal <- as.factor(d$Focal)
d$Focal <- as.integer(as.factor(d$Focal))

#Nursing model 1 (rank and ELA)
nursing1 <- map2stan(
  alist(
    
    bndur ~ dzagamma2(p, mu , scale ),
    
    logit(p) <- ap + ap_Focal[Focal] + bp_age*age  + bp_grpsize*grpsize + bp_momage*momage +
      bp_opuntia*ageop + bp_female*female + bp_challenge*challenge + bp_grass*grass + 
      bp_adv*advNorm + bp_rank*relrank  ,
    
    log(mu) <- am + am_Focal[Focal] + bm_age*age  + 
      bm_adv*advNorm + bm_rank*relrank +
      bm_challenge*challenge + bm_grass*grass + bm_grpsize*grpsize + bm_momage*momage +
      bm_opuntia*ageop + bm_female*female ,
    
    c(ap_Focal,am_Focal)[Focal] ~ dmvnormNC(sigma_focal,Rho_focal) ,
    c(ap,am,bp_age,bp_adv,bp_rank,bp_challenge,bp_grass,bp_grpsize,bp_momage,
      bp_opuntia,bp_female,
      
      bm_age,bm_adv,bm_rank,bm_challenge,bm_grass,bm_grpsize,bm_momage,
      bm_opuntia,bm_female) ~ dnorm(0,2) ,
    sigma_focal  ~ dcauchy(0,2)  ,
    Rho_focal ~ dlkjcorr(3) ,
    scale ~ dcauchy(0,2)  
    
  ),
  data=list(
    bndur=d$Bnminsrate,
    Focal=d$Focal,
    age = d$s_age,  
    female = d$female,
    challenge = d$s_challenges,
    grass = d$s_grass,
    grpsize = d$s_grpsize_cur,
    momage = d$s_momage,
    relrank = d$s_rank,
    ageop = d$s_ageopuntia,
    advNorm = d$s_adversity.cumul
  ),
  
  cores=3 , chains=2 , warmup=1500, iter=3000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

#Nursing model 2 (rank X ela)
nursing2 <- map2stan(
  alist(
    
    bndur ~ dzagamma2(p, mu , scale ),
    
    logit(p) <- ap + ap_Focal[Focal] + bp_age*age  + bp_grpsize*grpsize + bp_momage*momage + 
      bp_opuntia*ageop + bp_female*female + bp_challenge*challenge + bp_grass*grass +  
      bp_adv*advNorm + bp_rank*relrank +  
      bp_rank_adv*relrank*advNorm  ,
    
    log(mu) <- am + am_Focal[Focal] + bm_age*age  + 
      bm_adv*advNorm + bm_rank*relrank +
      bm_challenge*challenge + bm_grass*grass + bm_grpsize*grpsize + bm_momage*momage +
      bm_rank_adv*relrank*advNorm + 
      bm_opuntia*ageop + bm_female*female ,
    
    c(ap_Focal,am_Focal)[Focal] ~ dmvnormNC(sigma_focal,Rho_focal) ,
    c(ap,am,bp_age,bp_adv,bp_rank,bp_challenge,bp_grass,bp_grpsize,bp_momage,
      bp_rank_adv,
      bp_opuntia,bp_female,
      
      bm_age,bm_adv,bm_rank,bm_challenge,bm_grass,bm_grpsize,bm_momage,
      bm_rank_adv,
      bm_opuntia,bm_female) ~ dnorm(0,2) ,
    sigma_focal  ~ dcauchy(0,2)  ,
    Rho_focal ~ dlkjcorr(3) ,
    scale ~ dcauchy(0,2)  
    
  ),
  data=list(
    bndur=d$Bnminsrate,
    Focal=d$Focal,
    age = d$s_age,  
    female = d$female,
    challenge = d$s_challenges,
    grass = d$s_grass,
    grpsize = d$s_grpsize_cur,
    momage = d$s_momage,
    relrank = d$s_rank,
    ageop = d$s_ageopuntia,
    advNorm = d$s_adversity.cumul
  ),
  
  cores=3 , chains=2 , warmup=1500, iter=3000, WAIC=TRUE, types=list(adapt.delta=0.99)
)


###nursing models with separate ela components
nursing1_separate <- map2stan(
  alist(
    
    bndur ~ dzagamma2(p, mu , scale ),
    
    logit(p) <- ap + ap_Focal[Focal] + bp_age*age  + bp_momage*momage + bp_grpsize*grpsize +
      bp_opuntia*ageop + bp_female*female + bp_challenge*challenge + bp_grass*grass +  
      bp_rank*relrank + bp_grpsizeELA*grpsizeELA + bp_drought*drought + 
      bp_ibiELA*ibiELA + bp_momloss*momloss + bp_firstbirth*first ,
    
    log(mu) <- am + am_Focal[Focal] + bm_age*age  + bm_momage*momage + bm_grpsize*grpsize +
      bm_opuntia*ageop + bm_female*female + bm_challenge*challenge + bm_grass*grass + 
      bm_rank*relrank + bm_grpsizeELA*grpsizeELA + bm_drought*drought + 
      bm_ibiELA*ibiELA + bm_momloss*momloss + bm_firstbirth*first  ,
    
    c(ap_Focal,am_Focal)[Focal] ~ dmvnormNC(sigma_focal,Rho_focal) ,
    c(ap,am,bp_age,bp_momage,bp_rank,bp_challenge,bp_grass,bp_grpsize,
      bp_opuntia,bp_female,
      bp_grpsizeELA,bp_drought,bp_ibiELA,bp_momloss,bp_firstbirth,
      
      bm_age,bm_momage,bm_rank,bm_challenge,bm_grass,bm_grpsize,
      bm_opuntia,bm_female,
      bm_grpsizeELA,bm_drought,bm_ibiELA,bm_momloss,bm_firstbirth) ~ dnorm(0,2) ,
    sigma_focal  ~ dcauchy(0,2)  ,
    Rho_focal ~ dlkjcorr(3) ,
    scale ~ dcauchy(0,2)  
    
  ),
  data=list(
    bndur=d$Bnminsrate,
    Focal=d$Focal,
    age = d$s_age,  
    momage = d$s_momage,
    female = d$female,
    challenge = d$s_challenges,
    grass = d$s_grass,
    grpsize = d$s_grpsize_cur,
    relrank = d$s_rank,
    ageop = d$s_ageopuntia,
    grpsizeELA = d$s_grpsize,
    drought = d$s_drought,
    ibiELA = d$s_ibi,
    momloss = d$s_momlossadj,
    first = d$Primiparous
  ),
  
  cores=3 , chains=2 , warmup=1500, iter=3000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

#separate 2 (rank X ela)
nursing2_separate <- map2stan(
  alist(
    
    bndur ~ dzagamma2(p, mu , scale ),
    
    logit(p) <- ap + ap_Focal[Focal] + bp_age*age  + bp_momage*momage + bp_grpsize*grpsize +
      bp_opuntia*ageop + bp_female*female + bp_challenge*challenge + bp_grass*grass +  
      bp_rank*relrank + bp_grpsizeELA*grpsizeELA + bp_drought*drought + 
      bp_ibiELA*ibiELA + bp_momloss*momloss + bp_firstbirth*first +
      (bp_rank_grpsizeELA*grpsizeELA + bp_rank_drought*drought + 
         bp_rank_ibiELA*ibiELA + bp_rank_momloss*momloss + bp_rank_firstbirth*first)*relrank ,
    
    log(mu) <- am + am_Focal[Focal] + bm_age*age  + bm_momage*momage + bm_grpsize*grpsize +
      bm_opuntia*ageop + bm_female*female + bm_challenge*challenge + bm_grass*grass + 
      bm_rank*relrank + bm_grpsizeELA*grpsizeELA + bm_drought*drought + 
      bm_ibiELA*ibiELA + bm_momloss*momloss + bm_firstbirth*first +
      (bm_rank_grpsizeELA*grpsizeELA + bm_rank_drought*drought + 
         bm_rank_ibiELA*ibiELA + bm_rank_momloss*momloss + bm_rank_firstbirth*first)*relrank  ,
    
    c(ap_Focal,am_Focal)[Focal] ~ dmvnormNC(sigma_focal,Rho_focal) ,
    c(ap,am,bp_age,bp_momage,bp_rank,bp_challenge,bp_grass,bp_grpsize,
      bp_opuntia,bp_female,
      bp_grpsizeELA,bp_drought,bp_ibiELA,bp_momloss,bp_firstbirth,
      bp_rank_grpsizeELA,bp_rank_drought,bp_rank_ibiELA,bp_rank_momloss,bp_rank_firstbirth,
        
      
      bm_age,bm_momage,bm_rank,bm_challenge,bm_grass,bm_grpsize,
      bm_opuntia,bm_female,
      bm_grpsizeELA,bm_drought,bm_ibiELA,bm_momloss,bm_firstbirth,
      bm_rank_grpsizeELA,bm_rank_drought,bm_rank_ibiELA,bm_rank_momloss,bm_rank_firstbirth) ~ dnorm(0,2) ,
    sigma_focal  ~ dcauchy(0,2)  ,
    Rho_focal ~ dlkjcorr(3) ,
    scale ~ dcauchy(0,2)  
    
  ),
  data=list(
    bndur=d$Bnminsrate,
    Focal=d$Focal,
    age = d$s_age,  
    momage = d$s_momage,
    female = d$female,
    challenge = d$s_challenges,
    grass = d$s_grass,
    grpsize = d$s_grpsize_cur,
    relrank = d$s_rank,
    ageop = d$s_ageopuntia,
    grpsizeELA = d$s_grpsize,
    drought = d$s_drought,
    ibiELA = d$s_ibi,
    momloss = d$s_momlossadj,
    first = d$Primiparous
  ),
  
  cores=3 , chains=2 , warmup=1500, iter=3000, WAIC=TRUE, types=list(adapt.delta=0.99)
)


#carrying model 1 (rank and ela)

carrying1 <- map2stan(
  alist(
    
    bhdur ~ dzagamma2(p, mu , scale ),
    
    logit(p) <- ap + ap_Focal[Focal] + bp_age*age  + 
      bp_adv*advNorm + bp_rank*relrank +  
      bp_challenge*challenge + bp_grass*grass + bp_grpsize*grpsize + bp_momage*momage +
      bp_opuntia*ageop + bp_female*female ,
    
    log(mu) <- am + am_Focal[Focal] + bm_age*age  + 
      bm_adv*advNorm + bm_rank*relrank +
      bm_challenge*challenge + bm_grass*grass + bm_grpsize*grpsize + bm_momage*momage +
      bm_opuntia*ageop + bm_female*female ,
    
    c(ap_Focal,am_Focal)[Focal] ~ dmvnormNC(sigma_focal,Rho_focal) ,
    c(ap,am,bp_age,bp_adv,bp_rank,bp_challenge,bp_grass,bp_grpsize,bp_momage,
      bp_opuntia,bp_female,
      
      bm_age,bm_adv,bm_rank,bm_challenge,bm_grass,bm_grpsize,bm_momage,
      bm_opuntia,bm_female) ~ dnorm(0,2) ,
    sigma_focal  ~ dcauchy(0,2)  ,
    Rho_focal ~ dlkjcorr(3) ,
    scale ~ dcauchy(0,2)  
    
  ),
  data=list(
    bhdur=d$Bhminsrate,
    Focal=d$Focal,
    age = d$s_age,  
    female = d$female,
    challenge = d$s_challenges,
    grass = d$s_grass,
    grpsize = d$s_grpsize_cur,
    momage = d$s_momage,
    relrank = d$s_rank,
    ageop = d$s_ageopuntia,
    advNorm = d$s_adversity.cumul
  ),
  
  cores=3 , chains=2 , warmup=1500, iter=3000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

#carrying model 2 (rank X ela)

carrying2 <- map2stan(
  alist(
    
    bhdur ~ dzagamma2(p, mu , scale ),
    
    logit(p) <- ap + ap_Focal[Focal] + bp_age*age  + 
      bp_adv*advNorm + bp_rank*relrank +  
      bp_challenge*challenge + bp_grass*grass + bp_grpsize*grpsize + bp_momage*momage +  
      bp_rank_adv*relrank*advNorm + 
      bp_opuntia*ageop + bp_female*female ,
    
    log(mu) <- am + am_Focal[Focal] + bm_age*age  + 
      bm_adv*advNorm + bm_rank*relrank +
      bm_challenge*challenge + bm_grass*grass + bm_grpsize*grpsize + bm_momage*momage + 
      bm_rank_adv*relrank*advNorm + 
      bm_opuntia*ageop + bm_female*female ,
    
    c(ap_Focal,am_Focal)[Focal] ~ dmvnormNC(sigma_focal,Rho_focal) ,
    c(ap,am,bp_age,bp_adv,bp_rank,bp_challenge,bp_grass,bp_grpsize,bp_momage,
      bp_rank_adv,
      bp_opuntia,bp_female,
      
      bm_age,bm_adv,bm_rank,bm_challenge,bm_grass,bm_grpsize,bm_momage,
      bm_rank_adv,
      bm_opuntia,bm_female) ~ dnorm(0,2) ,
    sigma_focal  ~ dcauchy(0,2)  ,
    Rho_focal ~ dlkjcorr(3) ,
    scale ~ dcauchy(0,2)  
    
  ),
  data=list(
    bhdur=d$Bhminsrate,
    Focal=d$Focal,
    age = d$s_age,  
    female = d$female,
    challenge = d$s_challenges,
    grass = d$s_grass,
    grpsize = d$s_grpsize_cur,
    momage = d$s_momage,
    relrank = d$s_rank,
    ageop = d$s_ageopuntia,
    advNorm = d$s_adversity.cumul
  ),
  
  cores=3 , chains=2 , warmup=1500, iter=3000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

###########################################################################################
## Note: Binary ELA models are the same, but the advNorm variable = d$s_adversity.binary ##
###########################################################################################

# carrying models with separate ela components:
 
carrying1_separate <- map2stan(
  alist(
    
    bhdur ~ dzagamma2(p, mu , scale ),
    
    logit(p) <- ap + ap_Focal[Focal] + bp_age*age  + bp_momage*momage + bp_grpsize*grpsize +
      bp_opuntia*ageop + bp_female*female + bp_challenge*challenge + bp_grass*grass +  
      bp_rank*relrank + bp_grpsizeELA*grpsizeELA + bp_drought*drought + 
      bp_ibiELA*ibiELA + bp_momloss*momloss + bp_firstbirth*first ,
    
    log(mu) <- am + am_Focal[Focal] + bm_age*age  + bm_momage*momage + bm_grpsize*grpsize +
      bm_opuntia*ageop + bm_female*female + bm_challenge*challenge + bm_grass*grass + 
      bm_rank*relrank +bm_grpsizeELA*grpsizeELA + bm_drought*drought + 
      bm_ibiELA*ibiELA + bm_momloss*momloss + bm_firstbirth*first,
    
    c(ap_Focal,am_Focal)[Focal] ~ dmvnormNC(sigma_focal,Rho_focal) ,
    c(ap,am,bp_age,bp_momage,bp_rank,bp_challenge,bp_grass,bp_grpsize,
      bp_opuntia,bp_female,
      bp_grpsizeELA,bp_drought,bp_ibiELA,bp_momloss,bp_firstbirth,
      
      
      bm_age,bm_momage,bm_rank,bm_challenge,bm_grass,bm_grpsize,
      bm_opuntia,bm_female,
      bm_grpsizeELA,bm_drought,bm_ibiELA,bm_momloss,bm_firstbirth) ~ dnorm(0,2) ,
    sigma_focal  ~ dcauchy(0,2)  ,
    Rho_focal ~ dlkjcorr(3) ,
    scale ~ dcauchy(0,2)  
    
  ),
  data=list(
    bhdur=d$Bhminsrate,
    Focal=d$Focal,
    age = d$s_age,  
    momage = d$s_momage,
    female = d$female,
    challenge = d$s_challenges,
    grass = d$s_grass,
    grpsize = d$s_grpsize_cur,
    relrank = d$s_rank,
    ageop = d$s_ageopuntia,
    grpsizeELA = d$s_grpsize,
    drought = d$s_drought,
    ibiELA = d$s_ibi,
    momloss = d$s_momlossadj,
    first = d$Primiparous
  ),
  
  cores=3 , chains=2 , warmup=1500, iter=3000, WAIC=TRUE, types=list(adapt.delta=0.99)
)

# separate 2 (rank X ela)
carrying2_separate <- map2stan(
  alist(
    
    bhdur ~ dzagamma2(p, mu , scale ),
    
    logit(p) <- ap + ap_Focal[Focal] + bp_age*age  + bp_momage*momage + bp_grpsize*grpsize +
      bp_opuntia*ageop + bp_female*female + bp_challenge*challenge + bp_grass*grass +  
      bp_rank*relrank + bp_grpsizeELA*grpsizeELA + bp_drought*drought + 
      bp_ibiELA*ibiELA + bp_momloss*momloss + bp_firstbirth*first +
      (bp_rank_grpsizeELA*grpsizeELA + bp_rank_drought*drought + 
         bp_rank_ibiELA*ibiELA + bp_rank_momloss*momloss + bp_rank_firstbirth*first)*relrank ,
    
    log(mu) <- am + am_Focal[Focal] + bm_age*age  + bm_momage*momage + bm_grpsize*grpsize +
      bm_opuntia*ageop + bm_female*female + bm_challenge*challenge + bm_grass*grass + 
      bm_rank*relrank + bm_grpsizeELA*grpsizeELA + bm_drought*drought + 
      bm_ibiELA*ibiELA + bm_momloss*momloss + bm_firstbirth*first +
      (bm_rank_grpsizeELA*grpsizeELA + bm_rank_drought*drought + 
         bm_rank_ibiELA*ibiELA + bm_rank_momloss*momloss + bm_rank_firstbirth*first)*relrank  ,
    
    c(ap_Focal,am_Focal)[Focal] ~ dmvnormNC(sigma_focal,Rho_focal) ,
    c(ap,am,bp_age,bp_momage,bp_rank,bp_challenge,bp_grass,bp_grpsize,
      bp_opuntia,bp_female,
      bp_grpsizeELA,bp_drought,bp_ibiELA,bp_momloss,bp_firstbirth,
      bp_rank_grpsizeELA,bp_rank_drought,bp_rank_ibiELA,bp_rank_momloss,bp_rank_firstbirth,
      
      bm_age,bm_momage,bm_rank,bm_challenge,bm_grass,bm_grpsize,
      bm_opuntia,bm_female,
      bm_grpsizeELA,bm_drought,bm_ibiELA,bm_momloss,bm_firstbirth,
      bm_rank_grpsizeELA,bm_rank_drought,bm_rank_ibiELA,bm_rank_momloss,bm_rank_firstbirth) ~ dnorm(0,2) ,
    sigma_focal  ~ dcauchy(0,2)  ,
    Rho_focal ~ dlkjcorr(3) ,
    scale ~ dcauchy(0,2)  
    
  ),
  data=list(
    bhdur=d$Bhminsrate,
    Focal=d$Focal,
    age = d$s_age,  
    momage = d$s_momage,
    female = d$female,
    challenge = d$s_challenges,
    grass = d$s_grass,
    grpsize = d$s_grpsize_cur,
    relrank = d$s_rank,
    ageop = d$s_ageopuntia,
    grpsizeELA = d$s_grpsize,
    drought = d$s_drought,
    ibiELA = d$s_ibi,
    momloss = d$s_momlossadj,
    first = d$Primiparous
  ),
  
  cores=3 , chains=2 , warmup=1500, iter=3000, WAIC=TRUE, types=list(adapt.delta=0.99)
)




##################
###            ###
###   Plots    ###
###            ###
##################


# Nursing

a_foc_z <- matrix(0,1000,length(unique(d$Focal)))
age.seq=seq(min(d$s_adversity.cumul),max(d$s_adversity.cumul),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  relrank=rep(mean(d$s_rank),length(age.seq)),
  advNorm=age.seq
)
LM1 <- link(nursing1, n=1000 , data=d.pred1, replace=
              list(am_Focal=a_foc_z), WAIC=TRUE)
LM2 <- link(nursing2, n=1000 , data=d.pred1, replace=
              list(am_Focal=a_foc_z), WAIC=TRUE)

w <- compare(nursing1,nursing2 , sort=FALSE)@output$weight  ##extract weights from table
w <- w/sum(w)  ##make sure everything adds up to 1
idw <- round( w * 1000 ) ###round weights to nearest integer so samples we extract sum up to 1000

#generate predictions from each model that will be averaged
PM1 <- (1-LM1$p)*LM1$mu
PM2 <- (1-LM2$p)*LM2$mu

#make a vector or predictions sampled according to model weight
pred_1 <- PM1
pred_1[1:idw[2],] <- PM2[1:idw[2],]

pred1 <- pred_1[sample(nrow(pred_1)),] #randomly mixes vector
pred1.median <- apply(pred1, 2, median )
pred1.HPDI <-apply(pred1, 2, HPDI )

par(mfrow = c(1, 1), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))
xxa=age.seq
yya=as.vector(t(pred1[1:1000,]))

smoothScatter(rep(xxa,1000),yya,xlim=c(min(d$s_adversity.cumul),max(d$s_adversity.cumul)), colramp=colorRampPalette(c("white","#33CCFF")),nbin=200,transformation = function(x) x^.5,ylab='',cex=1.2,xaxt='n',ylim=c(0,.5),nrpoints=0)
points(Bnminsrate ~ s_adversity.cumul, data=dd , col=alpha("#33CCFF",0.6),pch=16, cex=0.9)
lines( age.seq , pred1.median ,lwd=1,col="black")
lines(age.seq, pred1.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1.HPDI[2,],lty=2,lwd=1)

lab.a=seq(1,2.5,.5)
lab.sa=(lab.a-mean(d$adversity.cumul))/sd(d$adversity.cumul)
axis(1,labels=NA, at=lab.sa, tck=-0.01)
mtext(lab.a,at=lab.sa,side=1, line=.5)

mtext( "Cumulative early life adversity score" , side=1 , line=2, outer=FALSE , cex=1.5)
mtext("Proportion of time observed nursing", side=2, line=2,cex=1.5, outer=TRUE)


### nursing interaction with rank

a_foc_z <- matrix(0,1000,length(unique(d$Focal)))
age.seq=seq(min(d$s_adversity.cumul),max(d$s_adversity.cumul),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  relrank=rep(-1.5,length(age.seq)),
  advNorm=age.seq
)

LM1a <- link(nursing1, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2a <- link(nursing2, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)

w <- compare(nursing1,nursing2 , sort=FALSE)@output$weight
w <- w/sum(w)
idw <- round( w * 1000 )

PM1a <- (1-LM1a$p)*LM1a$mu
PM2a <- (1-LM2a$p)*LM2a$mu


pred_1a <- PM1a
pred_1a[1:idw[2],] <- PM2a[1:idw[2],]

pred1a <- pred_1a[sample(nrow(pred_1a)),]
pred1a.median <- apply(pred1a, 2, median )
pred1a.HPDI <-apply(pred1a, 2, HPDI )

d.pred2 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  mom=rep(mean(d$Focal.),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  relrank=rep(1.5,length(age.seq)),
  advNorm=age.seq
)

LM1b <- link(nursing1, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2b <- link(nursing2, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)


PM1b <- (1-LM1b$p)*LM1b$mu
PM2b <- (1-LM2b$p)*LM2b$mu


pred_1b <- PM1b
pred_1b[1:idw[2],] <- PM2b[1:idw[2],]

pred1b <- pred_1b[sample(nrow(pred_1b)),]
pred1b.median <- apply(pred1b, 2, median )
pred1b.HPDI <-apply(pred1b, 2, HPDI )

par(mfrow = c(1, 2), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))

plot( Bnminsrate[d$s_rank==-1.5] ~ s_adversity.cumul[d$s_rank==-1.5] , data=dd , col=alpha("red",0.5), pch=16 ,ylim=c(0,.5),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_adversity.cumul),max(d$s_adversity.cumul)))
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

plot( Bnminsrate[d$s_rank==1.5] ~ s_adversity.cumul[d$s_rank==1.5] , data=dd , col=alpha("#33CCFF",0.5), pch=16 ,ylim=c(0,.5),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_adversity.cumul),max(d$s_adversity.cumul)))
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
mtext("Proportion of time observed nursing", side=2, line=2,cex=1.5, outer=TRUE)



# carrying
a_foc_z <- matrix(0,1000,length(unique(d$Focal)))
age.seq=seq(min(d$s_adversity.cumul),max(d$s_adversity.cumul),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  relrank=rep(mean(d$s_rank),length(age.seq)),
  advNorm=age.seq
)
LM1 <- link(carrying1, n=1000 , data=d.pred1, replace=
              list(am_Focal=a_foc_z), WAIC=TRUE)
LM2 <- link(carrying2, n=1000 , data=d.pred1, replace=
              list(am_Focal=a_foc_z), WAIC=TRUE)

w <- compare(carrying1,carrying2 , sort=FALSE)@output$weight
w <- w/sum(w)
idw <- round( w * 1000 )

PM1 <- (1-LM1$p)*LM1$mu
PM2 <- (1-LM2$p)*LM2$mu


pred_1 <- PM1
pred_1[1:idw[2],] <- PM2[1:idw[2],]

pred1 <- pred_1[sample(nrow(pred_1)),]
pred1.median <- apply(pred1, 2, median )
pred1.HPDI <-apply(pred1, 2, HPDI )

par(mfrow = c(1, 1), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))
xxa=age.seq
yya=as.vector(t(pred1[1:1000,]))
smoothScatter(rep(xxa,1000),yya,xlim=c(min(d$s_adversity.cumul),max(d$s_adversity.cumul)), colramp=colorRampPalette(c("white","#33CCFF")),nbin=200,transformation = function(x) x^.5,ylab='',cex=1.2,xaxt='n',ylim=c(0,.5),nrpoints=0)
points(Bhminsrate ~ s_adversity.cumul, data=dd , col=alpha("#33CCFF",0.6),pch=16, cex=0.9)
lines( age.seq , pred1.median ,lwd=1,col="black")
lines(age.seq, pred1.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1.HPDI[2,],lty=2,lwd=1)

lab.a=seq(1,2.5,.5)
lab.sa=(lab.a-mean(d$adversity.cumul))/sd(d$adversity.cumul)
axis(1,labels=NA, at=lab.sa, tck=-0.01)
mtext(lab.a,at=lab.sa,side=1, line=.5)

mtext( "Cumulative early life adversity score" , side=1 , line=2, outer=FALSE , cex=1.5)
mtext("Proportion of time observed carrying", side=2, line=2,cex=1.5, outer=TRUE)


#carrying interaction with rank
a_foc_z <- matrix(0,1000,length(unique(d$Focal)))
age.seq=seq(min(d$s_adversity.cumul),max(d$s_adversity.cumul),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  relrank=rep(-1.5,length(age.seq)),
  advNorm=age.seq
)

LM1a <- link(carrying1, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2a <- link(carrying2, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)

w <- compare(carrying1,carrying2 , sort=FALSE)@output$weight
w <- w/sum(w)
idw <- round( w * 1000 )

PM1a <- (1-LM1a$p)*LM1a$mu
PM2a <- (1-LM2a$p)*LM2a$mu


pred_1a <- PM1a
pred_1a[1:idw[2],] <- PM2a[1:idw[2],]

pred1a <- pred_1a[sample(nrow(pred_1a)),]
pred1a.median <- apply(pred1a, 2, median )
pred1a.HPDI <-apply(pred1a, 2, HPDI )

d.pred2 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  relrank=rep(1.5,length(age.seq)),
  advNorm=age.seq
)

LM1b <- link(carrying1, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2b <- link(carrying2, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)


PM1b <- (1-LM1b$p)*LM1b$mu
PM2b <- (1-LM2b$p)*LM2b$mu


pred_1b <- PM1b
pred_1b[1:idw[2],] <- PM2b[1:idw[2],]

pred1b <- pred_1b[sample(nrow(pred_1b)),]
pred1b.median <- apply(pred1b, 2, median )
pred1b.HPDI <-apply(pred1b, 2, HPDI )

par(mfrow = c(1, 2), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))

plot( Bhminsrate[d$s_rank==-1.5] ~ s_adversity.cumul[d$s_rank==-1.5] , data=dd , col=alpha("red",0.5), pch=16 ,ylim=c(0,.4),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_adversity.cumul),max(d$s_adversity.cumul)))
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

plot( Bhminsrate[d$s_rank==1.5] ~ s_adversity.cumul[d$s_rank==1.5] , data=dd , col=alpha("#33CCFF",0.5), pch=16 ,ylim=c(0,.4),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_adversity.cumul),max(d$s_adversity.cumul)))
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
mtext("Proportion of time observed carrying", side=2, line=2,cex=1.5, outer=TRUE)



#####################################################################
#####################################################################
### example of separate ela components by rank in supplementary:  ###
#####################################################################
#####################################################################

## carrying
a_foc_z <- matrix(0,1000,length(unique(d$Focal)))
age.seq=seq(min(d$s_grpsize),max(d$s_grpsize),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  grpsizeELA=age.seq,
  drought=rep(mean(d$s_drought),length(age.seq)),
  ibiELA=rep(mean(d$s_ibi),length(age.seq)),
  momloss=rep(mean(d$s_momlossadj),length(age.seq)),
  first=rep(mean(d$Primiparous),length(age.seq)),
  relrank=rep(-1.5,length(age.seq))
)

LM1a <- link(carrying1_separate, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2a <- link(carrying2_separate, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)

w <- compare(carrying1_separate,carrying2_separate , sort=FALSE)@output$weight
w <- w/sum(w)
idw <- round( w * 1000 )

PM1a <- (1-LM1a$p)*LM1a$mu
PM2a <- (1-LM2a$p)*LM2a$mu


pred_1a <- PM1a
pred_1a[1:idw[2],] <- PM2a[1:idw[2],]

pred1a <- pred_1a[sample(nrow(pred_1a)),]
pred1a.median <- apply(pred1a, 2, median )
pred1a.HPDI <-apply(pred1a, 2, HPDI )

d.pred2 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  grpsizeELA=age.seq,
  drought=rep(mean(d$s_drought),length(age.seq)),
  ibiELA=rep(mean(d$s_ibi),length(age.seq)),
  momloss=rep(mean(d$s_momlossadj),length(age.seq)),
  first=rep(mean(d$Primiparous),length(age.seq)),
  relrank=rep(1.5,length(age.seq))
)

LM1b <- link(carrying1_separate, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2b <- link(carrying2_separate, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)


PM1b <- (1-LM1b$p)*LM1b$mu
PM2b <- (1-LM2b$p)*LM2b$mu


pred_1b <- PM1b
pred_1b[1:idw[2],] <- PM2b[1:idw[2],]

pred1b <- pred_1b[sample(nrow(pred_1b)),]
pred1b.median <- apply(pred1b, 2, median )
pred1b.HPDI <-apply(pred1b, 2, HPDI )

par(mfrow = c(5, 2), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))

plot( Bhminsrate[d$s_rank==-1.5] ~ s_grpsize[d$s_rank==-1.5] , data=dd , col=alpha("red",0.5), pch=16 ,ylim=c(0,.4),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_grpsize),max(d$s_grpsize)))
pred1a.lines=pred1a[sample(nrow(pred1a),100,replace=F),]
for (i in 1:100){
  lines( age.seq , pred1a.lines[i,] ,lwd=2,col=alpha("red",0.1))
}
lines( age.seq , pred1a.median ,lwd=1,col="black")
lines(age.seq, pred1a.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1a.HPDI[2,],lty=2,lwd=1)
axis(2,cex.axis=.9)
mtext( "group size" , side=1 , line=.3, outer=FALSE , cex=1)
mtext( "low rank" , side=3 , line=1, outer=FALSE , cex=1.5)

plot( Bhminsrate[d$s_rank==1.5] ~ s_grpsize[d$s_rank==1.5] , data=dd , col=alpha("#33CCFF",0.5), pch=16 ,ylim=c(0,.4),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_grpsize),max(d$s_grpsize)))
pred1b.lines=pred1b[sample(nrow(pred1b),100,replace=F),]
for (i in 1:100){
  lines( age.seq , pred1b.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.1))
}
lines( age.seq , pred1b.median ,lwd=1,col="black")
lines(age.seq, pred1b.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1b.HPDI[2,],lty=2,lwd=1)

mtext( "group size" , side=1 , line=.3, outer=FALSE , cex=1)
mtext( "high rank" , side=3 , line=1, outer=FALSE , cex=1.5)
mtext("proportion of time observed carrying", side=2, line=2,cex=1.5, outer=TRUE)



a_foc_z <- matrix(0,1000,length(unique(d$Focal)))
age.seq=seq(min(d$s_drought),max(d$s_drought),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  grpsizeELA=rep(mean(d$s_grpsize),length(age.seq)),
  drought=age.seq,
  ibiELA=rep(mean(d$s_ibi),length(age.seq)),
  momloss=rep(mean(d$s_momlossadj),length(age.seq)),
  first=rep(mean(d$Primiparous),length(age.seq)),
  relrank=rep(-1.5,length(age.seq))
)

LM1a <- link(carrying1_separate, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2a <- link(carrying2_separate, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)


PM1a <- (1-LM1a$p)*LM1a$mu
PM2a <- (1-LM2a$p)*LM2a$mu


pred_1a <- PM1a
pred_1a[1:idw[2],] <- PM2a[1:idw[2],]

pred1a <- pred_1a[sample(nrow(pred_1a)),]
pred1a.median <- apply(pred1a, 2, median )
pred1a.HPDI <-apply(pred1a, 2, HPDI )

d.pred2 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  grpsizeELA=rep(mean(d$s_grpsize),length(age.seq)),
  drought=age.seq,
  ibiELA=rep(mean(d$s_ibi),length(age.seq)),
  momloss=rep(mean(d$s_momlossadj),length(age.seq)),
  first=rep(mean(d$Primiparous),length(age.seq)),
  relrank=rep(1.5,length(age.seq))
)

LM1b <- link(carrying1_separate, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2b <- link(carrying2_separate, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)


PM1b <- (1-LM1b$p)*LM1b$mu
PM2b <- (1-LM2b$p)*LM2b$mu


pred_1b <- PM1b
pred_1b[1:idw[2],] <- PM2b[1:idw[2],]

pred1b <- pred_1b[sample(nrow(pred_1b)),]
pred1b.median <- apply(pred1b, 2, median )
pred1b.HPDI <-apply(pred1b, 2, HPDI )

#par(mfrow = c(1, 2), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))

plot( Bhminsrate[d$s_rank==-1.5] ~ s_drought[d$s_rank==-1.5] , data=dd , col=alpha("red",0.5), pch=16 ,ylim=c(0,.4),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_drought),max(d$s_drought)))
pred1a.lines=pred1a[sample(nrow(pred1a),100,replace=F),]
for (i in 1:100){
  lines( age.seq , pred1a.lines[i,] ,lwd=2,col=alpha("red",0.1))
}
lines( age.seq , pred1a.median ,lwd=1,col="black")
lines(age.seq, pred1a.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1a.HPDI[2,],lty=2,lwd=1)
axis(2,cex.axis=.9)
mtext( "biomass" , side=1 , line=.3, outer=FALSE , cex=1)
#mtext( "low rank" , side=3 , line=-1, outer=FALSE , cex=1)

plot( Bhminsrate[d$s_rank==1.5] ~ s_drought[d$s_rank==1.5] , data=dd , col=alpha("#33CCFF",0.5), pch=16 ,ylim=c(0,.4),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_drought),max(d$s_drought)))
pred1b.lines=pred1b[sample(nrow(pred1b),100,replace=F),]
for (i in 1:100){
  lines( age.seq , pred1b.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.1))
}
lines( age.seq , pred1b.median ,lwd=1,col="black")
lines(age.seq, pred1b.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1b.HPDI[2,],lty=2,lwd=1)

mtext( "biomass" , side=1 , line=.3, outer=FALSE , cex=1)
#mtext( "high rank" , side=3 , line=-1, outer=FALSE , cex=1)
mtext("proportion of time observed carrying", side=2, line=2,cex=1.5, outer=TRUE)


a_foc_z <- matrix(0,1000,length(unique(d$Focal)))
age.seq=seq(min(d$s_ibi),max(d$s_ibi),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  grpsizeELA=rep(mean(d$s_grpsize),length(age.seq)),
  drought=rep(mean(d$s_drought),length(age.seq)),
  ibiELA=age.seq,
  momloss=rep(mean(d$s_momlossadj),length(age.seq)),
  first=rep(mean(d$Primiparous),length(age.seq)),
  relrank=rep(-1.5,length(age.seq))
)

LM1a <- link(carrying1_separate, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2a <- link(carrying2_separate, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)


PM1a <- (1-LM1a$p)*LM1a$mu
PM2a <- (1-LM2a$p)*LM2a$mu


pred_1a <- PM1a
pred_1a[1:idw[2],] <- PM2a[1:idw[2],]

pred1a <- pred_1a[sample(nrow(pred_1a)),]
pred1a.median <- apply(pred1a, 2, median )
pred1a.HPDI <-apply(pred1a, 2, HPDI )

d.pred2 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  grpsizeELA=rep(mean(d$s_grpsize),length(age.seq)),
  drought=rep(mean(d$s_drought),length(age.seq)),
  ibiELA=age.seq,
  momloss=rep(mean(d$s_momlossadj),length(age.seq)),
  first=rep(mean(d$Primiparous),length(age.seq)),
  relrank=rep(1.5,length(age.seq))
)

LM1b <- link(carrying1_separate, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2b <- link(carrying2_separate, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)


PM1b <- (1-LM1b$p)*LM1b$mu
PM2b <- (1-LM2b$p)*LM2b$mu


pred_1b <- PM1b
pred_1b[1:idw[2],] <- PM2b[1:idw[2],]

pred1b <- pred_1b[sample(nrow(pred_1b)),]
pred1b.median <- apply(pred1b, 2, median )
pred1b.HPDI <-apply(pred1b, 2, HPDI )

#par(mfrow = c(1, 2), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))

plot( Bhminsrate[d$s_rank==-1.5] ~ s_ibi[d$s_rank==-1.5] , data=dd , col=alpha("red",0.5), pch=16 ,ylim=c(0,.4),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_ibi),max(d$s_ibi)))
pred1a.lines=pred1a[sample(nrow(pred1a),100,replace=F),]
for (i in 1:100){
  lines( age.seq , pred1a.lines[i,] ,lwd=2,col=alpha("red",0.1))
}
lines( age.seq , pred1a.median ,lwd=1,col="black")
lines(age.seq, pred1a.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1a.HPDI[2,],lty=2,lwd=1)
axis(2,cex.axis=.9)
mtext( "ibi" , side=1 , line=.3, outer=FALSE , cex=1)
#mtext( "low rank" , side=3 , line=-1, outer=FALSE , cex=1)

plot( Bhminsrate[d$s_rank==1.5] ~ s_ibi[d$s_rank==1.5] , data=dd , col=alpha("#33CCFF",0.5), pch=16 ,ylim=c(0,.4),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_ibi),max(d$s_ibi)))
pred1b.lines=pred1b[sample(nrow(pred1b),100,replace=F),]
for (i in 1:100){
  lines( age.seq , pred1b.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.1))
}
lines( age.seq , pred1b.median ,lwd=1,col="black")
lines(age.seq, pred1b.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1b.HPDI[2,],lty=2,lwd=1)

mtext( "ibi" , side=1 , line=.3, outer=FALSE , cex=1)
#mtext( "high rank" , side=3 , line=-1, outer=FALSE , cex=1)
mtext("proportion of time observed carrying", side=2, line=2,cex=1.5, outer=TRUE)


a_foc_z <- matrix(0,1000,length(unique(d$Focal)))
age.seq=seq(min(d$s_momlossadj),max(d$s_momlossadj),length=1000)

d.pred1 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  grpsizeELA=rep(mean(d$s_grpsize),length(age.seq)),
  drought=rep(mean(d$s_drought),length(age.seq)),
  ibiELA=rep(mean(d$s_ibi),length(age.seq)),
  momloss=age.seq,
  first=rep(mean(d$Primiparous),length(age.seq)),
  relrank=rep(-1,length(age.seq))
)

LM1a <- link(carrying1_separate, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2a <- link(carrying2_separate, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)


PM1a <- (1-LM1a$p)*LM1a$mu
PM2a <- (1-LM2a$p)*LM2a$mu


pred_1a <- PM1a
pred_1a[1:idw[2],] <- PM2a[1:idw[2],]

pred1a <- pred_1a[sample(nrow(pred_1a)),]
pred1a.median <- apply(pred1a, 2, median )
pred1a.HPDI <-apply(pred1a, 2, HPDI )

d.pred2 <- list(
  Focal=rep(1,length(age.seq)),
  female=rep(mean(d$female),length(age.seq)),
  challenge=rep(mean(d$s_challenges),length(age.seq)),
  grass=rep(mean(d$s_grass),length(age.seq)),
  grpsize=rep(mean(d$s_grpsize_cur),length(age.seq)),
  ageop=rep(mean(d$s_ageopuntia),length(age.seq)),
  age=rep(mean(d$s_age),length(age.seq)),
  momage=rep(mean(d$s_momage),length(age.seq)),
  grpsizeELA=rep(mean(d$s_grpsize),length(age.seq)),
  drought=rep(mean(d$s_drought),length(age.seq)),
  ibiELA=rep(mean(d$s_ibi),length(age.seq)),
  momloss=age.seq,
  first=rep(mean(d$Primiparous),length(age.seq)),
  relrank=rep(1.5,length(age.seq))
)

LM1b <- link(carrying1_separate, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2b <- link(carrying2_separate, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)


PM1b <- (1-LM1b$p)*LM1b$mu
PM2b <- (1-LM2b$p)*LM2b$mu


pred_1b <- PM1b
pred_1b[1:idw[2],] <- PM2b[1:idw[2],]

pred1b <- pred_1b[sample(nrow(pred_1b)),]
pred1b.median <- apply(pred1b, 2, median )
pred1b.HPDI <-apply(pred1b, 2, HPDI )

#par(mfrow = c(1, 2), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))

plot( Bhminsrate[d$s_rank==-1] ~ s_momlossadj[d$s_rank==-1] , data=dd , col=alpha("red",0.5), pch=16 ,ylim=c(0,.4),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_momlossadj),max(d$s_momlossadj)))
pred1a.lines=pred1a[sample(nrow(pred1a),100,replace=F),]
for (i in 1:100){
  lines( age.seq , pred1a.lines[i,] ,lwd=2,col=alpha("red",0.1))
}
lines( age.seq , pred1a.median ,lwd=1,col="black")
lines(age.seq, pred1a.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1a.HPDI[2,],lty=2,lwd=1)
axis(2,cex.axis=.9)
mtext( "age at loss of mom (<4yrs)" , side=1 , line=.3, outer=FALSE , cex=1)
#mtext( "low rank" , side=3 , line=-1, outer=FALSE , cex=1)

plot( Bhminsrate[d$s_rank==1.5] ~ s_momlossadj[d$s_rank==1.5] , data=dd , col=alpha("#33CCFF",0.5), pch=16 ,ylim=c(0,.4),xlab='', xaxt='n' ,yaxt='n', cex.lab=.1 , ylab='', xlim=c(min(d$s_momlossadj),max(d$s_momlossadj)))
pred1b.lines=pred1b[sample(nrow(pred1b),100,replace=F),]
for (i in 1:100){
  lines( age.seq , pred1b.lines[i,] ,lwd=2,col=alpha("#33CCFF",0.1))
}
lines( age.seq , pred1b.median ,lwd=1,col="black")
lines(age.seq, pred1b.HPDI[1,],lty=2,lwd=1)
lines(age.seq, pred1b.HPDI[2,],lty=2,lwd=1)

mtext( "age at loss of mom (<4yrs)" , side=1 , line=.3, outer=FALSE , cex=1)
#mtext( "high rank" , side=3 , line=-1, outer=FALSE , cex=1)
mtext("proportion of time observed carrying", side=2, line=2,cex=1.5, outer=TRUE)


a_foc_z <- matrix(0,1000,length(unique(d$Focal)))

d.pred1 <- list(
  Focal=c(1,1),
  female=rep(mean(d$female),2),
  challenge=rep(mean(d$s_challenges),2),
  grass=rep(mean(d$s_grass),2),
  grpsize=rep(mean(d$s_grpsize_cur),2),
  ageop=rep(mean(d$s_ageopuntia),2),
  age=rep(mean(d$s_age),2),
  momage=rep(mean(d$s_momage),2),
  grpsizeELA=rep(mean(d$s_grpsize),2),
  drought=rep(mean(d$s_drought),2),
  ibiELA=rep(mean(d$s_ibi),2),
  momloss=rep(mean(d$s_momlossadj),2),
  first=c(0,1),
  relrank=c(-1,-1)
)

LM1a <- link(carrying1_separate, n=1000 , data=d.pred1, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2a <- link( carrying2_separate, n=1000 , data=d.pred1, replace=
                list(am_Focal=a_foc_z), WAIC=TRUE)

w <- compare(carrying1_separate,carrying2_separate  , sort=FALSE)@output$weight
w <- w/sum(w)  ##make sure everything adds up to 1
idw <- round( w * 1000 )


PM1a <- (1-LM1a$p)*LM1a$mu
PM2a <- (1-LM2a$p)*LM2a$mu



pred_1a <- PM1a
pred_1a[1:idw[2],] <- PM2a[1:idw[2],]

pred1a <- pred_1a[sample(nrow(pred_1a)),]
pred1a.median <- apply(pred1a, 2, median )
pred1a.HPDI <-apply(pred1a, 2, HPDI )

d.pred2 <- list(
  Focal=c(1,1),
  female=rep(mean(d$female),2),
  challenge=rep(mean(d$s_challenges),2),
  grass=rep(mean(d$s_grass),2),
  grpsize=rep(mean(d$s_grpsize_cur),2),
  ageop=rep(mean(d$s_ageopuntia),2),
  age=rep(mean(d$s_age),2),
  momage=rep(mean(d$s_momage),2),
  grpsizeELA=rep(mean(d$s_grpsize),2),
  drought=rep(mean(d$s_drought),2),
  ibiELA=rep(mean(d$s_ibi),2),
  momloss=rep(mean(d$s_momlossadj),2),
  first=c(0,1),
  relrank=c(1,1)
)

LM1b <- link(carrying1_separate, n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)
LM2b <- link(carrying2_separate , n=1000 , data=d.pred2, replace=
               list(am_Focal=a_foc_z), WAIC=TRUE)

PM1b <- (1-LM1b$p)*LM1b$mu 
PM2b <- (1-LM2b$p)*LM2b$mu


pred_1b <- PM1b
pred_1b[1:idw[2],] <- PM2b[1:idw[2],]

pred1b <- pred_1b[sample(nrow(pred_1b)),]
pred1b.median <- apply(pred1b, 2, median )
pred1b.HPDI <-apply(pred1b, 2, HPDI )

#par(mfrow = c(1, 2), cex=1.1, mar=c(0,0,0,0), oma=c(3,3,3,3))
dens(pred1a[,1], xlim=c(0,.5), xlab="m", ylim=c(0,20), ylab="",col="white" , xaxt='n', yaxt='n')
ll <- d$Bhminsrate[d$Primiparous==0]
points(ll, rep(-.01,length(ll)), pch=15 , col=col.alpha("orange1", alpha=0.1) , cex=0.75)

shade( density(pred1a[,1]) , lim= as.vector(HPDI(pred1a[,1], prob=0.9999)) , col = col.alpha("orange1", 0.5))
shade( density(pred1a[,2]) , lim= as.vector(HPDI(pred1a[,2], prob=0.9999)) , col = col.alpha("#33CCFF", 0.5))
ll <- d$Bhminsrate[d$Primiparous==1]

points(ll, rep(-.03, length(ll)), pch=17 , col=col.alpha("#33CCFF", alpha=0.1) , cex=0.75)
axis(1, at = seq(from=0 , to=1, by = .5) ,tck=-0.02 , labels=T )
axis(1, at = seq(from=0 , to=1, by = .1) ,tck=-0.01 , labels=F)
abline(v=median(pred1a[,1]) , lty=1)
abline(v=median(pred1a[,2]) , lty=2)

legend(.25,20, legend = c("", ""),
       col=c(1,1)  , lty=c(1,2),
       lw=1 , cex=.85, bty="n")

legend(.3,20,, legend = c("multip", "primip"),
       col=c(col.alpha("orange1", 0.5) , col.alpha("#33CCFF", 0.5) ) , pch=c(15,17),
       pt.cex=2 , cex=.85, bty="n")


dens(pred1b[,1], xlim=c(0,.5), xlab="", ylim=c(0,20), ylab="",col="white" , xaxt='n', yaxt='n')
ll <- d$Bhminsrate[d$Primiparous==0]
points(ll, rep(-.01,length(ll)), pch=15 , col=col.alpha("orange1", alpha=0.1) , cex=0.75)

shade( density(pred1b[,1]) , lim= as.vector(HPDI(pred1b[,1], prob=0.9999)) , col = col.alpha("orange1", 0.5))
shade( density(pred1b[,2]) , lim= as.vector(HPDI(pred1b[,2], prob=0.9999)) , col = col.alpha("#33CCFF", 0.5))
ll <- d$Bhminsrate[d$Primiparous==1]

points(ll, rep(-.03, length(ll)), pch=17 , col=col.alpha("#33CCFF", alpha=0.1) , cex=0.75)
axis(1, at = seq(from=0 , to=1, by = .5) ,tck=-0.02 , labels=T )
axis(1, at = seq(from=0 , to=1, by = .1) ,tck=-0.01 , labels=F)

abline(v=median(pred1b[,1]) , lty=1)
abline(v=median(pred1b[,2]) , lty=2)

legend(.25,20, legend = c("", ""),
       col=c(1,1)  , lty=c(1,2),
       lw=1 , cex=.85, bty="n")

legend(.3,20,, legend = c("multip", "primip"),
       col=c(col.alpha("orange1", 0.5) , col.alpha("#33CCFF", 0.5) ) , pch=c(15,17),
       pt.cex=2 , cex=.85, bty="n")




