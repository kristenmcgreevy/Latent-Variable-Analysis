##--------------------------------------------------------------------------##
## Source Functions Edited after meeting 7/10
##--------------------------------------------------------------------------##


expit <- function(x){
  exp(x)/(1 + exp(x))
}

parse_miss<-function(dat){
  
  miss_var_names <- colnames(dat)[ apply(dat, 2, function(x) sum(is.na(x)))>0 ]
  miss_dat <- dat[, miss_var_names]
  dat$comb <- apply(miss_dat, 1, function(x) paste0(is.na(x), collapse = '-') ) 
  map <- unique(dat$comb)
  map2 <- unique(dat$comb)
  map <- 1:length(map) 
  names(map) <- map2
  dat$comb <- map[dat$comb] 
  split_dat <- split(dat, dat$comb) 
  dl <- lapply(split_dat, as.list) 
  ll <- do.call('c', dl) 
  names(ll) <- paste0(substring(text = names(ll),
                                first = as.numeric(gregexpr("\\.",names(ll)))+1, 
                                last=nchar(names(ll)) ),
                      substring(text = names(ll),
                                first = 1, 
                                last=as.numeric(gregexpr("\\.",names(ll)))-1 ) )
  
  ind<-unlist(lapply(ll,function(l) length(l)!=sum(is.na(l)) ))
  ll<-ll[ ind ]
  return(ll)
}

count_miss<-function(d){
  miss_FEV1<-mean(is.na(d$FEV1_FVC))
  return(mFEV1=miss_FEV1)
}



jagsmod_COPD <- "
model{
for(i in 1:sum(Npat)){

logit(pCOPD[i]) <- codeCOPD_b_int + codeCOPD_b*latentCOPD[i] # COPD 
# medications
logit(pSABA[i]) <- codeSABA_b_int + codeSABA_b_dm*latentCOPD[i] # Short-acting beta-agonist prescription
logit(pSAAC[i]) <- codeSAAC_b_int + codeSAAC_b_dm*latentCOPD[i] # Short-acting anticholinergics prescription
logit(pLABA[i]) <- codeLABA_b_int + codeLABA_b_dm*latentCOPD[i] # Long-acting beta-agonists (LABA) prescription
logit(pLAMA[i]) <- codeLAMA_b_int + codeLAMA_b_dm*latentCOPD[i] # Long-acting anticholinergic (LAMA) prescription
logit(pSCS[i]) <- codeSCS_b_int + codeSCS_b_dm*latentCOPD[i]    # Systemic corticosteroids prescription
logit(pLSCS[i]) <- codeLSCS_b_int + codeLSCS_b_dm*latentCOPD[i] # Long-term Systemic corticosteroids prescription
logit(pTHEO[i]) <- codeTHEO_b_int + codeTHEO_b_dm*latentCOPD[i] # Theopylline prescription
logit(pOXYG[i]) <- codeOXYG_b_int + codeOXYG_b_dm*latentCOPD[i] # Oxygen 
# comorbidities 
logit(pcharlson_q3_ind[i]) <- codecharlson_q3_ind_int + codecharlson_q3_ind_b*latentCOPD[i]   # charlson_q3_ind
logit(pcharlson_q2_ind[i]) <- codecharlson_q2_ind_int + codecharlson_q2_ind_b*latentCOPD[i]   # charlson_q2_ind
logit(pcharlson_q4_ind[i]) <- codecharlson_q4_ind_int + codecharlson_q4_ind_b*latentCOPD[i]   # charlson_q4_ind

COPD_p6[i] ~ dbern(pCOPD[i])
meds1[i] ~ dbern(pSABA[i])
meds2[i] ~ dbern(pSAAC[i])
meds3[i] ~ dbern(pLABA[i])
meds4[i] ~ dbern(pLAMA[i])
meds6[i] ~ dbern(pSCS[i])
meds7[i] ~ dbern(pLSCS[i])
meds8[i] ~ dbern(pTHEO[i])
oxygen[i] ~ dbern(pOXYG[i])
charlson_q3_ind[i] ~ dbern(pcharlson_q3_ind[i]) 
charlson_q2_ind[i] ~ dbern(pcharlson_q2_ind[i])
charlson_q4_ind[i] ~ dbern(pcharlson_q4_ind[i])


logit(pMiss_FEV[i]) <- a0 + a1*smoke_miss[i] + a2*smoke_ind[i] 
+ a_COPD_miss*latentCOPD[i] + a4*primary_63[i] + a5*primary_miss[i] 
+ a6*specialty_12[i] + a7*specialty_miss[i] + a8*ins_dual_ind[i]
+ a9*ins_medicaid_ind[i] + a10*ins_medicare_ind[i] + a11*hbpc[i]
+ a12*pneumonia[i] + a13*influenza[i] + a14*pulmonary_rehab[i] 
+ a15*hispanic_ind[i] + a17*age_69[i] 
+ a18*age_79[i] + a19*age_80plus[i] + a20*race_amind[i] 
+ a21*race_asian[i] + a22*race_black[i] + a23*race_hawai[i] + a24*race_multi[i] 
+ a25*race_unkno[i] 
+ a26*bmi_under[i] + a27*bmi_over[i] + a28*bmi_obese[i] + a29*bmi_miss[i] 
+ a30*depress_miss[i] + a31*hf_miss[i] + a32*pulhyper_miss[i] 
+ a33*diab_miss[i] + a34*asthma_miss[i] + a35*bronch_miss[i] 
+ a36*depress_val[i] + a37*hf_val[i] + a38*pulhyper_val[i] 
+ a39*diab_val[i] + a40*asthma_val[i] + a41*bronch_val[i]

FEV1_miss_ind[i] ~ dbern(pMiss_FEV[i])

logit(rho_COPD[i]) <- r0.5[i] + r1.5*sex[i] + r2.5*smoke_ind[i] + r3.5*smoke_miss[i]
+ r4.5*age_69[i] + r5.5*age_79[i] + r6.5*age_80plus[i] + r7.5*hispanic_ind[i]
+ r8.5*ins_dual_ind[i] + r9.5*ins_medicaid_ind[i] + r10.5*ins_medicare_ind[i] 
+ r11.5*exercise_inactive[i] + r12.5*exercise_insufficient[i] + r13.5*exercise_miss[i]
+ r15.5*influenza[i] + r16.5*pneumonia[i] 
+ r17.5*pulmonary_rehab[i] + r18.5*marriage[i] + r19.5*marriage_miss[i] 
+ r20.5*primary_63[i] + r21.5*primary_miss[i] + r22.5*ipc[i] + r23.5*oppc[i] 
+ r121.5*bmi_under[i] + r122.5*bmi_over[i] + r123.5*bmi_obese[i] + r124.5*bmi_miss[i] 
+ r157.5*diab_val[i] 
+ r158.5*depress_val[i] + r159.5*anxiety_val[i] + r160.5*chronic_pain_val[i] 
+ r161.5*asthma_val[i] + r162.5*bronch_val[i] + r163.5*snf_ind_val[i] 
+ r164.5*pulm_ind_val[i]
}

# Loop over subjects with non-missing FEV1 # 
for (i in 1:(Npat[1])){
FEV1_FVC[i] ~ dnorm(FEV1_b_int + FEV1_b_dm*latentCOPD[i], FEV1_tau)
}

codeCOPD_b_int ~ dunif(-4,-2)        
codeCOPD_b ~ dnorm(0,0.1)           
codeSAAC_b_int ~ dunif(-5,-2.2)      
codeSAAC_b_dm ~ dnorm(1.5,0.15)         
codeSABA_b_int ~ dunif(-0.62, 0.62)  
codeSABA_b_dm ~ dnorm(.5,0.1)          
codeLABA_b_int ~ dunif(-1.75, -1.1) 
codeLABA_b_dm ~ dnorm(.5,0.1)          
codeLAMA_b_int ~ dunif(-5, -2.2)    
codeLAMA_b_dm ~ dnorm(1.8,0.15)          
codeSCS_b_int ~ dunif(-5, -2.2)      
codeSCS_b_dm ~ dnorm(1.4,0.15)          
codeLSCS_b_int ~ dunif(-7, -2.2)     
codeLSCS_b_dm ~ dnorm(.5,0.1)         
codeTHEO_b_int ~ dunif(-7, -2.95)    
codeTHEO_b_dm ~ dnorm(.25,0.1)        
codeOXYG_b_int ~ dunif(-2.2, -1.1)    
codeOXYG_b_dm ~ dnorm(0,0.1)         
codecharlson_q3_ind_int ~ dunif(-1.39, -0.85)    
codecharlson_q3_ind_b ~ dnorm(0,0.1)         
codecharlson_q2_ind_int ~ dunif(-1.39, -0.85)     
codecharlson_q2_ind_b ~ dnorm(0,0.1)          
codecharlson_q4_ind_int ~ dunif(-1.39, -0.85)      
codecharlson_q4_ind_b ~ dnorm(0,0.1)              

FEV1_b_int ~ dnorm(70,10) 
FEV1_b_dm ~ dnorm(-0.2,1) 
FEV1_tau ~ dgamma(0.001,1000) 
FEV1_sigma <- 1/FEV1_tau    # because we have to specify information matrix

a0 ~ dnorm(0,0.1) 
a1 ~ dnorm(0,0.1)
a2 ~ dnorm(0,0.1)
a4 ~ dnorm(0,0.1)
a5 ~ dnorm(0,0.1)
a6 ~ dnorm(0,0.1)
a7 ~ dnorm(0,0.1)
a8 ~ dnorm(0,0.1)
a9 ~ dnorm(0,0.1)
a10 ~ dnorm(0,0.1)
a11 ~ dnorm(0,0.1)
a12 ~ dnorm(0,0.1)
a13 ~ dnorm(0,0.1)
a14 ~ dnorm(0,0.1)
a15 ~ dnorm(0,0.1)
a17 ~ dnorm(0,0.1)
a18 ~ dnorm(0,0.1)
a19 ~ dnorm(0,0.1)
a20 ~ dnorm(0,0.1)
a21 ~ dnorm(0,0.1)
a22 ~ dnorm(0,0.1)
a23 ~ dnorm(0,0.1)
a24 ~ dnorm(0,0.1)
a25 ~ dnorm(0,0.1)
a26 ~ dnorm(0,0.1)
a27 ~ dnorm(0,0.1)
a28 ~ dnorm(0,0.1)
a29 ~ dnorm(0,0.1)
a30 ~ dnorm(0,0.1)
a31 ~ dnorm(0,0.1)
a32 ~ dnorm(0,0.1)
a33 ~ dnorm(0,0.1)
a34 ~ dnorm(0,0.1)
a35 ~ dnorm(0,0.1)
a36 ~ dnorm(0,0.1)
a37 ~ dnorm(0,0.1)
a38 ~ dnorm(0,0.1)
a39 ~ dnorm(0,0.1)
a40 ~ dnorm(0,0.1)
a41 ~ dnorm(0,0.1)

a_COPD_miss ~ dnorm(-1,0.1) 

r1.5 ~ dnorm(0,10)
r2.5 ~ dnorm(0,10)
r3.5 ~ dnorm(0,10)
r4.5 ~ dnorm(0,10)
r5.5 ~ dnorm(0,10)
r6.5 ~ dnorm(0,10)
r7.5 ~ dnorm(0,10)
r8.5 ~ dnorm(0,10)
r9.5 ~ dnorm(0,10)
r10.5 ~ dnorm(0,10)
r11.5 ~ dnorm(0,10)
r12.5 ~ dnorm(0,10)
r13.5 ~ dnorm(0,10)
r15.5 ~ dnorm(0,10)
r16.5 ~ dnorm(0,10)
r17.5 ~ dnorm(0,10)
r18.5 ~ dnorm(0,10)
r19.5 ~ dnorm(0,10)
r20.5 ~ dnorm(0,10)
r21.5 ~ dnorm(0,10)
r22.5 ~ dnorm(0,10)
r23.5 ~ dnorm(0,10)
r121.5 ~ dnorm(0,10)
r122.5 ~ dnorm(0,10)
r123.5 ~ dnorm(0,10)
r124.5 ~ dnorm(0,10)
r157.5 ~ dnorm(0,10)
r158.5 ~ dnorm(0,10)
r159.5 ~ dnorm(0,10)
r160.5 ~ dnorm(0,10)
r161.5 ~ dnorm(0,10)
r162.5 ~ dnorm(0,10)
r163.5 ~ dnorm(0,10)
r164.5 ~ dnorm(0,10)


for (i in 1:sum(Npat)){
latentCOPD[i] ~ dbern(rho_COPD[i])
r0.5[i] ~ dunif(-7,-1)
}
}
"

run_bayes_COPD<-function(your.data, jagsmod, Burnin, sample, method, adapt, Monitor){
  temp <- your.data
  pedsdata <- temp[, c("meds1","meds2","meds3","meds4", "meds6", "meds7", 
                       "meds8", "oxygen", "ipc",
                       "specialty_12", "specialty_miss", "hbpc", "oppc",
                       "marriage", "marriage_miss", 
                       "charlson_q2_ind", "charlson_q3_ind", "charlson_q4_ind",
                       "pneumonia", "influenza", "pulmonary_rehab", "hispanic_ind",
                       "ins_dual_ind", "ins_medicaid_ind", "ins_medicare_ind",
                       "smoke_ind", "smoke_miss", "age_69", "age_79", "age_80plus", 
                       "exercise_inactive", "exercise_insufficient", "exercise_miss",
                       "primary_63", "primary_miss", "hispanic_miss",
                       "bmi_under", "bmi_over", "bmi_obese", "bmi_miss", 
                       "race_amind", "race_asian", "race_black", 
                       "race_hawai", "race_multi", "race_unkno", 
                       "hf_miss", 
                       "pulhyper_miss", "diab_miss", "depress_miss", 
                       "anxiety_miss", "chronic_pain_miss", "asthma_miss", 
                       "bronch_miss", "snf_ind_miss", "pulm_ind_miss", 
                       "hf_val", "pulhyper_val", "diab_val", 
                       "depress_val", "anxiety_val", "chronic_pain_val", 
                       "asthma_val", "bronch_val", "snf_ind_val", 
                       "pulm_ind_val"
  )]
  
  pedsdata$COPD_p6 <- ifelse(temp$phenotype_6 == 1, 1, 0)
  pedsdata$sex <- ifelse((!is.na(temp$SEX)) & (temp$SEX=="female"), 1, 0)
  pedsdata$FEV1_miss_ind <- ifelse(is.na(temp$FEV1_FVC)==TRUE, 1, 0)
  
  misspat <- paste(pedsdata$FEV1_miss_ind)
  pedsdata <- pedsdata[order(misspat),]
  outdata <- data.frame(pedsdata,T1D = temp[order(misspat),"FEV1_FVC"])
  
  dat_list <- list()
  dat_list$`Npat` <- as.vector(table(misspat))
  dat_list <- c(dat_list, as.list(pedsdata))
  
  results<-run.jags(jagsmod_COPD, 
                    monitor=monitor, 
                    data = dat_list,
                    burnin = Burnin, sample = sample, n.chains = 3,
                    method = 'rjags', adapt = adapt, silent.jags = TRUE,  inits=NA) 
  # Format MCMC samples
  results_mcmc<-as.mcmc.list(results)
  mymodel.mcmc <- as.mcmc(results)
  mymodel.mat <- as.matrix(mymodel.mcmc)
  mymodel.dat <- as.data.frame(mymodel.mat)
  
  output <- list(mcmc_list=results_mcmc, outdata = outdata) 
  
  return(output)
  return(mymodel.dat)
}


