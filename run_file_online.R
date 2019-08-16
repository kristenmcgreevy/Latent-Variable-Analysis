##--------------------------------------------------------------------------##
## Apply Bayesian latent phenotyping model to COPD 
##
## UPDATED: by Kristen McGreevy 8/16/2019
##--------------------------------------------------------------------------##


library(runjags)
library(rjags)
library(beepr)
source("source_function_online.R") 

# if you want to test on a random subset of data
set.seed(3006) # random.org
index <- sample(1:nrow(your.data), 1000, replace=FALSE)
random_subset <- your.data[index, ]


##--------------------------------------------------------------------------##
## Run Bayesian Model
##--------------------------------------------------------------------------##

# may have to increase memory.limit in order for code to run. memory.limit(10000) seemed to work.

monitor <- c("rho_COPD","FEV1_b_int","FEV1_b_dm", "FEV1_tau", "codeCOPD_b_int", "codeCOPD_b",
             "codeSAAC_b_int", "codeSAAC_b_dm", "codeSABA_b_int", 
             "codeSABA_b_dm", "codeLABA_b_int", "codeLABA_b_dm", "codeLAMA_b_int", 
             "codeLAMA_b_dm", "codeSCS_b_int", 
             "codeSCS_b_dm", "codeLSCS_b_int", "codeLSCS_b_dm", "codeTHEO_b_int", 
             "codeTHEO_b_dm", "codeOXYG_b_int", "codeOXYG_b_dm", 
             "codecharlson_q3_ind_int", "codecharlson_q3_ind_b", 
             "codecharlson_q2_ind_int", "codecharlson_q2_ind_b", 
             "codecharlson_q4_ind_int", "codecharlson_q4_ind_b", "a_COPD_miss", 
             "a0", "a1", "a2", "a4", "a5", "a6", "a7", 
             "a8", "a9", "a10", "a11", "a12", "a13", "a14", "a15",
             "a17", "a18", "a19", "a20", "a21", "a22", "a23", "a24", "a25",
             "a26", "a27", "a28", "a29", "a30", "a31", "a32", "a33", "a34",
             "a35", "a36", "a37", "a38", "a39", "a40", "a41",
             "r1.5", "r2.5", 
             "r3.5", "r4.5", "r5.5", "r6.5", "r7.5", "r8.5", "r9.5", 
             "r10.5", "r11.5", "r12.5",
             "r13.5", "r15.5", "r16.5", "r17.5", "r18.5", "r19.5",
             "r20.5", "r21.5","r22.5","r23.5", 
             "r121.5", "r122.5", "r123.5", "r124.5", 
             "r157.5", "r158.5", "r159.5", "r160.5", 
             "r161.5", "r162.5", "r163.5", "r164.5"
)


start.time <- Sys.time()
COPD_bayes <- run_bayes_COPD(random_subset, jagsmod = jagsmod_COPD, Monitor = monitor, 
                                  Burnin = 2000, sample = 1000, method = 'simple', 
                                  adapt = 1000)
end.time <- Sys.time()
run.time <- end.time-start.time 
beep(4)

outdata <- COPD_bayes[[2]]
COPD_bayes2 <- COPD_bayes[[1]]

rho <- COPD_bayes2[[1]][,which(regexpr("rho_COPD",colnames(COPD_bayes2[[1]]))>0)]
parm <- COPD_bayes2[[1]][,which(regexpr("rho_COPD",colnames(COPD_bayes2[[1]]))<=0)]
parm[,c("a_COPD_miss")] <- exp(parm[,c("a_COPD_miss")])


postmeans <- apply(parm,2,mean)
postci   <-  apply(parm,2,quantile, probs = c(0.025, 0.975))

codesens <- apply(parm[, c("codeCOPD_b_int", "codeSAAC_b_int", "codeSABA_b_int", 
                               "codeLABA_b_int", "codeLAMA_b_int", 
                               "codeSCS_b_int", "codeLSCS_b_int", "codeTHEO_b_int", 
                               "codeOXYG_b_int", "codecharlson_q3_ind_int", 
                               "codecharlson_q2_ind_int", "codecharlson_q4_ind_int")] 
                  + parm[, c("codeCOPD_b", "codeSAAC_b_dm", "codeSABA_b_dm", 
                                 "codeLABA_b_dm", "codeLAMA_b_dm", 
                                 "codeSCS_b_dm", "codeLSCS_b_dm", "codeTHEO_b_dm", 
                                 "codeOXYG_b_dm", "codecharlson_q3_ind_b", 
                                 "codecharlson_q2_ind_b", "codecharlson_q4_ind_b")],
                  2, function(x){ mean(expit(x)) })
codesens.ci <- apply(parm[, c("codeCOPD_b_int", "codeSAAC_b_int", "codeSABA_b_int", 
                                  "codeLABA_b_int", "codeLAMA_b_int",  
                                  "codeSCS_b_int", "codeLSCS_b_int", "codeTHEO_b_int", 
                                  "codeOXYG_b_int", "codecharlson_q3_ind_int", 
                                  "codecharlson_q2_ind_int", "codecharlson_q4_ind_int")] 
                     + parm[, c("codeCOPD_b", "codeSAAC_b_dm", "codeSABA_b_dm", 
                                    "codeLABA_b_dm", "codeLAMA_b_dm", 
                                    "codeSCS_b_dm", "codeLSCS_b_dm", "codeTHEO_b_dm", 
                                    "codeOXYG_b_dm", "codecharlson_q3_ind_b", 
                                    "codecharlson_q2_ind_b", "codecharlson_q4_ind_b")],
                     2, function(x){quantile(expit(x), probs = c(0.025, 0.975))})
codespec <- apply(parm[,c("codeCOPD_b_int", "codeSAAC_b_int", "codeSABA_b_int", 
                              "codeLABA_b_int", "codeLAMA_b_int", 
                              "codeSCS_b_int", "codeLSCS_b_int", "codeTHEO_b_int", 
                              "codeOXYG_b_int", "codecharlson_q3_ind_int", 
                              "codecharlson_q2_ind_int", "codecharlson_q4_ind_int")],
                  2,function(x){mean(1-expit(x))})
codespec.ci <- apply(parm[,c("codeCOPD_b_int", "codeSAAC_b_int", "codeSABA_b_int", 
                                 "codeLABA_b_int", "codeLAMA_b_int",  
                                 "codeSCS_b_int", "codeLSCS_b_int", "codeTHEO_b_int", 
                                 "codeOXYG_b_int", "codecharlson_q3_ind_int", 
                                 "codecharlson_q2_ind_int", "codecharlson_q4_ind_int")],
                     2,function(x){quantile(1-expit(x), probs = c(0.025, 0.975))})

postCOPD <- apply(rho, 2, mean)
mean_overall <- mean(postCOPD) # estimated prevalence of COPD 
total_overall <- round(sum(postCOPD),0) # estimated number of people expected to have COPD


quantile_pred <- quantile(postCOPD, probs = mean_overall)
has_COPD <- ifelse(postCOPD <= quantile_pred, 1, 0)
latent_model_sensitivity = sum(has_COPD*random_subset$phenotype_6, na.rm=TRUE) / 
  sum(random_subset$phenotype_6, na.rm=TRUE)
latent_model_specificity = sum(has_COPD==0 & random_subset$phenotype_6==0, na.rm=TRUE) / 
  sum(random_subset$phenotype_6==0, na.rm=TRUE)

model_test <- apply_jags_model_sens(new_data=your.new.data, variable_labels=labels_table2, 
                                        training_data=random_subset, jags_model_output=COPD_bayes_run1)

library(ggplot2)

sens_1<- as.data.frame(t(codesens.ci))
sens_mean<- as.data.frame(codesens)

what <- as.data.frame(c("COPD", "SAAC", "SABA",  "LABA", "LAMA",  "SCS",  "LSCS", "THEO", 
                        "OXYG", "charlson_q3", "charlson_q2", "charlson_q4"), ncol=1)
sens_2 <- cbind(what, sens_1, sens_mean)
colnames(sens_2) <- c("parameters", "low", "high", "mean")
sens_2$parameters <- factor(sens_2$parameters, levels = sens_2$parameters[order(sens_2$mean)])

ggplot(sens_2, aes(x=parameters, y=mean)) +
  geom_pointrange(aes(ymin = low, ymax = high )) +
  labs(x="Parameters", y="Sensitivity") + theme_bw()


spec_1<- as.data.frame(t(codespec.ci))
spec_mean<- as.data.frame(codespec)

spec_2 <- cbind(what, spec_1, spec_mean)
colnames(spec_2) <- c("parameters", "low", "high", "mean")
spec_2$parameters <- factor(spec_2$parameters, levels = spec_2$parameters[order(spec_2$mean)])

ggplot(spec_2, aes(x=parameters, y=mean)) +
  geom_pointrange(aes(ymin = low, ymax = high )) +
  labs(x="Parameters", y="Specificity") + theme_bw()



library(MCMCvis)
library(dplyr)

## summary of estimated parameters arranged by convergence and magnitude of estimate. 
values<- MCMCsummary(COPD_bayes2, round = 3, excl = "rho_COPD")
small <- as.data.frame(values[, c("Rhat", "mean")])
colnames(small) <- c("Rhat", "mean")
param_names <- colnames(parm) %>% as.data.frame
colnames(param_names) <- "parameters"
small$magnitude <- abs(small$mean)
summary_table <- cbind(param_names, small) %>% arrange(Rhat, desc(magnitude))
summary_table[, c("parameters", "Rhat", "mean")]

P <- data.frame(
  parameters=c("a1", "a2",   "a3",   "a_COPD_miss", "a4", "a5",   "a6", "a7", "a8", 
               "a9", "a10", "a11",   "a12", "a13","a14", "a15", "a16", "a17", "a18", 
               "a19", "a20", "a21", "a22", "a23","a24", "a25", 
               "a26", "a27", "a28", "a29", "a30", "a31", "a32", "a33", "a34",
               "a35", "a36", "a37", "a38", "a39", "a40", "a41",
               "r1.5", "r2.5", 
               "r3.5","r4.5", "r5.5", "r6.5", "r7.5","r8.5", "r9.5", "r10.5", 
               "r11.5", "r12.5", "r13.5","r15.5", "r16.5", "r17.5", "r18.5", "r19.5", 
               "r20.5", "r21.5",  "r22.5", "r23.5", "r58.5","r121.5", "r122.5", 
               "r123.5","r124.5", "r125.5", "r126.5", "r127.5", "r128.5","r129.5", 
               "r130.5","r131.5", "r132.5", "r133.5",   "r134.5", "r135.5", "r136.5", 
               "r137.5", "r138.5",  "r139.5",   "r140.5","r141.5", "r142.5","r143.5", 
               "r144.5",  "r145.5", "r146.5","r147.5", "r148.5",  "r149.5", "r150.5", 
               "r151.5", "r152.5", "r153.5", "r154.5",  "r155.5", "r156.5", "r157.5",
               "r158.5", "r159.5", "r160.5", "r161.5", "r162.5", "r163.5","r164.5",
               "a0", "codecharlson_q2_ind_b", "codecharlson_q2_ind_int", "codecharlson_q3_ind_b", 
               "codecharlson_q3_ind_int", "codecharlson_q4_ind_b", "codecharlson_q4_ind_int", 
               "codeCOPD_b", "codeCOPD_b_int", "codeICS_b_dm", "codeICS_b_int", 
               "codeLABA_b_dm", "codeLABA_b_int", "codeLAMA_b_dm", "codeLAMA_b_int", 
               "codeLSCS_b_dm", "codeLSCS_b_int", "codeOXYG_b_dm", "codeOXYG_b_int", 
               "codeSAAC_b_dm", "codeSAAC_b_int", "codeSABA_b_dm", "codeSABA_b_int", 
               "codeSCS_b_dm", "codeSCS_b_int", "codeTHEO_b_dm", "codeTHEO_b_int", 
               "FEV1_b_dm", "FEV1_b_int"),
  Label =    c("Missing Smoking Status", "Smoker", "Sex Female","missing COPD Dx",
               "Primary Doc visits/year above avg",  "Missing Primary Doc visits",  
               "Specialty doctor visits/year above avg", "Missing Specialty Doc visits",  
               "Dual Insurance", "Medicaid Insurance", "Medicare Insurance", 
               "Home Based Palliative Care", "Pneumonia vaccination", "Influenza vaccination", 
               "Pulmonary rehabilitation", "Hispanic", "Missing Ethnicity",  "Age 60-69", 
               "Age 70-79", "Age 80 or over", "Race-American Indian",  "Race-Asian", "Race-Black", 
               "Race-Hawaiian or Pacific Islander", "Race-Mixed", "Race-Unknown", 
               "BMI-Underweight", "BMI-Overweight", "BMI-Obese", "BMI-Missing", 
               "Depression - Missing", "Heart Failure - Missing", "Pulmonary Hypertension - Missing",
               "Diabetes - Missing", "Asthma - Missing", "Bronchitis - Missing",
               "Depression", "Heart Failure", "Pulmonary Hypertension",  
               "Diabetes", "Asthma", "Bronchitis",
               "Sex Female", 
               "Smoker",  "Missing Smoking Status", "Age 60-69", "Age 70-79", "Age 80 or over", 
               "Hispanic", "Dual Insurance", "Medicaid Insurance", "Medicare Insurance",  
               "Exercise-inactive", "Exercise-insufficient", "Exercise-unknown", 
               "Influenza vaccination", "Pneumonia vaccination",  "Pulmonary rehabilitation", 
               "Partnered", "Partnering Unknown", "Primary Doc visits/year above avg",
               "Missing Primary Doc visits",  "Inpatient Palliative Care", 
               "Outpatient Palliative Care", "Missing Ethnicity", "BMI-Underweight", 
               "BMI-Overweight", "BMI-Obese", "BMI-Missing",  "Race-American Indian",  
               "Race-Asian", "Race-Black", "Race-Hawaiian or Pacific Islander", "Race-Mixed", 
               "Race-Unknown", "COPD Hospitalizations - Missing", "COPD Hospitalizations",  
               "COPD Emergency Department Visit - Missing",    "COPD Emergency Department Visit", 
               "COPD Observational visits - Missing",  "COPD Observational visits",  
               "Any Hospitalization - Missing", "Any Hospitalization",   
               "Any Emergency Department Visit - Missing",    "Any Emergency Department Visit", 
               "Any Observational visit - Missing",  "Any Observational visit", 
               "Any Urgent Care Visit - Missing", "Any Urgent Care Visit",   
               "Heart Failure - Missing", "Pulmonary Hypertension - Missing", "Diabetes - Missing", 
               "Depression - Missing",   "Anxiety - Missing", "Chronic pain - Missing", 
               "Asthma - Missing",  "Bronchitis - Missing", "SNF visit - Missing", 
               "Pulmonary visits - Missing", "Heart Failure", "Pulmonary Hypertension", 
               "Diabetes", "Depression", "Anxiety", "Chronic pain",  "Asthma", "Bronchitis", 
               "Any SNF Visit", "Any Pulmonary Visit",
               "Missing Model Intercept", "Charlson Q2 Slope", "Charlson Q2 Intercept", 
               "Charlson Q3 Slope", "Charlson Q3 Intercept", "Charlson Q4 Slope", 
               "Charlson Q4 Intercept", "Latent COPD Slope", "Latent COPD Intercept", 
               "ICS Slope", "ICS Intercept", "LABA slope", "LABA intercept", "LAMA slope", 
               "LAMA intercept", "LSCS slope", "LSCS intercept", "Oxygen slope", 
               "Oxygen intercept", "SAAC slope", "SAAC intercept", "SABA slope", "SABA intercept", 
               "SCS slope", "SCS intercept", "Theophylline slope", "Theophylline intercept", 
               "FEV1/FVC Missing slope", "FEV1/FVC Missing intercept" ))


# summary table with only the labels not coefficient names
sum_parm <- merge(P, summary_table, by = intersect(names(P), names(summary_table)),
                      all.y = TRUE, sort = TRUE, no.dups = TRUE)
sum <- sum_parm[, 2:5]
sum <- arrange(sum, sum$Rhat)

postmeans <- apply(parm,2,mean)
postci <- apply(parm,2,quantile, probs = c(0.025, 0.975))
postci <- as.data.frame(t(postci))
post_mean <- as.data.frame(postmeans)

mon <- as.data.frame(monitor, ncol=1)
monitor_norho <- mon[-c(mon[monitor=="rho_COPD",], mon[monitor=="FEV1_tau",]),]
post <- cbind(monitor_norho, postci, post_mean)
colnames(post) <- c("parameters", "low", "high", "mean")
post <- merge(P, post, by = intersect(names(P), names(post)),
                  all = TRUE, sort = TRUE, no.dups = TRUE)
post$parameters <- factor(post$parameters, levels = post$parameters[order(post$mean)])
post$overlaps_zero <- ifelse((post$low > 0 | post$high < 0),0,1) # max below zero, min above zero.


latent <- c("r1.5", "r2.5", 
            "r3.5","r4.5", "r5.5", "r6.5", "r7.5","r8.5", "r9.5", "r10.5", 
            "r11.5", "r12.5", "r13.5","r15.5", "r16.5", "r17.5", "r18.5", "r19.5", 
            "r20.5", "r21.5",  "r22.5", "r23.5", "r58.5","r121.5", "r122.5", 
            "r123.5","r124.5", "r125.5", "r126.5", "r127.5", "r128.5","r129.5", 
            "r130.5","r131.5", "r132.5", "r133.5",   "r134.5", "r135.5", "r136.5", 
            "r137.5", "r138.5",  "r139.5",   "r140.5","r141.5", "r142.5","r143.5", 
            "r144.5",  "r145.5", "r146.5","r147.5", "r148.5",  "r149.5", "r150.5", 
            "r151.5", "r152.5", "r153.5", "r154.5",  "r155.5", "r156.5", "r157.5",
            "r158.5", "r159.5", "r160.5", "r161.5", "r162.5", "r163.5","r164.5")
dependence <- c("FEV1_b_dm", "codeCOPD_b_int", "codeCOPD_b",
                "codeSAAC_b_int", "codeSAAC_b_dm", "codeSABA_b_int", 
                "codeSABA_b_dm", "codeLABA_b_int", "codeLABA_b_dm", "codeLAMA_b_int", 
                "codeLAMA_b_dm", "codeICS_b_int", "codeICS_b_dm", "codeSCS_b_int", 
                "codeSCS_b_dm", "codeLSCS_b_int", "codeLSCS_b_dm", "codeTHEO_b_int", 
                "codeTHEO_b_dm", "codeOXYG_b_int", "codeOXYG_b_dm", 
                "codecharlson_q3_ind_int", "codecharlson_q3_ind_b", 
                "codecharlson_q2_ind_int", "codecharlson_q2_ind_b", 
                "codecharlson_q4_ind_int", "codecharlson_q4_ind_b")
post$latent_var <- ifelse((post$parameters %in% latent),1,0)
latent_results <- post[post$parameters %in% latent,]
dependence_results <- post[post$parameters %in% dependence,]

missing_results <- post[!(post$parameters %in% latent) 
                            & !(post$parameters %in% dependence),]

# remove variables that were not in the analysis
included_latent <-latent_results[!is.na(latent_results$low),]
included_dependence <-dependence_results[!is.na(dependence_results$low),]
included_missing <-missing_results[!is.na(missing_results$low),]



library(ggplot2)
library(gplots)
today <- Sys.Date()
figurename <- paste(">location<Figure_Results_1000_obs_phen6def", today, ".pdf", sep="")
pdf(figurename) 

## only latent variables
ggplot(included_latent, aes(x=parameters, y=mean, color=factor(overlaps_zero))) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  ggtitle("Latent Model Variables \n Phenotype 6 COPD definition") +
  labs(x=" ", y="95% Credible Intervals") + theme_bw() +
  theme(axis.text.x = element_text(size=10), plot.title = element_text(lineheight=.8, hjust = 0.5))+
  coord_flip() +
  scale_x_discrete(labels = included_latent$Label) +
  scale_colour_manual(values=c("#1c9099", "#99CCCC"))

## only prior specified variables
ggplot(included_dependence, aes(x=parameters, y=mean, color=factor(overlaps_zero))) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  ggtitle("Dependent Factors \n Phenotype 6 COPD definition") +
  labs(x=" ", y="95% Credible Intervals") + theme_bw() +
  theme(plot.title = element_text(lineheight=.8, hjust = 0.5)) +
  coord_flip()  +
  scale_x_discrete(labels = included_dependence$Label) +
  scale_colour_manual(values=c("#453781FF", "#20A387FF"))

## only missing variables
ggplot(included_missing, aes(x=parameters, y=mean, color=factor(overlaps_zero))) +
  geom_pointrange(aes(ymin = low, ymax = high)) +
  ggtitle("Missing Model \n Phenotype 6 COPD definition") +
  labs(x=" ", y="95% Credible Intervals") + theme_bw() +
  theme(plot.title = element_text(lineheight=.8, hjust = 0.5)) +
  coord_flip()  +
  scale_x_discrete(labels = included_missing$Label) +
  scale_colour_manual(values=c("#453781FF", "#20A387FF"))

MCMCtrace(COPD_bayes2,  ISB = FALSE, excl = 'rho', pdf = FALSE)

dev.off() 









