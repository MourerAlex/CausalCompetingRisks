#install.packages("dplyr")
#install.packages("zoo")
#install.packages("Hmisc")


library(dplyr)
library(zoo)
library(Hmisc)

####
#### IMPORT DATA 
####
setwd/Users/stensrud/Documents/SPRINT_2020a/SPRINT-POP/data
getwd()
#Import data from the SPRINT trial ▸ ⁨data⁩
baseline <- read.csv("/Users/stensrud/Documents/SPRINT_2020a/SPRINT-POP/data/baseline.csv")
bp <- read.csv("/Users/stensrud/Documents/SPRINT_2020a/SPRINT-POP/data/bp.csv")
outcomes <- read.csv("/Users/stensrud/Documents/SPRINT_2020a/SPRINT-POP/data/outcomes.csv")
safety <- read.csv("/Users/stensrud/Documents/SPRINT_2020a/SPRINT-POP/data/safety.csv")

####
#### DATA CLEANING AND FORMATTING
####

head(bp)
# Convert to months as unit
days_month <- 365/12

# Find subjects 75 years or older
id_75 <- baseline$MASKID[baseline$AGE >= 75]
length(id_75)

# Time points (months) to include in the long format data frame 
cut_times <- 0:(12*6)

# Create a variable for mean arterial pressure (MAP). 
bp$MAP <- 1/3 * bp$SBP + 2/3 * bp$DBP;

# Create an integer denoting the month at which blood pressure is measured
bp$MONTH <- as.integer(sapply(strsplit(as.character(bp$VISITCODE), split='M', fixed=TRUE), function(x) (x[1])))
bp$MONTH[bp$VISITCODE=="RZ"] = 0
bp <- bp[order(bp$MASKID,bp$MONTH),]

# Create a discrete time indicator of the month at which death or / and AKI occur
safety$AKI_SAE_MONTH <- round(safety$AKI_SAE_DAYS/days_month)
outcomes$DEATH_MONTH <- round(outcomes$T_DEATH/days_month)

# Create a new data frame which includes the variables of interest (only)
compact_frame <- data.frame(baseline$MASKID,baseline$SMOKE_3CAT,baseline$SUB_CKD,baseline$STATIN,baseline$BMI,baseline$INTENSIVE,
                   baseline$AGE, baseline$SUB_CVD,baseline$FEMALE,
                   outcomes$DEATH_MONTH,outcomes$EVENT_DEATH,
                   safety$AKI_SAE_MONTH, safety$AKI_SAE_EVNT)

# Create appropriate column names
colnames(compact_frame) <- sapply(strsplit(as.character(colnames(compact_frame)), split='.', fixed=TRUE), function(x) (x[2]))

# Extra column needed for conversion to long format (to include first time interval)
compact_frame$Tstart <- -0.1 

# Create a smaller frame only including subjects 75 years or older with complete case measurements
elderly_frame <- na.omit(compact_frame[compact_frame$MASKID %in% id_75,])
# NA check: dim(elderly_frame); dim(na.omit(elderly_frame)). There are 49 subjects with missing values

# Ensure that subjects are at risk until latest event / censoring time (latest time of death and kidney failure)
elderly_frame$AKI_SAE_MONTH_C <- elderly_frame$AKI_SAE_MONTH 
elderly_frame$AKI_SAE_MONTH_C[elderly_frame$AKI_SAE_EVNT==0] <- apply(select(elderly_frame[elderly_frame$AKI_SAE_EVNT==0,],AKI_SAE_MONTH, DEATH_MONTH),1, max)
elderly_frame$CENSORING_EVENT <- elderly_frame$AKI_SAE_EVNT==0 & elderly_frame$EVENT_DEATH==0
elderly_frame$DEATH_BF_AKI_EVENT <- elderly_frame$AKI_SAE_EVNT==0 & elderly_frame$EVENT_DEATH==1

# Create a data frame in long format
long_frame <- survSplit(data=elderly_frame,cut=cut_times,start="Tstart",end="AKI_SAE_MONTH_C",event="AKI_SAE_EVNT")
# Add extra column to indicate censoring
long_frame_censoring <- survSplit(data=elderly_frame,cut=cut_times,start="Tstart",end="AKI_SAE_MONTH_C",event="CENSORING_EVENT")
# Overwrite CENSORING_EVENT in long_frame such that CENSORING_EVENT indicates time at which censoring occurs
long_frame$CENSORING_EVENT <- long_frame_censoring$CENSORING_EVENT

# Create an extra column to indicate death. Again, use max time of censoring and death
long_frame_death <- survSplit(data=elderly_frame,cut=cut_times,start="Tstart",end="AKI_SAE_MONTH_C",event="DEATH_BF_AKI_EVENT")
# Create an indicator of death BEFORE kidney failure in long_frame
long_frame$DEATH_BF_AKI_EVENT <- long_frame_death$DEATH_BF_AKI_EVENT
# Sanity check: sum(long_frame$DEATH_BF_AKI_EVENT)+sum(long_frame$AKI_SAE_EVNT)+sum(long_frame$CENSORING_EVENT)==dim(elderly_frame)[1]

# Match times in long_frame and bp, and only select subjects in id_75
matched_frame <- merge(long_frame[,colnames(long_frame) %in% c("MASKID","AKI_SAE_MONTH_C")], bp, by.x = c("MASKID", "AKI_SAE_MONTH_C"), by.y = c("MASKID", "MONTH") )
matched_frame <- matched_frame[order(matched_frame$MASKID,matched_frame$AKI_SAE_MONTH_C),]

# Create long_matched_frame, which includes blood pressure measurements in the long_frame. 
long_matched_frame <- left_join(long_frame,matched_frame,by=c("MASKID","AKI_SAE_MONTH_C"))
# Sanity check:  dim(long_matched_frame); dim(long_frame)
# Sanity check: Check that all subjects have blood pressure measurements at time 0
# Remark: 2 subjects do not have blood pressure at first visit. Set their BP to -1, but all of these subjects had a measurement after 1 month.
# sum(long_matched_frame$DEATH_MONTH==0 & long_matched_frame$AKI_SAE_MONTH_C==1 & is.na(long_matched_frame$MAP))
long_matched_frame$SBPcomplete <- long_matched_frame$SBP
long_matched_frame$SBPcomplete[long_matched_frame$AKI_SAE_MONTH_C==0 & is.na(long_matched_frame$SBP)] <- -1

long_matched_frame$DBPcomplete <- long_matched_frame$DBP
long_matched_frame$DBPcomplete[long_matched_frame$AKI_SAE_MONTH_C==0 & is.na(long_matched_frame$DBP)] <- -1


long_matched_frame$MAPcomplete <- long_matched_frame$MAP
long_matched_frame$MAPcomplete[long_matched_frame$AKI_SAE_MONTH_C==0 & is.na(long_matched_frame$MAP)] <- -1

# Carry forward blood pressure measurements until a new blood pressure measurement is reported.
long_matched_frame$SBPcomplete <- na.locf(long_matched_frame$SBPcomplete)
long_matched_frame$DBPcomplete <- na.locf(long_matched_frame$DBPcomplete)
long_matched_frame$MAPcomplete <- na.locf(long_matched_frame$MAPcomplete)

# Re-define missing initial values as NAs
long_matched_frame$MAPcomplete[long_matched_frame$MAPcomplete==(-1)] = NA
long_matched_frame$SBPcomplete[long_matched_frame$SBPcomplete==(-1)] = NA
long_matched_frame$DBPcomplete[long_matched_frame$DBPcomplete==(-1)] = NA

# Now, 2 subjects have NA's for BP the first visit
# Reduce frame by removing the subjects who doesn't have measurements of either DBP or SBP at first visit
long_matched_frame = long_matched_frame[!is.na(long_matched_frame$MAPcomplete),]

# Characteristics of the sample (excluded with comments)
# str(long_matched_frame[long_matched_frame$AKI_SAE_MONTH_C==0 & long_matched_frame$INTENSIVE==1,])
# str(long_matched_frame[long_matched_frame$AKI_SAE_MONTH_C==0 & long_matched_frame$INTENSIVE==0,])
# sum(elderly_frame$CENSORING_EVENT==1 & elderly_frame$AKI_SAE_MONTH_C<30)

# Create a binary MAP variable that is a candidate for L )have
long_matched_frame$low_MAP <- long_matched_frame$MAPcomplete < 75 #

# Fit pooled logistic model for the event of interest and the competing event
# using separate models for treatment INTENSIVE GROUP and NOT INTENSIVE GROUP. Rationale: Rely on less parametric assumptions
# Select follow-up times (which may be smaller than cut_times). Here we use 30 months
max_time <- 30 # months
follow_up_times <- 1:max_time

# Sanity: No subject with death and AKI same month sum(long_matched_frame$AKI_SAE_EVNT==1 & long_matched_frame$DEATH_BF_AKI_EVENT==1)

####
#### Here follows functions that are needed to obtain effect estimates
####
# Function to calculate a non-parametric outcome model
discrete_cuminc <- function(weight_vector, inputdata, grp=0, outcome_y=TRUE,follow_up=1:max_time){
  event_vec <- rep(NA, length.out=length(follow_up))
  counter <- 1 
  # count number of individuals in grp (that is, we cound those who were present at baseline)
  n_grp <- sum(inputdata$AKI_SAE_MONTH_C==0 & inputdata$INTENSIVE==grp)
  for(i in follow_up){
    if(outcome_y){
      indices <- inputdata$AKI_SAE_MONTH_C==i & inputdata$INTENSIVE == grp & inputdata$CENSORING_EVENT==0 & inputdata$DEATH_BF_AKI_EVENT==0 
      eventIndicator <- indices & inputdata$AKI_SAE_EVNT==1 
    }else{
      indices <- inputdata$AKI_SAE_MONTH_C==i & inputdata$INTENSIVE == grp & inputdata$CENSORING_EVENT==0 
      eventIndicator <- indices & inputdata$DEATH_BF_AKI_EVENT==1 
    }
    event_vec[counter] <- sum(weight_vector[eventIndicator]) / n_grp
    counter <- counter+1
  }
  output_cuminc <- cumsum(event_vec)
  return(output_cuminc)
}

##############
##### Imported and adapted code, which is similar to the one used in Stensrud et al, JASA, 2020. 
##############
non_para_cum_haz <- function(weightVector, inputdata, grp, outcomeKidney=TRUE){
  outputHazards <- rep(NA, length.out=length(follow_up_times))
  counter <- 1 
  for(i in follow_up_times){
    if(outcomeKidney){
      indices <- inputdata$AKI_SAE_MONTH_C==i & inputdata$INTENSIVE == grp & inputdata$CENSORING_EVENT==0 & inputdata$DEATH_BF_AKI_EVENT==0 
      eventIndicator <- indices & inputdata$AKI_SAE_EVNT==1 
    }else{
      indices <- inputdata$AKI_SAE_MONTH_C==i & inputdata$INTENSIVE == grp & inputdata$CENSORING_EVENT==0 
      eventIndicator <- indices & inputdata$DEATH_BF_AKI_EVENT==1 
    }
    outputHazards[counter] <- sum(weightVector[eventIndicator]) / sum(weightVector[indices])
    counter <- counter+1
  }
  return(outputHazards)
}

# Unchanged from pgform_ipwcs.R used in Young et al, Statistics in Medicine, 2020. 
nonParametricCumInc <- function(hazard1,hazard2,competing=FALSE){
  inc <- rep(NA, length.out=length(follow_up_times))
  cumulativeSurvival <- c(1, cumprod( (1-hazard1) * (1-hazard2) ))
  counter <- 1 
  for(i in 1:length(follow_up_times)){
    if(!competing){
      inc[i] <- hazard1[i] * (1-hazard2[i]) * cumulativeSurvival[i]
    }else{
      inc[i] <- hazard1[i] * cumulativeSurvival[i]
    }
  }
  cumInc <- cumsum(inc)
  return(cumInc)
}

####
####  The function bootSeparableEffectsLK fits pooled logistic models for discrete time hazard functions of interest
####
bootSeparableEffectsLK <- function(input_frame){
  plr_fit_AKI_TREAT <- glm(AKI_SAE_EVNT ~ AKI_SAE_MONTH_C + I(AKI_SAE_MONTH_C^2) + I(AKI_SAE_MONTH_C^3) 
                           + SMOKE_3CAT + SUB_CKD + STATIN 
                           + SUB_CVD + FEMALE + SBPcomplete + SBPcomplete^2 +  DBPcomplete +SBPcomplete*AKI_SAE_MONTH_C + SBPcomplete*SUB_CKD,
                           data = input_frame[input_frame$INTENSIVE==1 & input_frame$AKI_SAE_MONTH_C < max_time,], family = binomial())
  plr_fit_AKI_STD <- glm(AKI_SAE_EVNT ~ AKI_SAE_MONTH_C +  I(AKI_SAE_MONTH_C^2) + I(AKI_SAE_MONTH_C^3) 
                         + SMOKE_3CAT + SUB_CKD + STATIN 
                         + SUB_CVD + FEMALE + SBPcomplete + SBPcomplete^2 + DBPcomplete +SBPcomplete*AKI_SAE_MONTH_C + SBPcomplete*SUB_CKD,
                         data = input_frame[input_frame$INTENSIVE==0 & input_frame$AKI_SAE_MONTH_C < max_time,], family = binomial())
  plr_fit_CENS <- glm(CENSORING_EVENT ~ AKI_SAE_MONTH_C +  I(AKI_SAE_MONTH_C^2) + I(AKI_SAE_MONTH_C^3) 
                      + INTENSIVE + INTENSIVE*AKI_SAE_MONTH_C + SMOKE_3CAT + SUB_CKD + STATIN 
                      + AGE + SUB_CVD + FEMALE + SBPcomplete + DBPcomplete,
                      data = input_frame[input_frame$AKI_SAE_MONTH_C < max_time,], family = binomial())
  # Use the models to predict probabilities
  # Consider first the weigths based on modeling the event of interest (Y) and censoring
  # Remove to consider a common model
  pred_AKI_TREAT <- 1-predict(plr_fit_AKI_TREAT, newdata = input_frame, type = 'response')
  pred_AKI_STD <- 1-predict(plr_fit_AKI_STD, newdata = input_frame, type = 'response')
  pred_CENS <- 1-predict(plr_fit_CENS , newdata = input_frame, type = 'response')
  cum_pred_AKI_TREAT <- unlist(aggregate(pred_AKI_TREAT~input_frame$MASKID,FUN = cumprod)$pred_AKI_TREAT,use.names = F) #unlist because aggregate creates list
  cum_pred_AKI_STD <- unlist(aggregate(pred_AKI_STD~input_frame$MASKID,FUN = cumprod)$pred_AKI_STD,use.names = F) #unlist because aggregate creates list
  cum_pred_CENS <- unlist(aggregate(pred_CENS~input_frame$MASKID,FUN = cumprod)$pred_CENS,use.names = F) #unlist because aggregate creates list
  # Ready to create IPW: A_y=0,A_d=1
  # Create empty weight vector
  ipw_y <- rep(NA, length.out = length(pred_CENS))
  # index time zeros
  t0 <- input_frame$AKI_SAE_MONTH_C==0
  ipw_y[t0] <- (1-cum_pred_AKI_STD[t0])/(1-cum_pred_AKI_TREAT[t0])
  ipw_y[!t0] <- (1-pred_AKI_STD[!t0])/(1-pred_AKI_TREAT[!t0]) * 
    cum_pred_AKI_STD[which(!t0)-1]/cum_pred_AKI_TREAT[which(!t0)-1]
  ipw_cens <- 1/cum_pred_CENS
  ipw_sep_eff <- ipw_y/cum_pred_CENS
  ipw_sep_eff_inv <- 1/(ipw_y * cum_pred_CENS)
  # quantile(ipw_sep_eff, c(0.01, 0.99))# Sanity check of the weights. 
  # Calculate cumulative hazard
  cuminc_y0_d1 <- discrete_cuminc(ipw_sep_eff, input_frame,grp=1,follow_up = follow_up_times)
  cuminc_y1_d0 <- discrete_cuminc(ipw_sep_eff_inv, input_frame,grp=0,follow_up = follow_up_times)
  # Using the same models as for the sep eff. to estimate total eff.
  cuminc_y0_d0 <- discrete_cuminc(ipw_cens, input_frame,grp=0,follow_up = follow_up_times)
  cuminc_y1_d1 <- discrete_cuminc(ipw_cens, input_frame, grp=1,follow_up = follow_up_times)
  cuminc_d_y0_d0 <- discrete_cuminc(ipw_cens, input_frame,grp=0,follow_up = follow_up_times,outcome_y=FALSE)
  cuminc_d_y1_d1 <- discrete_cuminc(ipw_cens, input_frame, grp=1,follow_up = follow_up_times,outcome_y=FALSE)
  
  # Using weighted aalen johansen to estimate sep eff. 
  # first call the  function non_para_cum_haz,
  # which returns cumulative cause-sepcific hazards.
  cumhaz_y0_aj <- non_para_cum_haz(ipw_cens, inputdata=input_frame, grp=0, outcomeKidney=TRUE)
  cumhaz_d0_aj <- non_para_cum_haz(ipw_cens, inputdata=input_frame, grp=0, outcomeKidney=FALSE)
  cuminc_y0_d0_aj <- nonParametricCumInc(cumhaz_y0_aj,cumhaz_d0_aj)
  cumhaz_y1_aj <- non_para_cum_haz(ipw_cens, inputdata=input_frame, grp=1, outcomeKidney=TRUE)
  cumhaz_d1_aj <- non_para_cum_haz(ipw_cens, inputdata=input_frame, grp=1, outcomeKidney=FALSE)
  cuminc_y1_d1_aj <- nonParametricCumInc(cumhaz_y1_aj,cumhaz_d1_aj)
  return(data.frame(cuminc_y0_d0,cuminc_y1_d1,cuminc_y0_d1,cuminc_d_y0_d0,
                    cuminc_d_y1_d1,cuminc_y1_d0,cuminc_y0_d0_aj,cuminc_y1_d1_aj))
}

#####
##### Create bootstrap samples by calling the function bootSeparableEffectsLK
#####

#Fix seed for bootstrap. 
set.seed(123)

# Find point estimates from the full sample. 
estimates_full_data <- bootSeparableEffectsLK(long_matched_frame)

boot_samples <- 1 #499 to get uncertainty estimates.
n_cols <- ncol(estimates_full_data) # The number of estimands: number of columns in the output from RiskEstimatesFromBootstrapSample
n_rows <- length(follow_up_times) # The number of time points
subjects <- unique(long_matched_frame$MASKID)
n_subjects <- 1:length(subjects)

output_array <- array(NA,dim=c(n_rows, n_cols, boot_samples+1))
output_array[,,1] <- as.matrix(estimates_full_data)

# For loop that obtains estimates from each bootstrap sample and
# stores the output in a 3 dimensional array (time points of interests, parameter, bootstrap sample)
for(i in (1:boot_samples+1)){
  boot_index <- sample(subjects,size = length(subjects), replace=T)
  # Create a vector that contains indices for rows in long_matched_frame that 
  # are to be used in the bootstrap data. 
  list_indices <- sapply(boot_index, function(x) which(long_matched_frame$MASKID==x))
  dim_list_indices <- lengths(list_indices)
  # NEW_ID is a unique "bootstrap" individual identifier. 
  new_id <- rep(1:length(subjects),dim_list_indices)
  long_frame_indices <- unlist(list_indices)
  # Sanity check: 
  # str(unlist(long_frame_indices)); str(unique(unlist(long_frame_indices))) 
  boot_frame <- long_matched_frame[long_frame_indices,]
  boot_frame$MASKID = new_id
  boot_estimates <- as.matrix(bootSeparableEffectsLK(boot_frame)) #Make bootstrap estimates
  output_array[,,i] <- boot_estimates
}

# Save bootstrap samples: save(output_array,file = "SeparableEffectsBootstrapOutputAugustXX.Rdata")
# Load boostrap samples: load("SeparableEffectsBootstrapOutputAugustXX.Rdata")

time_of_interest <- which(follow_up_times==30) # 30 months

# Total effect on kidney failure
tot_effect_boot <- apply(output_array[time_of_interest,,],MARGIN=2,function(x) x[1]-x[2])
confidence_intervals_tot_eff <- quantile(tot_effect_boot,probs=c(0.025,0.975))

# Total effect on mortality
tot_effect_death_boot <- apply(output_array[time_of_interest,,],MARGIN=2,function(x) x[4]-x[5])
confidence_intervals_death_tot_eff <- quantile(tot_effect_death_boot,probs=c(0.025,0.975))

# Ay separable effect
ay_sep_effect_boot <- apply(output_array[time_of_interest,,],MARGIN=2,function(x) x[3]-x[2])
confidence_intervals_ay_eff <- quantile(ay_sep_effect_boot,probs=c(0.025,0.975))
head(ay_sep_effect_boot)

# Ad separable effect
ad_sep_effect_boot <- tot_effect_boot-ay_sep_effect_boot
confidence_intervals_ay_eff <- quantile(ad_sep_effect_boot,probs=c(0.025,0.975))
head(ad_sep_effect_boot)

# Obtain non-parametric confidence intervals from the bootstrap samples.
confidence_intervals <- apply(output_array[time_of_interest,,],MARGIN=1,quantile,probs=c(0.025,0.975))

# Find point estimates of various cumulative incidences
# _aj denotes 
cuminc_y0_d0 <- estimates_full_data[,1]
cuminc_y1_d1 <- estimates_full_data[,2]
cuminc_y0_d1 <- estimates_full_data[,3]
cuminc_d_y0_d0  <- estimates_full_data[,4]
cuminc_d_y1_d1 <- estimates_full_data[,5]
cuminc_y1_d0 <- estimates_full_data[,6]
cuminc_y0_d0_aj <- estimates_full_data[,7] 
cuminc_y1_d1_aj <- estimates_full_data[,8] 
head(estimates_full_data[,1])
total_effect <- cuminc_y1_d1-cuminc_y0_d0
ay_sep_effect <-  cuminc_y1_d1-cuminc_y0_d1 


x_times <- c(0,follow_up_times)
#setwd("/Users/stensrud/Dropbox/statistikk/r scripts/competing risks/aki sprint")
pdf("plot_cum_inc_both_events_2021.pdf")
plot(x_times, c(0,cuminc_y1_d1), type="s",ylim=c(0,0.1), ylab=expression('Cumulative incidence'), 
     xlab="Month (k+1)",lty=1,lwd=1,xlim=c(0,max_time), cex.lab =1.5, cex.axis = 1.5)
lines(x_times,c(0,cuminc_y0_d0),type="s", col=2,ylim=c(0,1),lwd=1)
lines(x_times,c(0,cuminc_d_y1_d1),type="s", lty=2,col=1,ylim=c(0,1),lwd=1)
lines(x_times,c(0,cuminc_d_y0_d0),type="s", lty=2,col=2,ylim=c(0,1),lwd=1)
legend(x=0,y=0.1, c(expression('AKI (Intensive, Pr('*'Y'['k+1']^{"a=1"}*'=1' *')'), expression('AKI (Standard, Pr('*'Y'['k+1']^{"a=0"}*'=1' *')'),expression('Death (Intensive, Pr('*'D'['k+1']^{"a=1"} *'=1'*')'), expression('Death (Standard, Pr('*'D'['k+1']^{"a=0"}*'=1' *')')),bty = "n",
       col=c(1,2,1,2), lty=c(1,1,2,2),cex=1.5,pt.cex=1.5,lwd=1.5)
dev.off()

pdf("plot_cum_inc_ay_2021.pdf")
par(mar=c(5.1,4.3,4.1,2.1))
plot(x_times,c(0,cuminc_y1_d1), type="s",ylim=c(0,0.1), ylab=expression('Cumulative incidence'), 
     xlab="Month (k+1)",lty=1,lwd=1,xlim=c(0,max_time), cex.lab =1.5, cex.axis = 1.5)
lines(x_times,c(0,cuminc_y0_d0), type="s",col=2,ylim=c(0,1),lwd=1)
lines(x_times,c(0,cuminc_y0_d1), type="s",col=3,ylim=c(0,1),lwd=1)
#lines(follow_up_times,cuminc_y1_d0, type="s",col=4,ylim=c(0,1),lwd=1)
legend(x=0,y=0.1, c(expression('AKI (Intenisve, Pr('*'Y'['k+1']^{"a"['Y']*'=1' * ' , ' * "a"['D']*'=1'} *'=1' * ')'* ')' ),
                    expression('AKI (Standard, Pr('*'Y'['k+1']^{"a"['Y']*'=0' * ' , ' * "a"['D']*'=0'} *'=1' * ')'* ')' ),
                    expression('AKI (Modified treatment, Pr('*'Y'['k+1']^{"a"['Y']*'=0' * ' , ' * "a"['D']*'=1'} *'=1' * ')'* ')')),bty = "n",
             col=c(1,2,3,3), lty=c(1,1,1,2),cex=1.5,pt.cex=1.5,lwd=1.5)
#legend(x=0,y=0.1, c("(Ay=1,Ad=1)","(Ay=0,Ad=0)","(Ay=0,Ad=1)"),bty = "n",
#      col=c(1,2,3,3), lty=c(1,1,1,2),cex=1,pt.cex=1,lwd=1)
dev.off()

pdf("effect_estim.pdf")
plot(total_effect, type="s",ylim=c(0,0.1), ylab="Difference in cumulative incidence", 
     xlab="Months",lty=1,lwd=1,xlim=c(0,max_time))
lines(ay_sep_effect, type="s",col=1,ylim=c(0,1),lwd=1,lty=2)
legend(x=0,y=0.1, c("Total effect",expression('A'['Y']*' separable effect')),bty = "n",
       col=c(1,1), lty=c(1,2),cex=1,pt.cex=1,lwd=1)
#legend(x=0,y=0.1, c("(Ay=1,Ad=1)","(Ay=0,Ad=0)","(Ay=0,Ad=1)"),bty = "n",
#      col=c(1,2,3,3), lty=c(1,1,1,2),cex=1,pt.cex=1,lwd=1)
dev.off()

