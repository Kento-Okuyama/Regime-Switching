set.seed(7)

library(mclust)

# parameter options for GMM (shape of Gaussians) 
modelNames <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "VEI", "EEE", 
                "VEE", "EVE", "VVE", "EEV", "VEV", "EVV", "VVI", "VVV")

############################################
############################################
# Option1: apply GMM on observed variables #
############################################
############################################

#########################################
# Clustering applied by slicing by time #
#########################################

gmmCodes <- list()
gmm <- list()

# take all usable time points (34, 51-52 were not successfully fit)
ts <- c(2:33, 35:50) 

############################################
# apply GMM clustering for all time points #
############################################

for (t in ts) {
  x <- yw[complete.cases(yw[,t,]),t,1:18]
  gmmCodes[[t]] <- Codes[complete.cases(yw[,t,])]
  gmm[[t]] <- Mclust(x, G=2, modelNames=modelNames)
}

##########################################
##########################################
# Option2: apply GMM on latent variables #
##########################################
##########################################

############################################
# Confirmatory Factor Analysis with lavaan #
############################################

library(lavaan)

#######
# CFA #
#######

model_cfa <- '

# latent variables
subImp =~ Av1_state + Iv1_state + Uv1_state
cost =~ Co1_state + Co2_state
understanding =~ Leist_verstehen_state + Leist_bearbeiten_state
stress =~ Leist_stress_state + Leist_ueberfordert_state
afraid =~ Angst_abbruch_state + Angst_scheitern_state
pa =~ PANP01_state + PANP05_state + PANP08_state
na =~ PANN01_state + PANN05_state + PANN09_state

'

fit_cfas <- list()
etas <- list()

# take all usable time points (25, 34, 28, 34, 49, 51-52 were not successfully fit)
ts <- c(2:24, 26:27, 29:33, 35:48, 50) 

for (t in ts) {
  x <- yw[complete.cases(yw[,t,]),t,1:18]
  colnames(x) <- colnames(y_w[2:19])
  fit_cfas[[t]] <- cfa(model_cfa, data=x)
  etas[[t]] <- lavPredict(fit_cfa, method = "Bartlett")
  gmmCodes[[t]] <- Codes[complete.cases(yw[,t,])]
  gmm[[t]] <- Mclust(etas[[t]], G=2, modelNames=modelNames)
}

########################################################
########################################################
## What follows is not necessary for the main results ##
########################################################
########################################################

##################################################
# apply GMM clustering for a specific time point #
##################################################

# time at 10 (43 missing out of 117)
# sum(is.na(yw[,10,]))/17
x <- yw[complete.cases(yw[,10,]),10,1:18]
gmm_10 <- Mclust(x, G=2, modelNames=modelNames)

# summary(gmm_10)
# gmm_10$parameters$mean
# gmm_10$classification
# gmm_10$parameters$pro
# gmm_10$z
# gmm_10$BIC


# time at 20 (65 missing out of 117)
# sum(is.na(yw[,20,1:18]))/17
x <- yw[complete.cases(yw[,20,]),20,1:18]
gmm_20 <- Mclust(x, G=2, modelNames=modelNames)

# summary(gmm_20)
# gmm_20$parameters$mean
# gmm_20$classification
# gmm_20$BIC

# time at 30 (83 missing out of 117)
# sum(is.na(yw[,30,1:18]))/17
x <- yw[complete.cases(yw[,30,]),30,1:18]
gmm_30 <- Mclust(x, G=2, modelNames=modelNames)

# summary(gmm_30)
# gmm_30$parameters$mean
# gmm_30$classification
# gmm_30$BIC

# time at 40 (90 missing out of 117)
# sum(is.na(yw[,40,1:18]))/17
x <- yw[complete.cases(yw[,40,]),40,1:18]
gmm_40 <- Mclust(x, G=2, modelNames=modelNames)

# summary(gmm_40)
# gmm_40$parameters$mean
# gmm_40$classification
# gmm_40$BIC

################################################
# apply GMM clustering by pooling all the data #
################################################

# pooled (xx missing out of 6089)
# sum(is.na(yw_pool))/17
x <- yw_pool[complete.cases(yw_pool),2:18]
Codes_pool <- Codes_pool[complete.cases(yw_pool)]

gmm_pool <- Mclust(x, G=2, modelNames=modelNames)

summary(gmm_pool)
gmm_pool$parameters$mean
gmm_pool$classification
gmm_pool$BIC

#####################
# Check for NRIBE05 #
#####################

gmm_pool$classification[Codes_pool=="NRIBE05"] # the trajectory of the person throughout time
yw_pool[complete.cases(yw_pool),1:19][Codes_pool=="NRIBE05",c(1,19)] # the person did not drop out

#####################
# Check for AMACH11 #
#####################

gmm_pool$classification[Codes_pool== "AMACH11"] # the trajectory of the person throughout time
yw_pool[complete.cases(yw_pool),1:19][Codes_pool=="AMACH11",c(1,19)] # the person dropped out

#####################
# Check for SANSO03 #
#####################

gmm_pool$classification[Codes_pool== "SANSO03"] # the trajectory of the person throughout time
yw_pool[complete.cases(yw_pool),1:19][Codes_pool=="SANSO03",c(1,19)] # the person did not drop out

#####################
# Check for DJOJA05 #
#####################

gmm_pool$classification[Codes_pool== "DJOJA05"] # the trajectory of the person throughout time
yw_pool[complete.cases(yw_pool),1:19][Codes_pool=="DJOJA05",c(1,19)] # the person did not drop out, but not completed the survey

