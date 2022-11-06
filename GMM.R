set.seed(7)

library(mclust)

# parameter options for GMM (shape of Gaussians) 
modelNames <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "VEI", "EEE", 
                "VEE", "EVE", "VVE", "EEV", "VEV", "EVV", "VVI", "VVV")

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
