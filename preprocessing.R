# pi0 <- 0.95
# sit1: Nt = Nt + Nt2, N2 = N2
# scope: eta prediction and state prediction

data <- read.table("sam_1718_t_C.csv", header=TRUE, sep= ";", dec=",", na="-99")
names(data)[2] <- "tage.num"
names(data)

cols <- c("Code", "tage.num", "abi_note", "m_vertief", "m_vorbereit", "finanz_eltern", "sprache", "fachsem", "sex", "fw_pkt", "gesamt_iq", "Av", "Iv", "Se", "Uv", 
             "Co", "aist_r", "aist_i", "aist_a", "aist_s", "aist_e", "aist_c", "Ikont", "Ekont", "bfi_ex", "bfi_ve", "bfi_ge", "bfi_ne", "Av1_state", "Iv1_state", 
             "Uv1_state", "Co1_state", "Co2_state", "Pa", "PANP01_state", "PANP05_state", "PANP08_state", "Na", "PANN01_state", "PANN05_state", "PANN09_state", 
             "Angst_abbruch_state", "Angst_scheitern_state", "Leist_verstehen_state", "Leist_bearbeiten_state", "Leist_stress_state", "Leist_ueberfordert_state", 
             "Summe", "event", "Code", "Gruppe", "ProzentT1", "NoteT1", "ProzentT2", "NoteT2", "ProzentP", "NoteP", "EGKG")

dim(data)

#################################### 
# between-level observed variables #
####################################

################################### 
# within-level observed variables #
###################################

cols_w <- c("Av1_state", "Iv1_state", "Uv1_state", "Co1_state", "Co2_state", "Leist_verstehen_state", "Leist_bearbeiten_state", 
              "Leist_stress_state", "Leist_ueberfordert_state", "Angst_abbruch_state", "Angst_scheitern_state", "PANP01_state", 
              "PANP05_state", "PANP08_state", "PANN01_state", "PANN05_state", "PANN09_state")
y_w <- data[, c("Code", "tage.num", cols_w, "event")]
dim(y_w) # 4063 x 20
y_w$Code <- as.factor(y_w$Code)
y_w$tage.num <- as.factor(y_w$tage.num)

Codes <- levels(y_w$Code)
Nt <- length(levels(y_w$tage.num))
N <- length(levels(y_w$Code))
nVar <- length(c("tage.num", cols_w, "event"))

yw <- array(0, c(N, Nt, nVar)) 
dim(yw) # 122 x 52 x 19

count <- 0
for (i in 1:N){
  yw_i <- y_w[y_w$Code==levels(y_w$Code)[i],]
  yw[i,,1] <- as.numeric(levels(y_w$tage.num))
  

  for (t in 1:Nt){
    # if more than one response for a given day, average the responses
    if(length(yw_i$tage.num==levels(yw_i$tage.num)[t])==1 | length(yw_i$tage.num==levels(yw_i$tage.num)[t])==0){
      for (var in 2:nVar)
        yw[i,t,var] <- yw_i[yw_i$tage.num==levels(yw_i$tage.num)[t],c("tage.num", cols_w, "event")[var]]
      }
    else{
      for (var in 2:nVar)
        yw[i,t,var] <- mean(yw_i[yw_i$tage.num==levels(yw_i$tage.num)[t], c("tage.num", cols_w, "event")[var]],na.rm=T)
      count <- count + 1
    }
  }
}

yw[is.nan(yw)] <- NA

# remove persons whose event column contains non-binary entries 
yw[,,19][is.na(yw[,,19])] <- 99
yw[,,19][yw[,,19]==0.5] <- 1
Codes <- Codes[rowSums(yw[,,19]<0)==0]
yw <- yw[rowSums(yw[,,19]<0)==0,,]
yw[,,19][yw[,,19]==99] <- NA



# substitute NA values after the dropout with 1s 
N <- length(yw[,1,1])
for (i in 1:N){
  for (t in 2:Nt){
    if (is.na(yw[i,t,19]))
      yw[i,t,19] <- yw[i,t-1,19]
  }
  # delete switch back
  for (t in 1:(Nt-1)){
    if (yw[i,t,19]>0)
      yw[i,t:(t+1),19] <- c(1,1)
  }
}

yw_pool <- array(0, c(N * Nt, nVar)) 
Codes_pool <- list()
for (t in 1:Nt){
  yw_pool[(t-1)*N+(1:N),1:nVar] <- yw[1:N,t,1:nVar]
  Codes_pool[(t-1)*N+(1:N)] <- Codes[1:N]
}
dim(yw_pool)
length(Codes_pool)
summary(yw_pool)

# z_trans <- function(x)
#   (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) 
# yw_pool[,1:18] <- apply(yw_pool[,1:18], 2, z_trans)

colnames(yw_pool) <- colnames(y_w[2:20])
summary(yw_pool)

library(lavaan)

#######
# cfa #
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

fit_cfa <- cfa(model_cfa, data=yw_pool)
parTable(fit_cfa)
summary(fit_cfa, fit.measures=TRUE, standardized=TRUE)

eta_pool <- lavPredict(fit_cfa, method = "Bartlett")



##########
# lavaan #
##########

model_lavaan <- '

# latent variables
subImp =~ Av1_state + Iv1_state + Uv1_state
cost =~ Co1_state + Co2_state
understanding =~ Leist_verstehen_state + Leist_bearbeiten_state
stress =~ Leist_stress_state + Leist_ueberfordert_state
afraid =~ Angst_abbruch_state + Angst_scheitern_state
pa =~ PANP01_state + PANP05_state + PANP08_state
na =~ PANN01_state + PANN05_state + PANN09_state

# factor variances
subImp ~~ subImp
cost ~~ cost
understanding ~~ understanding
stress ~~ stress
afraid ~~ afraid
pa ~~ pa
na ~~ na
'

fit_lavaan <- lavaan(model_lavaan, data=yw_pool, auto.cov.lv.x=TRUE)
parTable(fit_lavaan)


