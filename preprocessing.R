
# import data
data <- read.table("sam_1718_t_C.csv", header=TRUE, sep= ";", dec=",", na="-99")
names(data)[2] <- "tage.num"
names(data)

# select columns
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

# select intra-individual variables
cols_w <- c("Av1_state", "Iv1_state", "Uv1_state", "Co1_state", "Co2_state", "Leist_verstehen_state", "Leist_bearbeiten_state", 
              "Leist_stress_state", "Leist_ueberfordert_state", "Angst_abbruch_state", "Angst_scheitern_state", "PANP01_state", 
              "PANP05_state", "PANP08_state", "PANN01_state", "PANN05_state", "PANN09_state")
y_w <- data[, c("Code", "tage.num", cols_w, "event")]
dim(y_w) # 4063 x 20
y_w$Code <- as.factor(y_w$Code)
y_w$tage.num <- as.factor(y_w$tage.num)

Codes <- levels(y_w$Code)

# nDays (NOT equally spaced)
Nt <- length(levels(y_w$tage.num))
# n persons
N <- length(levels(y_w$Code))
# n variables
nVar <- length(c("tage.num", cols_w, "event"))


#######################################
# reshape the ts data in 3 dimensions #
#######################################

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

#################################################################
# remove persons whose event column contains non-binary entries #
#################################################################

yw[,,19][is.na(yw[,,19])] <- 99
yw[,,19][yw[,,19]==0.5] <- 1
Codes <- Codes[rowSums(yw[,,19]<0)==0]
yw <- yw[rowSums(yw[,,19]<0)==0,,]
yw[,,19][yw[,,19]==99] <- NA


N <- length(yw[,1,1])

##################################################
# 1. substitute post-dropout NA values with ones #
# 2. delete switch back                          #
##################################################
for (i in 1:N){
  for (t in 2:Nt){
    if (is.na(yw[i,t,19]))
      yw[i,t,19] <- yw[i,t-1,19]
  }
  
  for (t in 1:(Nt-1)){
    if (yw[i,t,19]>0)
      yw[i,t:(t+1),19] <- c(1,1)
  }
}