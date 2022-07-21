# Code by Robert Eugene del Carlo, Ph.D.
# This code processes a toxicological dose response assay of Tetrodotoxin on Skeletal Muscles dissected from Thamnophis garter snakes.
# Using the output from code03_TTX_dose_response_12_October_2020.R, specifically the file p4C4PwithTTX-ContrAmpl.N.g..csv, this code will run a Boltzmann Sigmoidal fitting routine and report parameters for dose response.

library(tidyverse)
library(minpack.lm)

setwd("~/Desktop/R_del_carlo_R_Sample_Code/")

# The dataframe containing information about individual snakes will be read in from the file named 'sinf' short for snake_info.
# The data contained in this file are organized in columns of the following categories:
# Snake/Species/Genotype/MAMU/County/Longitude/Latitude/SVLmm/BodMassg/Sex/Date_Experimented/Date_Collected/Days_in_Between/Rater_(M1)/Rater_(M2)/Rater_(M3)/Rater_(M4)/Rater_(M5)/Rater_(M6)/Rater_(M7)/Rater_(M8)/Rater_(M9)/Rater_(M10)/Rater_(M11)/Rater_(M12)/Rater_(M13)
sinf = read.csv("SnakeInfo-09.30.2020.csv")

# The dataframe containing the mass of the skeletal muscles used in this experiment will be read in from the file named 'smm' short for Snake_Muscle_Masses
smm = read.csv("SnakeSkeletalMuscleMasses-9.28.2020.csv")
# For experiments where a skeletal muscle mass is not recorded, apply a dummy value of 0.999 (grams)
# reset all missing muscle mass values to -1.0smm
smm[is.na(smm)] <- 0.999

# Specify the directory housing the input data, called dname short for directory_name.
dname = "~/Desktop/R_del_carlo_R_Sample_Code/TTX/output"

setwd(dname)

# Create a directory to house all the sigmoidal fitting outputs
dir.create("./sigmoidal")
# Create a directory for report outputs
dir.create("./sigmoidal/output/")
# Create a directory for graphical outputs
dir.create("./sigmoidal/graphical_output/")

smm_gathered <- smm %>% gather(Muscle, MusMassg, M1:M13) %>% select(-Date) %>% rename(Snake = SnakeID)

# Call the output file from code03 that contains the average maximal force from all the peaks in each file, organized in rows relative to dose of TTX in (nM)
df <- read.csv("./p4C4PwithTTX-ContrAmpl.N.g..csv", row.names = 2)
df <- df[1:nrow(df),2:ncol(df)]
df <- sweep(df,2,as.numeric(df[1,]),"/")
df <- df[,which(colSums(!is.na(df)) >= 4)]
head(df)

# Define a Boltzmann Sigmoidal Function.
# These variables are critical to describing the dose response of these experiments using a mathematical formuala.
# x is the variable, the concentration of TTX in nanomolar (nM) at the time of experiment.
# A1 is the amplitude of the response at 0 nM (control)
# A2 is the amplitude of the response at the maximal TTX concentration applied. It also serves as an offset value.
# x0, different than x, is the point of half maximal inhibition (that is, the point at which the toxin has inhibited half of the contractile force)
# dx is the factor describing the rate of change or slope between A1 and A2.

sigmoidal <- function(x,A1,A2,x0,dx){
  (A1-A2)/(1 + exp((log(x)-log(x0))/dx)) + A2
}
sigmoidal <- Vectorize(sigmoidal)

# This variation on the sigmoidal function, relies on an exponential formula. Depending on your preference, you can modify this component elsewhere in the code.
varsigmoidal <- function(x,x0,dx){
  (1)/(1 + exp((log(x)-log(x0))/dx))
}
varsigmoidal <- Vectorize(varsigmoidal)

# This residual function allows for an optimization of the fitting routine.
residFun <- function(parS,observed,indices){
  sigmoidal(as.numeric(rownames(df))[indices],parS$A1,parS$A2,parS$x0,parS$dx) - observed
}

# This residual function relies on the variant of the sigmoidal function to incorporate its variable sets.
varresidFun <- function(parS,observed,indices){
  varsigmoidal(as.numeric(rownames(df))[indices],parS$x0,parS$dx) - observed
}

# Appropriate starting values for intial stages of the fitting routine.
# The value of x0 is low enough to adequately capture TTX-sensitive animals and high enough to be useful to quickly fit TTX-resistant animals.
# Tailor the value of x0 based on an eyeballing of the TTX concentration at half-maximal inhibition for your own data.
parStart <- list(A2 = 0, A1 = 1, x0 = 450, dx = 2)
fitParams <- function(x){
  nls.out <- nls.lm(par = parStart, fn = varresidFun, control = nls.lm.control(maxiter = 100), observed = df[!is.na(df[,x]),x], indices = !is.na(df[,x]))
  unlist(nls.out$par[3:4])
}

result <- Vectorize(fitParams)(1:ncol(df))
colnames(result) <- colnames(df)

# The derivative of the sidmoidal function offers key insight to where a toxin begins to take effect and when its effects begin to saturate.
sigmoidalDeriv <- function(x,A1,A2,x0,dx){
  -((A1 - A2) * exp((log(x) + log(x0))/dx))/(dx * x * (exp(log(x)/dx) + exp(log(x0)/dx))^2)
}
sigmoidalDeriv <- Vectorize(sigmoidalDeriv)

max.abs <- function(x,...){
  sign(x[which.max(abs(x))])*max(abs(x),...)
}

# The maximum rate of decline once the toxin reaches a subclinical concentration may be proportional to the animal's overall resistance.
maxslope <- apply(result,2,function(p){max.abs(sigmoidalDeriv(as.numeric(rownames(df)),1,0,p[1],p[2]), na.rm = T)})

# Solve a sigmoidal function based on the toxin effective range, forcing it to start within the top 10% of the contaction and end in the bottom 10%.
solveSigmoidal <- function(y,A1,A2,x0,dx){
  exp(log(-(A2 - A1)/(y - A2) - 1) * dx + log(x0))
}

solveSigmoidal <- Vectorize(solveSigmoidal)
range.10.90 <- abs(solveSigmoidal(0.9,1,0,result[1,],result[2,]) - solveSigmoidal(0.1,1,0,result[1,],result[2,]))

result <- rbind(result,maxslope,range.10.90)
write.csv(t(result), "./sigmoidal/output/p4C4PwithTTX-amplitude-sigmoidal.csv")

# Read in the results file just created above.
rdf <- read.csv("./sigmoidal/output/p4C4PwithTTX-amplitude-sigmoidal.csv")
x <- strsplit(as.character(rdf$X),"_")
rdf$Snake <- unlist(lapply(x,'[[',1))
rdf$Muscle <- unlist(lapply(x,'[[',2))
rdf$Species <- sinf$Species[match(rdf$Snake, sinf$Snake)]
rdf$Genotype <- sinf$Genotype[match(rdf$Snake, sinf$Snake)]
rdf$MAMU <- sinf$MAMU[match(rdf$Snake, sinf$Snake)]
rdf$County <- sinf$COUNTY[match(rdf$Snake, sinf$Snake)]
rdf$Long <- sinf$Longitude[match(rdf$Snake, sinf$Snake)]
rdf$Lat <- sinf$Latitude[match(rdf$Snake, sinf$Snake)]
rdf$SVLmm <- as.character(sinf$SVLmm[match(rdf$Snake, sinf$Snake)])
rdf <- rdf %>% left_join(smm_gathered)
rdf$BodMassg <- as.character(sinf$BodMassg[match(rdf$Snake, sinf$Snake)])
rdf$Sex <- as.character(sinf$Sex[match(rdf$Snake, sinf$Snake)])
rdf$DtExp <- as.character(sinf$Date_Experimented[match(rdf$Snake, sinf$Snake)])
rdf$DtColl <- as.character(sinf$Date_Collected[match(rdf$Snake, sinf$Snake)])
rdf$DInBet <- as.character(sinf$Days_in_Between[match(rdf$Snake, sinf$Snake)])
write.csv(rdf, "./sigmoidal/output/TTX-Amplitude-Sigmoidal-rpt.csv")

# Open the most important summary file from the analysis and use it to run preliminary visualization.
system('open ./sigmoidal/output/TTX-Amplitude-Sigmoidal-rpt.csv')

# Provide some preliminary analysis by comparing each genotype within species using a boxplot.
summary_output = read.csv("./sigmoidal/output/TTX-Amplitude-Sigmoidal-rpt.csv")
dataspecies <- unique(summary_output[,9])
for (each in dataspecies) {
  label=paste0("Thamnophis ",each)
  species_summary_subset <- summary_output[summary_output$Species == each, ]
  figname=paste0("./sigmoidal/graphical_output/Thamnophis_",each,"-IC50.jpg")
  quartz(height=6, width=6)
  par(mar=c(5,5,1,1))
  boxplot(log10(species_summary_subset$x0)~species_summary_subset$Genotype, xlab=label, ylab="TTX IC50 (nM)")
  quartz.save(figname, type = "jpg")
  graphics.off()
}
