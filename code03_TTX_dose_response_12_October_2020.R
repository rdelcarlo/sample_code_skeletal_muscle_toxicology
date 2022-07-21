# Code by Robert Eugene del Carlo, Ph.D.
# This code processes a toxicological dose response assay of Tetrodotoxin on Skeletal Muscles dissected from Thamnophis garter snakes.
# This code takes skeletal muscle contraction data, in .csv format, extracts the maximum value for each peak, and creates an average for the file.
# The average twitch force of in each file is a proxy for the effect of the toxin at that concentration, from 0nM (control) up to as high as 120,000nM.
#
# P4C4PwithTTX - Data analysis
#     Protocol 4: TTX Dose Response Curve
#     1st stimulus occurs at 5s (seconds)
#     System delay is 4ms (miliseconds)
#     Pulsewidth is 500.0us (microseconds)
#     Stimulus rate is 0.2 Hz
#     Last stimulus is at 20s
#
# Processing:
#   Contraction Amplitude (N/g)
#     This will then be plotted against the concentration of TTX that was on the muscle at that time.
# Once all amplitudes and concentration (for one individual/muscle) are plotted, a sigmoidal curve must be generated and the 4 values must be extracted:
#   The x-value at 50% y (x0)
#   The rate of change (dx)
#   The width of the curve from 10% to 90%
#   The peak first derivative of that sigmoidal
#
library(reshape2)

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
dname = "~/Desktop/R_del_carlo_R_Sample_Code/TTX"

setwd(dname)

# Create a directory for report outputs
dir.create("./output/")
# Create a directory for graphical outputs
dir.create("./graphical_output/")

x = strsplit(dname, "/")
# Write / overwrite Column Headers for output summary file
ofsum = paste("./output/", x[[1]][lengths(x)],  ".csv", sep = "")

hdrs <- c(
  "Seq",
  "File",
  "Snake",
  "Muscle",
  "Mass",
  "Dose",
  "Pulse",
  "BaseF(N/g)",
  "ContrAmpl(N/g)",
  "To100pct(ms)",
  "To10pct(ms)",
  "To50pct(ms)",
  "10to50pct(ms)",
  "FMaxRateOfChg(N/g/s)",
  "FMinRateOfChg(N/g/s)",
  "ToFChgMax(ms)",
  "ToFChgMin(ms)",
  "DiffFChgMaxToMin(ms)"
)
write.table(
  t(hdrs),
  file = ofsum,
  append = FALSE,
  quote = TRUE,
  sep = ",",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = FALSE,
  col.names = FALSE
)

# Output Force File
ofF <- "./output/p4C4PTTX-force.csv"
ofFhdr <- c(
  "Species",
  "GenoType",
  "MAMU",
  "COUNTY",
  "LONGITUDE",
  "LATITUDE",
  "SVLmm",
  "MusMassg",
  "BodMassg",
  "Sex",
  "DATE_Experimented",
  "Date_Collected",
  "Days_In_Between",
  "Pulse",
  "Snake",
  "Muscle",
  "Dose(nM)",
  "Rater",
  format((1:1500) / 10000, scientific = FALSE)
)
write.table(
  t(ofFhdr),
  file = ofF,
  append = FALSE,
  quote = TRUE,
  sep = ",",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = FALSE,
  col.names = FALSE
)

# Output Force first derivative File
ofF1d <- "./output/p4c4pTTX-force1d.csv"
ofF1dhdr <- c(
  "Species",
  "GenoType",
  "MAMU",
  "COUNTY",
  "LONGITUDE",
  "LATITUDE",
  "SVLmm",
  "MusMassg",
  "BodMassg",
  "Sex",
  "DATE_Experimented",
  "Date_Collected",
  "Days_In_Between",
  "Pulse",
  "Snake",
  "Muscle",
  "Dose(nM)",
  "Rater",
  format((1:29) / 200, scientific = FALSE)
)
write.table(
  t(ofF1dhdr),
  file = ofF1d,
  append = FALSE,
  quote = TRUE,
  sep = ",",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = FALSE,
  col.names = FALSE
)

# Identify all data files
files = list.files(path = ".", pattern = "csv")
q <- strsplit(files, "-")

fTTX <- paste("./output/", "P4C4PwithTTXFiles",  ".csv", sep = "")
write.table(
  t(c(
    "FileName", "Snake", "Dose(nM)", "Muscle", "Rater", "Mass(g)"
  )),
  file = fTTX,
  append = FALSE,
  quote = TRUE,
  sep = ",",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = FALSE,
  col.names = FALSE
)

fnum <- 0
for (file in files) {
  fnum <- fnum + 1
  y <- strsplit(file, "-")
  # Define the SnakeID from the filename.
  # Pull out the species information for the current snake.
  snake <- trimws(y[[1]][1])
  # From the filename, parse out the muscle number (M#) describing which muscle was dissected from the current snake.
  muscle <- trimws(y[[1]][2])
  # From the filename, parse out the dose of TTX (in nanomolar) and subsequently process it to retrieve a numerical value.
  dose <- trimws(y[[1]][3])
  dose <- sub(".dat.csv", "", dose)
  dose <- as.numeric(gsub("\\D", "", dose))
  dose <- (dose - 1) / 10
  # Using the SnakeID as a key, parsed from the filename, extract relevant information from the snake_info (sinf) and snake_muscle_mass (smm) dataframes.
  sSpecies <- as.character(sinf$Species[which(sinf$Snake == snake)])
  # Pull out the sodium channel genotype information for the current snake.
  sGenotype <- as.character(sinf$Genotype[which(sinf$Snake == snake)])
  # Pull out the organism-level TTX resistance information for the current snake.
  sMAMU <- as.character(sinf$MAMU[which(sinf$Snake == snake)])
  # Pull out the county sampling information for the current snake.
  sCounty <- as.character(sinf$COUNTY[which(sinf$Snake == snake)])
  # Pull out the Longitudinal information for the current snake.
  sLong <- as.character(sinf$LONGITUDE[which(sinf$Snake == snake)])
  # Pull out the Latitudinal information for the current snake.
  sLat <- as.character(sinf$LATITUDE[which(sinf$Snake == snake)])
  # Pull out the Snout-Vent Length information for the current snake.
  sSVLmm <- as.character(sinf$SVLmm[which(sinf$Snake == snake)])
  # Pull out the Organismal Body Mass information for the current snake.
  sBodMassg <- as.character(sinf$BodMassg[which(sinf$Snake == snake)])
  # Pull out the gential sex information for the current snake.
  sSex <- as.character(sinf$Sex[which(sinf$Snake == snake)])
  # Pull out the date of the experiment for the current snake.
  sDtExp <- as.character(sinf$DATE.Experimented[which(sinf$Snake == snake)])
  # Pull out the sampling date for the current snake.
  sDtColl <- as.character(sinf$Date.Collected[which(sinf$Snake == snake)])
  # Pull out the days in captivity information for the current snake.
  sDInBet <- as.character(sinf$Days.in.Between[which(sinf$Snake == snake)])
  # Get the mass of the current muscle
  muscMass <- smm[which(smm$SnakeID == snake), muscle] / 1000
  sMusMassg <- muscMass
  # Get Rater from Snake Info file for the current snake and current muscle
  rater <- as.character(sinf[which(sinf$Snake==snake), paste0("Rater..", muscle, ".")])

  write.table(
    t(c(file, snake, dose, muscle, rater, muscMass)),
    file = fTTX,
    append = TRUE,
    quote = TRUE,
    sep = ",",
    eol = "\n",
    na = "NA",
    dec = ".",
    row.names = FALSE,
    col.names = FALSE
  )
  
  ofdata <- paste("./output/", snake, "-", muscle, ".csv", sep = "")
  raw = read.csv(file)
  # Select all rows with stimulus > 1
  stmRows <- which(raw$Stimulus > 1.0)
  if (length(stmRows) < 4) {
    cat("File: ", file, "; Stim Rows: ", length(stmRows), "\n")
    next
  }
  # Find the last value of each stimulus pulse to establish the beginning of the contraction.
  categories <- cutree(hclust(dist(stmRows)), k = 4)
  stmRows <- aggregate(stmRows, list(cats = categories), "max")[, 2]
  if (length(stmRows) < 4) {
    cat("File: ", file, "; Stim Rows: ", length(stmRows), "\n")
    next
  }
  
  # Define the mean baseline force prior to the stimulus onset
  i <- 0
  for (sr in stmRows) {
    i <- i + 1
    meanF = mean(raw$Force..g.[sr - 1000:sr])
    rspF <-
      (raw$Force..g.[stmRows[1]:(stmRows[1] + 1499)] - meanF) * (0.00980665 / muscMass)
    oFline <- c(
      sSpecies,
      sGenotype,
      sMAMU,
      sCounty,
      sLong,
      sLat,
      sSVLmm,
      sMusMassg,
      sBodMassg,
      sSex,
      sDtExp,
      sDtColl,
      sDInBet,
      i,
      snake,
      muscle,
      dose,
	  rater,
      rspF
    )
    write.table(
      t(oFline),
      file = ofF,
      append = TRUE,
      quote = TRUE,
      sep = ",",
      eol = "\n",
      na = "NA",
      dec = ".",
      row.names = FALSE,
      col.names = FALSE
    )
    #  Low pass 200Hz filter -> pick every 30th entry; for computational expediency provided minimal graphing.
    rspF1d <- diff(rspF[1:30 * 50]) * 200
    oF1dline <- c(
      sSpecies,
      sGenotype,
      sMAMU,
      sCounty,
      sLong,
      sLat,
      sSVLmm,
      sMusMassg,
      sBodMassg,
      sSex,
      sDtExp,
      sDtColl,
      sDInBet,
      i,
      snake,
      muscle,
      dose,
	  rater,
      rspF1d
    )
    write.table(
      t(oF1dline),
      file = ofF1d,
      append = TRUE,
      quote = TRUE,
      sep = ",",
      eol = "\n",
      na = "NA",
      dec = ".",
      row.names = FALSE,
      col.names = FALSE
    )
    maxF <- max(rspF)
    t100p <- min(which(rspF >= maxF))
    # convert all times to ms from tenth of a ms
    t10p <- min(which(rspF > 0.1 * maxF))
    t50p <- t100p + min(which(rspF[t100p:length(rspF)] <= 0.5 * maxF))
    maxF1d <- max(rspF1d)
    minF1d <- min(rspF1d)
    t1dmax <- min(which(rspF1d >= maxF1d)) * 30
    t1dmin <- min(which(rspF1d <= minF1d)) * 30
    sumLine = c(
      fnum,
      file,
      snake,
      muscle,
      muscMass,
      dose,
      i,
      meanF,
      maxF,
      t100p / 10,
      t10p / 10,
      t50p / 10,
      (t50p - t10p) / 10,
      maxF1d,
      minF1d,
      t1dmax / 10,
      t1dmin / 10,
      (t1dmin - t1dmax) / 10
    )

    write.table(
      t(sumLine),
      file = ofsum,
      append = TRUE,
      quote = TRUE,
      sep = ",",
      eol = "\n",
      na = "NA",
      dec = ".",
      row.names = FALSE,
      col.names = FALSE
    )
  }
}

df = read.csv("./output/p4C4PTTX-force.csv")

write.csv(t(df[order(df$Pulse, df$Snake, df$Muscle, df$Dose), ]), "./output/p4C4PTTX-force-rpt.csv")

df = read.csv("./output/p4C4PTTX-force1d.csv")
write.csv(t(df[order(df$Pulse, df$Snake, df$Muscle, df$Dose), ]), "./output/p4C4PTTX-force1d-rpt.csv")

df = read.csv(ofsum)
agg <-
  aggregate(
    df,
    by = list(
      Dose = df$Dose,
      Muscle = df$Muscle,
      Snake = df$Snake
      ,       File = df$File
    ),
    FUN = mean,
    simplify = TRUE
  )
a1 <- agg[, c(4, 3, 2, 9:10, 12:22)]
write.csv(a1, paste("./output/", x[[1]][lengths(x)],  "-rpt.csv",sep=""))

# Create reports for each parameter by dose
library(reshape2)
a1_cast <- list()
for (var in colnames(a1)[6:ncol(a1)]) {
  a1_cast[[var]] <-
    dcast(a1,
          Dose ~ Snake + Muscle,
          value.var = var,
          fun.aggregate = mean)
}

for (var in names(a1_cast)) {
  write.csv(a1_cast[[var]], paste0("./output/p4C4PwithTTX-", var, ".csv"))
}

# Open the most important summary file from the analysis for quality control checking by hand.
system('open ./p4C4PwithTTX-ContrAmpl.N.g..csv')
