# Code by Robert Eugene del Carlo, Ph.D.
# This code takes skeletal muscle contraction data, in .csv format, to extract key performance indicators such as baseline tension, maximum tetanic force, tetanic duration, and first derivatives.
# This code was built to take data from the tetanus (sustained stimulus) protocol for data collected on an Aurora Scientific Myograph.
# TetanusP3 - Data analysis
#   Protocol 3: One 2sec pulse
#   1st stimulus occurs at 0.5s (seconds)
#   System delay is 4ms (miliseconds)
#   Pulsewidth is 500.0us (microseconds)
#   Stimulus rate is 1000 Hz
#   Last stimulus is at 2.5s
#
# Processing:
#   1. Pulse on if stimulus > 1.0
#   2. Base force = average of prior 1000 readings
#   3. Data of interest - 30,000 (3 secs) from start of pulse
#   4. Normalized force = (measuredF - baseF)*9.80665/(Muscle Mass in grams)
#   5. Added Rater

# To avoid errors, I define a function for evaluating minima
errorless_min <- function(x) {if (length(x)>0) min(x) else Inf}

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
dname = "~/Desktop/R_del_carlo_R_Sample_Code/Tetanus"
setwd(dname)

# Create a directory for report outputs
dir.create("./output/")
# Create a directory for graphical outputs
dir.create("./graphical_output/")

x = strsplit(dname, "/")

# Write / overwrite Column Headers for output summary file
ofsum = paste("./output/", x[[1]][lengths(x)],  ".csv", sep = "")
hdrs <-
  c(
    "Species", # Species of the snake, often of the genus Thamnophis; species include atratus, sirtalis, couchii, elegans and hammondii
    "GenoType", # Genotype referring to voltage gated sodium channel NaV1.4 that are either wildtype/TTX-sensitive/WT or mutant/TTX-resistant/P/EPN/V/VA/LVNV
    "Snake", # Snake ID, a combination of the collector ID and the animal number, e.g. CRF2345
    "Muscle", # the muscle that was extracted from the animal (i.e. the first muscle of the experiment is M1)
    "Rater", # This value describes the initials of the individual who performed the dissection and experiment (e.g. REdC = Robert Eugene del Carlo)
    "Mass(g)", # If available, this is the mass of the animal which donated a muscle for this experiment.
    "BaseF(N/g)",  # This value is calculated as the average force in the 1.5 seconds prior to the electrical stimulus eliciting a contractile response. It is recorded in Newtons per gram muscle mass.
    "ContrAmpl(N/g)", # This is the amplitude of the force produced after stimulus (i.e. the maximum force recorded within 1.5 seconds after the stimulus less the BaseF(N/g))
    "To10pct(ms)", # This is the time required to reach 10% maximal force in milliseconds
    "To50pct(ms)", # This is the time required to reach 50% relaxation *after* the maximum, recorded in milliseconds.
    "10to50pct(ms)", # This is To50pct(ms) less To10pct(ms), a proxy for contraction duration.
    "FMaxRateOfChg(N/g/s)", # This is the maximum rate of force development, that is the positive peak of the first derivative.
    "FMinRateOfChg(N/g/s)", # This is the maximum rate of relaxation following stimulation, that is the negative peak of the first derivative.
    "ToFChgMax(ms)", # This is the time between the stimulus and the positive peak of the first derivative.
    "ToFChgMin(ms)", # This is the time between the stimulus and the negative peak of the first derivative.
    "DiffFChgMaxToMin(ms)", # This is ToFChgMin(ms) less ToFChgMax(ms) - a proxy for contraction duration.
    "MAMU", # This is the 'Mass Adjusted Mouse Unit' - a measure of whole animal toxin resistance. 1 MAMU is the amount of TTX needed to kill a mouse. Consequently, a snake carrying 2 MAMU is resistant to twice the lethal mammalian dose.
    "County", # This is the county from whith the animal was sampled.
    "Longitude", # This is the longitude from which the animal was sampled.
    "Latitude", # This was the latitude from which the animal was sampled.
    "SVLmm", # This is the snout-vent length, a physical measure of the animals body that correlates to its chronological age.
    "MusMassg", # This is the mass of the skeletal muscle that generated the contractions recorded in such a file
    "BodMassg", # This is a duplicate column of Mass for ease of viewing.
    "Sex", # If available, this is the genital determined sex of the snake (M or F)
    "Date_Experimented", # This is the date of the experiment, recorded from lab notebooks
    "Date_Collected", # This is the date of the animals sampling from the environment.
    "Days_in_Between" # This is Date_Experimented less Date_Collected.
  )

# Create a table to house the processed data using the column headers above.
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
# Generate an output force file; definitions of columns correspond to those above.
ofF <- "./output/p3Tetanus-force.csv"
ofFhdr <-
  c(
    "Species",
    "GenoType",
    "Snake",
    "Muscle",
    "MAMU",
    "County",
    "Longitude",
    "Latitude",
    "SVLmm",
    "MusMassg",
    "BodMassg",
    "Sex",
    "Date_Experimented",
    "Date_Collected",
    "Days_In_Between",
    "Rater",
    format((1:30000) / 10000, scientific = FALSE)
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
# Create a file for the first derivative (1d) of the input force.
ofF1d <- "./output/p3Tetanus-force1d.csv"
ofF1dhdr <-
  c(
    "Species",
    "GenoType",
    "Snake",
    "Muscle",
    "MAMU",
    "County",
    "Longitude",
    "Latitude",
    "SVLmm",
    "MusMassg",
    "BodMassg",
    "Sex",
    "Date_Experimented",
    "Date_Collected",
    "Days_In_Between",
    "Rater",
    (1:30000) / 10000
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

fTet <- paste("./output/", "P3-TetanusFiles",  ".csv", sep = "")
write.table(
  t(c("File", "Snake", "Muscle", "Rater", "Mass(g)")),
  file = fTet,
  append = FALSE,
  quote = TRUE,
  sep = ",",
  eol = "\n",
  na = "NA",
  dec = ".",
  row.names = FALSE,
  col.names = FALSE
)

# This for loop makes a large conglomerate file that is used by later code to make the subsetted files.
for (file in files) {
  y <- strsplit(file, "-")
  snake <- trimws(y[[1]][1])
  # Using the SnakeID as a key, parsed from the filename, extract relevant information from the snake_info (sinf) and snake_muscle_mass (smm) dataframes.
  # Pull out the species information for the current snake.
  sSpecies <- as.character(sinf$Species[which(sinf$Snake == snake)])
  # Pull out the sodium channel genotype information for the current snake.
  sGenotype <- as.character(sinf$Genotype[which(sinf$Snake == snake)])
  # Pull out the organism-level TTX resistance information for the current snake.
  sMAMU <- as.character(sinf$MAMU[which(sinf$Snake == snake)])
  # Pull out the county sampling information for the current snake.
  sCounty <- as.character(sinf$County[which(sinf$Snake == snake)])
  # Pull out the Longitudinal information for the current snake.
  sLong <- as.character(sinf$Longitude[which(sinf$Snake == snake)])
  # Pull out the Latitudinal information for the current snake.
  sLat <- as.character(sinf$Latitude[which(sinf$Snake == snake)])
  # Pull out the Snout-Vent Length information for the current snake.
  sSVLmm <- as.character(sinf$SVLmm[which(sinf$Snake == snake)])
  # Pull out the Organismal Body Mass information for the current snake.
  sBodMassg <- as.character(sinf$BodMassg[which(sinf$Snake == snake)])
  # Pull out the gential sex information for the current snake.
  sSex <- as.character(sinf$Sex[which(sinf$Snake == snake)])
  # Pull out the date of the experiment for the current snake.
  sDtExp <-as.character(sinf$Date_Experimented[which(sinf$Snake == snake)])
  # Pull out the sampling date for the current snake.
  sDtColl <-as.character(sinf$Date_Collected[which(sinf$Snake == snake)])
  # Pull out the days in captivity information for the current snake.
  sDInBet <-as.character(sinf$Days_in_Between[which(sinf$Snake == snake)])
  # From the filename, parse out the muscle number (M#) describing which muscle was dissected from the current snake.
  muscle <- trimws(y[[1]][2])
  # Get Rater from Snake Info file for the current snake and current muscle
  rater <- as.character(sinf[which(sinf$Snake==snake), paste0("Rater..", muscle, ".")])
  # Get the mass of the current muscle
  muscMass <- smm[which(smm$SnakeID == snake), muscle] / 1000
  sMusMassg <- muscMass
  
  write.table(
    t(c(file, snake, muscle, rater, muscMass)),
    file = fTet,
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
  stmRows <- as.integer(rownames(raw[raw$Stimulus > 1.0, ]))
  if (length(stmRows) < 99) {
    cat("File: ", file, "; Stim Rows: ", length(stmRows), "\n")
    next
  }
  # Define the mean baseline force prior to the stimulus onset
  sr <- stmRows[1]
  meanF = mean(raw$Force..g.[sr - 1000:sr])
  rspF <-
    (raw$Force..g.[stmRows[1]:(stmRows[1] + 29999)] - meanF) * (0.00980665 / muscMass)
  oFline <- c(
    sSpecies,
    sGenotype,
    snake,
    muscle,
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
  # Smooth spline and take first derivative
  dt = (1:30000) / 10000
  spline <- smooth.spline(x= dt, y=rspF)
  prediction <- predict(spline)
  Force_Prime <- diff(rspF)/diff(dt)
  prediction <- predict(spline, deriv=1)
  rspF1d <- prediction$y
  oF1dline <- c(
    sSpecies,
    sGenotype,
    snake,
    muscle,
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
  tmf <- max(round(rspF, 0)) - 10
  tmmin = min(which(rspF >= tmf))
  tmmax <- max(which(rspF >= tmf))
  maxF <- mean(rspF[tmmin:tmmax])
  t100p <- round((tmmin + tmmax) / 2)
  # convert all times to ms from tenth of a ms
  t10p <- min(which(rspF > 0.1 * maxF))
  t50p <- t100p + min(which(rspF[t100p:length(rspF)] <= 0.5 * maxF))
  maxF1d <- max(rspF1d)
  minF1d <- min(rspF1d)
  t1dmax <- min(which(rspF1d >= maxF1d)) # * 50
  t1dmin <- min(which(rspF1d <= minF1d)) # * 50
  sumLine = c(
    sSpecies,
    sGenotype,
    snake,
    muscle,
	  rater,
    muscMass,
    meanF,
    maxF,
    t10p / 10,
    t50p / 10,
    (t50p - t10p) / 10,
    maxF1d,
    minF1d,
    t1dmax / 10,
    t1dmin / 10,
    (t1dmin - t1dmax) / 10,
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
    sDInBet
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
  
  # Create an image for each analysis processed for quick quality control.
  figname=paste0("./graphical_output/",snake,"-",muscle,".jpg")
  quartz(height=6, width=12)
  par(mar=c(5,5,1,1), mfrow=c(1,2))
  plot(spline,
       xlab="Time (ms)",
       ylab="Force (N/g-tissue)",
       cex.lab=1.25,
       las=1,
       pch=21,
       bg="gray",
       cex=1.5,
       lwd=1.4,
       plot=FALSE
  ) ; box(lwd=2) %>%
    title(main = snake, font.main= 4)
  plot(rspF1d,
       xlab="Time (ms)",
       ylab="Rate of Force (N/g-tissue/second)",
       cex.lab=1.25,
       las=1,
       pch=21,
       bg="gray",
       cex=0.5,
       lwd=1.4,
       plot=FALSE
  ) ; box(lwd=2)
  quartz.save(figname, type = "jpg")
  graphics.off()
  
} # End of Loops


df = read.csv("./output/p3tetanus-force.csv")
df <- df[order(df$GenoType, df$Species, df$Snake, df$Muscle),]
for (s in unique(df$Species)) {
  write.csv(df[which(df$Species == s), ], paste0("./output/p3Tetanus-Force-", s, ".csv"))
}

df = read.csv("./output/p3tetanus-force1d.csv")
df <- df[order(df$GenoType, df$Species, df$Snake, df$Muscle),]
for (s in unique(df$Species)) {
  write.csv(df[which(df$Species == s), ], paste0("./output/p3Tetanus-Force1d-", s, ".csv"))
}

setwd("~/Desktop/R_del_carlo_R_Sample_Code/Tetanus/output")

# For end users who intend to open this in Excel, this step transposes the above outputs.
# Excel has a limit on its number of columns, but rows are less limited.
syslist <- list.files(path = '.')
totranspose <- grep('p3Tetanus-', syslist, value = TRUE)
for (report in totranspose) {
  inputname <- report
  outputname <- paste0('transposed-',report)
  system(paste0("ruby -rcsv -e 'puts CSV.parse(STDIN).transpose.map &:to_csv' <", inputname, ">", outputname))
}

# Open the most important summary file from the analysis and use it to run preliminary visualization.
system('open ./Tetanus.csv')

# Provide some preliminary analysis by comparing each genotype within species using a boxplot.
summary_output = read.csv("./Tetanus.csv")
dataspecies <- unique(summary_output[,1])
for (each in dataspecies) {
  label=paste0("Thamnophis ",each)
  species_summary_subset <- summary_output[summary_output$Species == each, ]
  boxplot(log10(species_summary_subset$ContrAmpl.N.g.)~species_summary_subset$GenoType, xlab=label, ylab="Peak Force Ouput (Newtons / gram muscle)")
}