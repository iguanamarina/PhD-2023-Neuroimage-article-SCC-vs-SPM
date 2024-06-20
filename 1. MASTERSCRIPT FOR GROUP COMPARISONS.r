############################## ################### ############## ### 
##
## Script name: 1. MASTERSCRIPT FOR GROUP COMPARISONS.r
##
## Purpose of script: Use SCCs for the estimation of regions with different
##                    between-group activity levels in PET images with simulated ROIs 
##                    and compare the results with SPM.
##
##                    Basically search for differences between Control and Alzheimer 
##                    patients with SCCs and SPM. The true regions are previously known.
##                    So we compare SCCs and SPM against our true known regions and compute 
##                    metrics of accuracy.
##
## Notes: There are several time-consuming tasks which can be skipped by using the
##        load() function which appears commented at the end of the code chunk. Read
##        the code before running and you might save some time.
##
## Date Created: 2022-04-19
## Last Update: 2024-06-20
##
## Author: Juan A. Arias (M.Sc.)
## Email: juanantonio.arias.lopez@usc.es
## Webpage: https://juan-arias.xyz
##   
############################## ################### ############## ### 

####
# 1) PREAMBLE  ---- 
####

#* Working directory ----
setwd("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM")

#* Options ----
options(scipen = 6, digits = 4) # View outputs in non-scientific notation
memory.limit(30000000)     # This is needed on some PCs to increase memory allowance

#* Install packgs ----

packgs <- c("gamair", "oro.nifti", "memisc", "devtools", "remotes", "readr", 
            "imager", "itsadug", "fields", "BPST", "triangulation", "ImageSCC", 
            "tidyr", "dplyr", "stringr", "threadr", "memisc", "ggplot2")

for(packgs in packgs) {
  if (!requireNamespace(packgs, quietly = TRUE)) {
    install.packages(packgs)
  } else {
    message(paste("Package", packgs, "is already installed"))
  }
}

remotes::install_github("skgrange/threadr")

# Then load them:
lapply(packgs, library, character.only = TRUE); library(threadr)

#* Package neuroSCC ----
remotes::install_github("iguanamarina/neuroSCC")
library(neuroSCC)

#* Other parameters ----

# Brain slice to analyze
param.z = 35

####  
# 2) CONTOURS OF NEURO-DATA ----------
####

# First thing we need is the contours of our data so that we can perform 
# Delaunay triangulations which are essential for the calculation of SCCs. 
# This is one of several steps in which appropiate pre-processing of the 
# neuroimage files is essential.

#* Basic Template ----

# Load the template
template = neuroSCC::neuroCleaner("Auxiliary Files/new_mask") 

# Keep the relevant slice
template <- subset(template, template$z == param.z)

# Get limits of the file structure 
x <- max(template$x) 
y <- max(template$y)  
xy <- x*y

# Change structure and keep just PET data 
template <- t(as.matrix(template[ , "pet"]))

# If there are NAs, replace them with zeros
template[is.nan(template)] <- 0

#* Assign Coordinates ---- 
x <- rep(1:x, each = y, length.out = xy)
y <- rep(1:y, length.out = xy)
z <- cbind(as.matrix(x), as.matrix(y))
dat <- as.data.frame(cbind(z, t(template)))
dat[is.na(dat)] <- 0
rownames(dat) <- NULL; rm(x, y)

#* Get neuroContour ----
# Get contour for the area where values change from 0 to 1 
library(contoureR); library(ggplot2)
contour = contoureR::getContourLines(dat, levels = c(0)) 

# Visually test that we got the contours
ggplot(contour, aes(x, y)) + geom_path() 


####
# 3) CREATE SCC MATRIXES FOR CONTROL GROUP ------
####

# We are going to compare two groups of images: Control vs Simulated Alzheimer
# But first we need to do some in-between-steps so that the data is in the correct 
# format for Functional Data Analysis (SCC is a FDA technique).

#* Load CN data and Create DB ----
setwd("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/PETimg_masked for simulations")

# All the control files follow this name pattern (w00) and only a numeric value changes
pattern <- "^masked_swwwC\\d+_tripleNormEsp_w00_rrec_OSEM3D_32_it1.nii"

# Get the list of files with that name pattern
files <- list.files(pattern = pattern, full.names = FALSE)

# Create the database for CN files
database_CN <- data.frame(CN_number = integer(), z = integer(), x = integer(), y = integer(), pet = integer())

# Loop to process each file
for(i in 1:length(files)) {
  # Use neuroSCC::neuroCleaner to process the file
  temporal <- neuroSCC::neuroCleaner(files[i])
  # Extract the number from the file name using a regular expression
  CN_number <- sub("masked_swwwC(\\d+)_.*", "\\1", basename(files[i]))
  # Print the current control number being processed
  print(paste("Processing Control Nº", CN_number))
  # Create a column with the extracted number
  CN_number <- rep(CN_number, length.out = nrow(temporal))
  # Add the column to the temporary dataframe
  temporal <- cbind(CN_number, temporal)
  # Append the processed data to the main dataframe
  database_CN <- rbind(database_CN, temporal)
}

#* Create SCC Matrix ----

# Working on a functional data setup requires for the data to be in a concrete format which
# usually implies a long line of data points so that each row represents a function

# Preallocate the matrix
SCC_CN <- matrix(nrow = length(files), ncol = xy)

# Loop through the files
for(i in seq_along(files)) {
  
  # Get the CN_number corresponding to the current file
  CN_number <- sub("masked_swwwC(\\d+)_.*", "\\1", basename(files[i]))
  print(paste("Converting Control Nº", CN_number))
  # Subset the dataframe according to parameters
  subset_data <- database_CN[database_CN$CN_number == CN_number & database_CN$z == param.z, ]
  Y <- subset_data[1:xy, "pet"]
  # Convert to matrix and transpose
  Y <- t(as.matrix(Y))
  # Replace NaN values with 0
  Y[is.nan(Y)] <- 0
  # Assign the values to the matrix
  SCC_CN[i, ] <- Y
  
}

# Now it should be in matrix format with every row representing a Control file 
# setwd(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", as.character(param.z)))
# save(SCC_CN, file = "SCC_CN.RData") # SCC matrix for Controls

# Load results to save time:
# load("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z35/SCC_CN.RData")

