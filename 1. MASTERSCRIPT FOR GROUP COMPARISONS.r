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

#* Install CRAN packgs ----

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

# Set environment to not stop installations due to warning messages
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE)

#* Install GitHub packages ----
packagesGithub <- c("funstatpackages/BPST", "funstatpackages/Triangulation", "funstatpackages/ImageSCC")

# Install and load GitHub packages
for (pkg in packagesGithub) {
  remotes::install_github(pkg)
}

library(BPST); library(Triangulation); library(ImageSCC)

#* Install package neuroSCC ----
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

# Get contour for the area where values change from 0 to 1 using neuroSCC::neuroContour()
contourCoordinates <- neuroSCC::neuroContour(dat, levels = c(0))

# Test the results
plot(contourCoordinates[[1]])       # External boundaries
if (length(contourCoordinates) > 1) {
  for (j in 2:length(contourCoordinates)) {
    points(contourCoordinates[[j]]) # Holes or internal contours if present
  }
}

#* Get coordinates in pckg Triangulation format: ----

VT = Triangulation::TriMesh(contourCoordinates[[1]], n = 8) 

# n = Triangulation degree of fineness (8 is recommended)
# However, higher values can be used as Arias-López et al. (2021) suggests 
# computing times for higher n values are still sensible.

#* Save these contour coordinates ----
# Define the directory path
directoryPath <- paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", 
                        as.character(param.z))

# Check if the directory exists and create it if it does not
if (!dir.exists(directoryPath)) {
  message("This folder does not exist and will be created: ", directoryPath)
  dir.create(directoryPath, recursive = TRUE)
} else {
  message("This folder already exists: ", directoryPath)
}

# Set the working directory and save
setwd(directoryPath)
save(VT, file = paste0("contour", as.character(param.z), ".RData"))


####
# 3) CREATE SCC MATRIXES FOR CONTROL GROUP ------
####

# We are going to compare two groups of images: Control vs Simulated Alzheimer
# But first we need to do some in-between-steps so that the data is in the correct 
# format for Functional Data Analysis (SCC is a FDA technique).

#* Create CN DataBase ----

# Set the working directory for image files
setwd("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/PETimg_masked for simulations")

# Create the database using a specific file pattern 
pattern <- "^masked_swwwC\\d+_tripleNormEsp_w00_rrec_OSEM3D_32_it1.nii"

# Use pattern as parameter for neuroSCC::databaseCreator
database_CN <- neuroSCC::databaseCreator(pattern)


#* Create CN Matrix ----

# Assuming that 'database' is for Controls and that 'pattern', 'param.z', and 'xy' are defined in the script
SCC_CN <- neuroSCC::matrixCreator(database_CN, pattern, param.z, xy)

# Now it should be in matrix format with every row representing a Control file 
# setwd(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", as.character(param.z)))
# save(SCC_CN, file = "SCC_CN.RData") # SCC matrix for Controls

# Load results to save time:
# load("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z35/SCC_CN.RData")


####
# 4) CREATE SCC MATRIXES FOR PATHOLOGICAL GROUP ------
####

# Set initial working directory
base_dir <- "~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/PETimg_masked for simulations"
setwd(base_dir)

# Get the list of files matching the general pattern and extract maximum value
files <- list.files(pattern = "^masked_swwwC\\d+_tripleNormEsp_.*\\.nii$", full.names = TRUE)
max_number <- max(as.integer(sub("masked_swwwC(\\d+)_.*", "\\1", basename(files))), na.rm = TRUE)

#* Define parameters ----

# First the number of original patient
number <- paste0("C", 1:max_number)
# The selected region
region <- c("roiAD", "w32", "w79", "w214", "w271", "w413")
# And then the ROI
roi <- c(1, 2, 4, 6, 8)

#* Create AD DataBase and Matrix ----

for (i in 1:length(region)) {
  
  for (j in 1:length(roi)) {
    
    setwd(base_dir)
    
    # Tuning parameters
    REGION = region[i]
    ROI = roi[j]
    pattern <- paste0("^masked_swwwC\\d+_tripleNormEsp_", REGION, "_0_", ROI, "_rrec_OSEM3D_32_it1\\.nii")
    
    # Create the database for pathological data
    database_AD <- neuroSCC::databaseCreator(pattern, control = FALSE)
    
    # Create SCC Matrix
    SCC_matrix <- neuroSCC::matrixCreator(database_AD, pattern, param.z, xy)
    
    # Define result directory
    result_dir <- paste0(dirname(base_dir), "/Results/z", as.character(param.z))
    if (!dir.exists(result_dir)) {
      dir.create(result_dir, recursive = TRUE)
    }
    
    setwd(result_dir)
    
    file_path <- paste0("SCC_", REGION, "_", ROI, ".RData")
    
    # Option 1: Check if the file exists before saving
    if (!file.exists(file_path)) {
      save(SCC_matrix, file = file_path)
      message("File saved: ", file_path)
    } else {
      message("File already exists: ", file_path)
    }
    
    # Option 2: Always overwrite the file (commented out)
    # save(SCC_matrix, file = file_path)
    # message("File saved (overwritten): ", file_path)
  }
}

# In case we wanted to check one of these matrices you can custom use this code.
# Otherwise, there is no need to check them now as they will be called from 
# a loop later on:

# # Set REGION and ROI variables you want to check
# REGION <- "w32" # Options: "roiAD", "w32", "w79", "w214", "w271", "w413"
# ROI <- 1        # Options: 1, 2, 4, 6, 8
# 
# # Define the result directory
# result_dir <- paste0(dirname(base_dir), "/Results/z", as.character(param.z))
# 
# # Define the file path based on REGION and ROI
# file_path <- paste0(result_dir, "/SCC_", REGION, "_", ROI, ".RData")
# 
# # Load the SCC matrix as SCC_matrix
# load(file_path)


####
# 5) TRIANGULATION AND OTHER PARAMETERS ------
####

# Based on the results from section "2 - Get coordinates in pckg Triangulation format"
setwd(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", as.character(param.z)))

# Load previously calculated contour as VT 
load(paste0("contour", as.character(param.z), ".RData"))

# In order to be consistent with other packages and to avoid nomenclature problems
# we use common names Brain.V and Brain.Tr. 
# From here onwards most of the names follow the ones provided by Wang et al (2019)

Brain.V <- VT[[1]]
Brain.Tr <- VT[[2]]

V.est = as.matrix(Brain.V)
Tr.est = as.matrix(Brain.Tr)
V.band = as.matrix(Brain.V)
Tr.band = as.matrix(Brain.Tr) 


####
# 6) SCC ESTIMATIONS ------
####

# Set initial working directory
base_dir <- "~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM"
param_z_dir <- paste0(base_dir, "/Results/z", as.character(param.z))
result_dir <- paste0(param_z_dir, "/results")

# Define regions and ROIs
region <- c("roiAD", "w32", "w79", "w214", "w271", "w413")
roi <- c(1, 2, 4, 6, 8)

#* 'The Loop': Fasten your seat belts ----
for (i in 1:length(region)) {
  for (j in 1:length(roi)) {
    setwd(param_z_dir)
    
    # Define file names for CN and AD SCC data
    name_CN <- paste0(param_z_dir, "/SCC_CN.RData")
    name_AD <- paste0(param_z_dir, "/SCC_", region[i], "_", roi[j], ".RData")
    
    # Load SCC data
    SCC_CN <- threadr::read_rdata(name_CN)
    SCC_AD <- threadr::read_rdata(name_AD)
    
    #** Mean Average Normalization ----
    SCC_CN <- neuroSCC::meanNormalization(SCC_CN)
    SCC_AD <- neuroSCC::meanNormalization(SCC_AD)
    
    #** Parameters for SCC computation ----
    d.est <- 5 # degree of spline for mean function
    d.band <- 2 # degree of spline for SCC
    r <- 1 # smoothing parameter
    lambda <- 10^{seq(-6, 3, 0.5)} # penalty parameters
    alpha.grid <- c(0.10, 0.05, 0.01) # vector of confidence levels
    
    #** Construction of SCCs ----
    setwd(result_dir)
    
    result_file <- paste0("SCC_COMP_", region[i], "_", roi[j], ".RData")
    
    if (file.exists(result_file)) {
      print(paste0("Nice! The file ", as.character(result_file), " already exists.")) # Skip computation if file exists
    } else {
      SCC_COMP <- ImageSCC::scc.image(Ya = SCC_AD, Yb = SCC_CN, Z = z, 
                                      d.est = d.est, d.band = d.band, r = r,
                                      V.est.a = V.est, Tr.est.a = Tr.est,
                                      V.band.a = V.band, Tr.band.a = Tr.band,
                                      penalty = TRUE, lambda = lambda, alpha.grid = alpha.grid,
                                      adjust.sigma = TRUE)
      save(SCC_COMP, file = result_file)
    }
  }
}
