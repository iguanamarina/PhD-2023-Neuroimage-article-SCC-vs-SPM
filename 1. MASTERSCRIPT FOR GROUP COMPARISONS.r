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
rownames(dat) <- NULL

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

# SCC_CN <- neuroSCC::matrixCreator(database_CN, pattern, param.z, xy)

# Now it should be in matrix format with every row representing a Control file 
# setwd(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", as.character(param.z)))
# save(SCC_CN, file = "SCC_CN.RData") # SCC matrix for Controls

# Load results to save time:
load("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z35/SCC_CN.RData")


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


####
# 7) SCC EVALUATION ----
####

# This section requires loading the TRUE points with differences in PET activity to test
# against them. These are called ROI_something. Basically, we will get the TRUE POINTS 
# according to the creators of these simulation files at Santiago de Compostela's Hospital
# and then test them against SCC and SPM and get some metrics on each method's efficiency.

#* Theoretical ROIs ----

# Define the base directory (can be skipped if you ran the whole script)
base_dir <- "~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM"

# Define the regions to be processed (same as previous)
regions <- c("w32", "w79", "w214", "w271", "w413", "wroiAD")

# Define the number of patients (max_number is previously calculated)
numbers <- 1:max_number  

# Call processROIs function from the neuroSCC package
neuroSCC::processROIs(base_dir, regions, numbers)

# Load ROI data
library(tidyverse)

# Define the table ROI directory 
roi_dir <- "~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/roisNormalizadas/tables"

# Load all the ROI tables
ROI_data <- lapply(seq_along(regions), function(i) {
  lapply(numbers, function(k) {
    readRDS(file.path(roi_dir, paste0("ROItable_", regions[i], "_", k, ".RDS")))
  })
})

# Filter and combine ROI data by z coordinate and pet value
T_points <- list() # Initialize an empty list to store T_points
for (i in seq_along(regions)) {
  for (j in seq_along(numbers)) {
    df <- ROI_data[[i]][[j]] %>%
      dplyr::filter(z == param.z & pet == 1) %>%
      dplyr::select(y, x) %>%
      tidyr::unite(newcol, y, x, remove = TRUE)
    T_points[[paste0(as.character(regions[i]), "_C", as.character(numbers[j]))]] <- df
  }
}

# Flatten the list of T_points
# This is the list of TRUE points where the simulation has been carried out
T_points <- unlist(T_points, recursive = FALSE)

# Remove unwanted part from the names of the list elements
names(T_points) <- sub("\\.newcol$", "", names(T_points))

# Clean up
rm(ROI_data)


#* Hypothetical points according to SCCs ---- 

# Hypothetical points and sens/esp for different methods go all together
# First, you need to perform the analysis with SPM and export results to a 
# binary file named "binary.nii"
# Pay attention to the folders so that you understand where things come from.

# Prelocate data.frames:

SCC_vs_SPM_complete <- data.frame(
  method = integer(),
  region = integer(),
  roi = integer(),
  number = integer(),
  sens = integer(),
  esp = integer(),
  PPV = integer(),
  NPV = integer()
)

SCC_vs_SPM <- data.frame(
  method = integer(),
  region = integer(),
  roi = integer(),
  sensMEAN = integer(),
  sensSD = integer(),
  espMEAN = integer(),
  espSD = integer(),
  ppvMEAN = integer(),
  ppvSD = integer(),
  npvMEAN = integer(),
  npvSD = integer()
)

# Load the template and automatically detect the limits of x and y
setwd("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM")
template <- neuroSCC::neuroCleaner("Auxiliary Files/new_mask")

# Keep the relevant slice
template <- subset(template, template$z == param.z)

# Get the limits of the file structure
x <- max(template$x)
y <- max(template$y)
xy <- x * y

# Calculate the combinations of coordinates present
x_coords <- rep(1:x, each = y, length.out = xy)
y_coords <- rep(1:y, length.out = xy)
total_coords <- data.frame(y = y_coords, x = x_coords)
total_coords <- tidyr::unite(total_coords, newcol, c(y, x), remove = TRUE)

# Clean up temporary variables
rm(x_coords); rm(y_coords)

# Main loop to process SCC and SPM results (fasten your seatbelts)
for (k in 1:length(roi)) {
  
  regions <- c("w32", "w79", "w214", "w271", "w413", "roiAD")
  
  H_points <- list()
  
  # Set working directory to results folder
  setwd(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", 
               as.character(param.z), "/results"))
  
  # Load SCC results for each region
  for (i in 1:length(regions)) {
    load(paste0("SCC_COMP_", regions[i], "_", roi[k], ".RData"))
    H_points[[as.character(regions[i])]] <- neuroSCC::getPoints(SCC_COMP)[[1]]
    H_points[[as.character(regions[i])]] <- tidyr::unite(
      as.data.frame(H_points[[i]]), newcol, c(row, col), remove = TRUE
    )
  }
  
  rm(list = ls(pattern = "^SCC_COMP"))
  
  # Ensure the correct name for wroiAD region
  regions <- c("w32", "w79", "w214", "w271", "w413", "wroiAD")
  names(H_points)[6] <- "wroiAD"
  
  # Calculate sensitivity and specificity for SCC
  for (i in 1:length(regions)) {
    SCC_sens_esp <- data.frame(region = integer(), group = integer(), 
                               sens = integer(), esp = integer(), 
                               ppv = integer(), npv = integer())
    
    for (j in 1:length(number)) {
      inters <- dplyr::inner_join(
        H_points[[as.character(regions[i])]], 
        data.frame(newcol = T_points[[paste0(as.character(regions[i]), "_", 
                                             as.character(number[j]))]]), 
        by = "newcol"
      )
      
      # Sensitivity = TP / (TP + FN)
      sensibilitySCC <- nrow(inters) / length(T_points[[paste0(as.character(regions[i]), "_", 
                                                               as.character(number[j]))]]) * 100
      
      true_neg <- dplyr::anti_join(
        total_coords,
        data.frame(newcol = T_points[[paste0(as.character(regions[i]), "_", 
                                             as.character(number[j]))]]),
        by = "newcol"
      )
      hypo_neg <- dplyr::anti_join(total_coords, H_points[[as.character(regions[i])]], by = "newcol")
      
      anti_inters <- dplyr::inner_join(true_neg, hypo_neg, by = "newcol")
      
      # Specificity = TN / (TN + FP)
      specificitySCC <- nrow(anti_inters) / nrow(true_neg) * 100
      
      # PPV = TP / (TP + FP) -> probability of having the disease after a positive test result
      FalsePositive <- dplyr::inner_join(
        H_points[[as.character(regions[i])]], true_neg, by = "newcol"
      )
      FalseNegative <- dplyr::inner_join(
        hypo_neg, data.frame(newcol = T_points[[paste0(as.character(regions[i]), "_", 
                                                       as.character(number[j]))]]), 
        by = "newcol"
      )
      
      PPV <- (nrow(inters) / (nrow(inters) + nrow(FalsePositive))) * 100
      
      # NPV = TN / (FN + TN) -> probability of not having the disease after a negative test result
      NPV <- (nrow(anti_inters) / (nrow(anti_inters) + nrow(FalseNegative))) * 100
      
      temp <- data.frame(region = regions[i], group = number[j], sens = sensibilitySCC, 
                         esp = specificitySCC, ppv = PPV, npv = NPV)
      SCC_sens_esp <- rbind(SCC_sens_esp, temp)
    }
    
    means <- data.frame(region = regions[i], group = "MEAN", sens = mean(SCC_sens_esp$sens, na.rm = TRUE), 
                        esp = mean(SCC_sens_esp$esp, na.rm = TRUE),
                        ppv = mean(SCC_sens_esp$ppv, na.rm = TRUE),
                        npv = mean(SCC_sens_esp$npv, na.rm = TRUE))
    sds <- data.frame(region = regions[i], group = "SD", sens = sd(SCC_sens_esp$sens, na.rm = TRUE), 
                      esp = sd(SCC_sens_esp$esp, na.rm = TRUE),
                      ppv = sd(SCC_sens_esp$ppv, na.rm = TRUE),
                      npv = sd(SCC_sens_esp$npv, na.rm = TRUE))
    SCC_sens_esp <- rbind(SCC_sens_esp, means, sds)
    
    setwd(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", 
                 as.character(param.z), "/results", "/ROI", roi[k]))
    readr::write_csv(SCC_sens_esp, paste0("sens_esp_SCC_", regions[i], "_", roi[k], ".csv"), 
                     na = "NA", append = FALSE)
  }
  
  setwd(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", 
               as.character(param.z), "/results", "/ROI", roi[k]))
  
  sens_esp_SCC_w32 <- readr::read_csv(paste0("sens_esp_SCC_w32_", roi[k], ".csv"), 
                                      na = "NA", show_col_types = FALSE)
  sens_esp_SCC_w79 <- readr::read_csv(paste0("sens_esp_SCC_w79_", roi[k], ".csv"), 
                                      na = "NA", show_col_types = FALSE)
  sens_esp_SCC_w214 <- readr::read_csv(paste0("sens_esp_SCC_w214_", roi[k], ".csv"), 
                                       na = "NA", show_col_types = FALSE)
  sens_esp_SCC_w271 <- readr::read_csv(paste0("sens_esp_SCC_w271_", roi[k], ".csv"), 
                                       na = "NA", show_col_types = FALSE)
  sens_esp_SCC_w413 <- readr::read_csv(paste0("sens_esp_SCC_w413_", roi[k], ".csv"), 
                                       na = "NA", show_col_types = FALSE)
  sens_esp_SCC_wroiAD <- readr::read_csv(paste0("sens_esp_SCC_wroiAD_", roi[k], ".csv"), 
                                         na = "NA", show_col_types = FALSE)
  
  # Sensitivity and Specificity for SPM
  for (i in 1:length(regions)) {
    setwd(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", 
                 as.numeric(param.z), "/SPM", "/ROI", i, "_", regions[i], "_0", roi[k]))
    binary <- neuroSCC::neuroCleaner("binary.nii")
    H_points_SPM <- binary[binary$z == as.numeric(param.z) & binary$pet == 1, 2:3]
    H_points_SPM <- tidyr::unite(as.data.frame(H_points_SPM), newcol, c(y, x), remove = TRUE)
    
    SPM_sens_esp <- data.frame(region = integer(), group = integer(), sens = integer(), 
                               esp = integer(), ppv = integer(), npv = integer())
    
    for (j in 1:length(number)) {
      inters <- dplyr::inner_join(
        H_points_SPM, 
        data.frame(newcol = T_points[[paste0(as.character(regions[i]), "_", 
                                             as.character(number[j]))]]), 
        by = "newcol"
      )
      
      # Sensitivity = TP / (TP + FN)
      sensibilitySPM <- nrow(inters) / length(T_points[[paste0(as.character(regions[i]), "_", 
                                                               as.character(number[j]))]]) * 100
      
      true_neg <- dplyr::anti_join(
        total_coords, 
        data.frame(newcol = T_points[[paste0(as.character(regions[i]), "_", 
                                             as.character(number[j]))]]), 
        by = "newcol"
      )
      
      hypo_neg <- dplyr::anti_join(total_coords, H_points_SPM, by = "newcol")
      
      anti_inters <- dplyr::inner_join(true_neg, hypo_neg, by = "newcol")
      
      # Specificity = TN / (TN + FP)
      specificitySPM <- nrow(anti_inters) / nrow(true_neg) * 100
      
      # PPV = TP/(TP+FP) -> probability of having the disease after a positive test result
      FalsePositive <- dplyr::inner_join(H_points_SPM, true_neg, by = "newcol")
      FalseNegative <- dplyr::inner_join(
        hypo_neg, 
        data.frame(newcol = T_points[[paste0(as.character(regions[i]), "_", as.character(number[j]))]]), 
        by = "newcol"
      )
      
      PPV <- (nrow(inters) / (nrow(inters) + nrow(FalsePositive))) * 100
      
      # NPV = TN/(FN+TN) -> probability of not having the disease after a negative test result
      NPV <- (nrow(anti_inters) / (nrow(anti_inters) + nrow(FalseNegative))) * 100
      
      temp <- data.frame(region = regions[i], group = number[j], sens = sensibilitySPM, 
                         esp = specificitySPM, ppv = PPV, npv = NPV)
      SPM_sens_esp <- rbind(SPM_sens_esp, temp)
    }
    
    means <- data.frame(region = regions[i], group = "MEAN", sens = mean(SPM_sens_esp$sens, na.rm = TRUE), 
                        esp = mean(SPM_sens_esp$esp, na.rm = TRUE),
                        ppv = mean(SPM_sens_esp$ppv, na.rm = TRUE),
                        npv = mean(SPM_sens_esp$npv, na.rm = TRUE))
    sds <- data.frame(region = regions[i], group = "SD", sens = sd(SPM_sens_esp$sens, na.rm = TRUE), 
                      esp = sd(SPM_sens_esp$esp, na.rm = TRUE),
                      ppv = sd(SPM_sens_esp$ppv, na.rm = TRUE),
                      npv = sd(SPM_sens_esp$npv, na.rm = TRUE))
    SPM_sens_esp <- rbind(SPM_sens_esp, means, sds)
    
    setwd(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", 
                 as.character(param.z), "/results", "/ROI", roi[k]))
    
    readr::write_csv(SPM_sens_esp, paste0("sens_esp_SPM_", regions[i], "_", roi[k], ".csv"), 
                     na = "NA", append = FALSE)
  }
  
  setwd(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", 
               as.character(param.z), "/results", "/ROI", roi[k]))
  
  sens_esp_SPM_w32 <- readr::read_csv(paste0("sens_esp_SPM_w32_", roi[k], ".csv"), 
                                      na = "NA", show_col_types = FALSE)
  sens_esp_SPM_w79 <- readr::read_csv(paste0("sens_esp_SPM_w79_", roi[k], ".csv"), 
                                      na = "NA", show_col_types = FALSE)
  sens_esp_SPM_w214 <- readr::read_csv(paste0("sens_esp_SPM_w214_", roi[k], ".csv"), 
                                       na = "NA", show_col_types = FALSE)
  sens_esp_SPM_w271 <- readr::read_csv(paste0("sens_esp_SPM_w271_", roi[k], ".csv"), 
                                       na = "NA", show_col_types = FALSE)
  sens_esp_SPM_w413 <- readr::read_csv(paste0("sens_esp_SPM_w413_", roi[k], ".csv"), 
                                       na = "NA", show_col_types = FALSE)
  sens_esp_SPM_wroiAD <- readr::read_csv(paste0("sens_esp_SPM_wroiAD_", roi[k], ".csv"), 
                                         na = "NA", show_col_types = FALSE)
  
  #* Create Complete Lists: ----
  listSCC <- ls(pattern = "^sens_esp_SCC")
  
  for (i in 1:length(listSCC)) {
    data <- get(listSCC[[i]])[1:max_number, ]
    method <- rep("SCC", times = max_number)
    Roi <- rep(as.character(roi[k]), times = max_number) 
    tempSCC <- cbind(as.data.frame(method), data[, 1], as.data.frame(Roi), data[, 2:6])  
    SCC_vs_SPM_complete <- rbind(SCC_vs_SPM_complete, tempSCC)
  }
  
  listSPM <- ls(pattern = "^sens_esp_SPM")
  
  for (i in 1:length(listSPM)) {
    data <- get(listSPM[[i]])[1:max_number, ]
    method <- rep("SPM", times = max_number)
    Roi <- rep(as.character(roi[k]), times = max_number)  
    tempSPM <- cbind(as.data.frame(method), data[, 1], as.data.frame(Roi), data[, 2:6])  
    SCC_vs_SPM_complete <- rbind(SCC_vs_SPM_complete, tempSPM)
  }
  
  #* Create a Final Reduced List: ----
  for (i in 1:length(listSCC)) {
    data <- get(listSCC[[i]])[(max_number + 1):(max_number + 2), ]
    temp <- data.frame(method = "SCC",
                       region = as.character(data$region[1]), 
                       roi = roi[k], 
                       sensMEAN = data$sens[1], 
                       sensSD = data$sens[2], 
                       espMEAN = data$esp[1], 
                       espSD = data$esp[2],
                       ppvMEAN = data$ppv[1],
                       ppvSD = data$ppv[2],
                       npvMEAN = data$npv[1],
                       npvSD = data$npv[2])    
    SCC_vs_SPM <- rbind(SCC_vs_SPM, temp)
  }
  
  for (i in 1:length(listSPM)) {
    data <- get(listSPM[[i]])[(max_number + 1):(max_number + 2), ]
    temp <- data.frame(method = "SPM",
                       region = as.character(data$region[1]), 
                       roi = roi[k], 
                       sensMEAN = data$sens[1], 
                       sensSD = data$sens[2], 
                       espMEAN = data$esp[1], 
                       espSD = data$esp[2],
                       ppvMEAN = data$ppv[1],
                       ppvSD = data$ppv[2],
                       npvMEAN = data$npv[1],
                       npvSD = data$npv[2]) 
    SCC_vs_SPM <- rbind(SCC_vs_SPM, temp)
  }
}

# View(SCC_vs_SPM)
# View(SCC_vs_SPM_complete)

# setwd(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", as.numeric(param.z), "/results"))
# saveRDS(SCC_vs_SPM, file = "SCC_vs_SPM.RDS")
# saveRDS(SCC_vs_SPM_complete, file = "SCC_vs_SPM_complete.RDS")

#* Export as LaTeX table code ----

setwd("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Article")

library(dplyr)

SCC_vs_SPM <- readRDS(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", 
                             as.numeric(param.z), "/results/SCC_vs_SPM.RDS"))

SCC_vs_SPM <- SCC_vs_SPM %>% dplyr::arrange(method)

SCC_vs_SPM$region <-   factor(SCC_vs_SPM$region,      
                              levels = c("w32", "w79", "w214", "w271", "w413", "wroiAD"))

SCC_vs_SPM$roi <- as.numeric(SCC_vs_SPM$roi)
SCC_vs_SPM$roi <- as.character(as.numeric(SCC_vs_SPM$roi) * 10)

latex <- SCC_vs_SPM %>%
  mutate(Sensibility = paste(round(sensMEAN, 2), round(sensSD, 2), sep = "*")) %>%
  mutate(Specificity = paste(round(espMEAN, 2), round(espSD, 2), sep = "*")) %>%
  mutate(PPV = paste(round(ppvMEAN, 2), round(ppvSD, 2), sep = "*")) %>%
  mutate(NPV = paste(round(npvMEAN, 2), round(npvSD, 2), sep = "*")) %>%
  arrange(method, roi, region) %>%
  dplyr::select(method, roi, region, Sensibility, Specificity, PPV, NPV)

latex

# install.packages("xtable")
library(xtable)
print(xtable(latex, type = "latex"), 
      file = "table.tex",
      include.rownames = FALSE)


####  
# 8) VISUALIZATIONS----------
####  

# Load data if not already in environment
referencia <- readRDS(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", 
                             as.numeric(param.z), "/results/SCC_vs_SPM.RDS"))
table <- readRDS(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", 
                        as.numeric(param.z), "/results/SCC_vs_SPM_complete.RDS"))

# Load necessary packages
library(tidyverse)
library(lemon)
library(ggplot2)
library(gridExtra)
# library(envalysis) # This package allows for the publish_theme() but it 
                     # might not work on your version, that's why we stick to minimal_theme()

# Modify data structure
table$region <- factor(table$region,
                       levels = c("w32", "w79", "w214", "w271", "w413", "wroiAD"))
table$Roi <- as.character(as.numeric(table$Roi) * 10)

# Merge tables
table2 <- rbind(table[table$method == "SCC", ], table[table$method == "SPM", ])

# Set working directory for saving figures
setwd(paste0("~/Documents/GitHub/PhD-2023-Neuroimage-article-SCC-vs-SPM/Results/z", 
             as.numeric(param.z), "/Figures")) 

#* Sensitivity and Specificity for wroiAD ----
graph1 <- ggplot(data = table2[table2$region == "wroiAD", ], aes(x = Roi, y = sens)) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_boxplot(aes(fill = method)) +
  xlab("Level of Induced Hypoactivity (%)") +
  ylab("Sensitivity (%)") +
  guides(fill = guide_legend(title = "Legend")) +
  scale_fill_brewer(palette = "Set1")

graph2 <- ggplot(data = table2[table2$region == "wroiAD", ], aes(x = Roi, y = esp)) +
  geom_boxplot(aes(fill = method, col = method)) +
  guides(col = FALSE) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  xlab("Level of Induced Hypoactivity (%)") +
  ylab("Specificity (%)") +
  guides(fill = guide_legend(title = "Legend")) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1")

# Arrange graphs side by side
combined_graph <- grid.arrange(graph1, graph2, ncol = 2)

# Save combined graph as PNG
ggsave(filename = paste0("sens_esp_wroiAD_", as.numeric(param.z), ".png"), 
       plot = combined_graph, 
       width = 24, 
       height = 18, 
       units = "cm",
       dpi = 600)

#* Sensitivity for all regions and ROIs ----
graph3 <- ggplot(data = table2, aes(x = Roi, y = sens)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_boxplot(outlier.colour = NULL, aes(fill = method), outlier.size = 1, lwd = 0.25) +
  xlab("Hypoactivity (%)") +
  ylab("Sensitivity (%)") +
  guides(fill = guide_legend(title = "Legend")) +
  facet_wrap(~region, ncol = 2) +
  facet_rep_wrap(~region, repeat.tick.labels = TRUE) +
  coord_capped_cart(bottom = "both", left = "both") +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(legend.text = element_text(size = 14)) +
  scale_fill_brewer(palette = "Set1")

# Arrange graph
plot3 <- grid.arrange(graph3)

# Save sensitivity graph for all regions and ROIs
ggsave(filename = paste0("sens_ALL_", as.numeric(param.z), ".png"), 
       plot = plot3, 
       width = 28.95, 
       height = 18.3, 
       units = "cm",
       dpi = 600)

graph4 <- ggplot(data = table2, aes(x = Roi, y = esp)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_boxplot(outlier.colour = NULL, aes(fill = method, col = method), outlier.size = 1, lwd = 0.25) +
  guides(col = FALSE) +
  xlab("Hypoactivity (%)") +
  ylab("Specificity (%)") +
  coord_cartesian(ylim = c(0, 100)) +
  guides(fill = guide_legend(title = "Legend")) +
  facet_wrap(~region, nrow = 3) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1")

# Arrange graph
plot4 <- grid.arrange(graph4)

# Save specificity graph for all regions and ROIs
ggsave(filename = paste0("esp_ALL_", as.numeric(param.z), ".png"), 
       plot = plot4, 
       width = 22, 
       height = 26, 
       units = "cm",
       dpi = 600)

# Arrange sensitivity and specificity graphs side by side
plot_together <- grid.arrange(plot3, plot4, ncol = 2)

# Save combined sensitivity and specificity graph
ggsave(filename = paste0("sens_esp_ALL_", as.numeric(param.z), ".png"), 
       plot = plot_together, 
       width = 22, 
       height = 26, 
       units = "cm",
       dpi = 600)

#* PPV and NPV for wroiAD ----
graph5 <- ggplot(data = table2, aes(x = Roi, y = ppv)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_boxplot(aes(fill = method)) +
  xlab("Hypoactivity (%)") +
  ylab("Positive Predictive Value (%)") +
  guides(fill = guide_legend(title = "Legend")) +
  facet_wrap(~region, ncol = 2) +
  facet_rep_wrap(~region, repeat.tick.labels = TRUE) +
  coord_capped_cart(bottom = "both", left = "both") +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(0, 100))

graph6 <- ggplot(data = table2, aes(x = Roi, y = npv)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  coord_cartesian(ylim = c(0, 100)) +
  geom_boxplot(aes(fill = method)) + 
  xlab("Level of Induced Hypoactivity (%)") + 
  ylab("Negative Predictive Value (%)") + 
  guides(fill = guide_legend(title = "Legend")) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1")

# Arrange graph
plot_ppv <- grid.arrange(graph5)
plot_npv <- grid.arrange(graph6)

# Save PPV graph for wroiAD
ggsave(filename = paste0("PPV_wroiAD_", as.numeric(param.z), ".png"), 
       plot = plot_ppv, 
       width = 22, 
       height = 26, 
       units = "cm",
       dpi = 600)

# Save NPV graph for wroiAD
ggsave(filename = paste0("NPV_wroiAD_", as.numeric(param.z), ".png"), 
       plot = plot_npv, 
       width = 22, 
       height = 26, 
       units = "cm",
       dpi = 600)

#* PPV and NPV for all regions and ROIs ----
graph7 <- ggplot(data = table2, aes(x = Roi, y = ppv)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_boxplot(outlier.colour = NULL, aes(fill = method), outlier.size = 1, lwd = 0.25) +
  xlab("Hypoactivity (%)") +
  ylab("Positive Predictive Value (%)") +
  guides(fill = guide_legend(title = "Legend")) +
  facet_wrap(~region, ncol = 3) +
  facet_rep_wrap(~region, repeat.tick.labels = TRUE) +
  coord_capped_cart(bottom = "both", left = "both") +
  theme(panel.border = element_blank(), axis.line = element_line()) +
  theme(legend.text = element_text(size = 14)) +
  scale_fill_brewer(palette = "Set1")

# Arrange graph
plot7 <- grid.arrange(graph7)

# Save PPV graph for all regions and ROIs
ggsave(filename = paste0("ppv_ALL_", as.numeric(param.z), ".png"), 
       plot = plot7, 
       width = 28.95, 
       height = 18.3, 
       units = "cm",
       dpi = 600)

graph8 <- ggplot(data = table2, aes(x = Roi, y = npv)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_boxplot(outlier.colour = NULL, aes(fill = method, col = method), outlier.size = 1, lwd = 0.25) +
  xlab("Hypoactivity (%)") +
  ylab("Negative Predictive Value (%)") +
  coord_cartesian(ylim = c(0, 100)) +
  guides(fill = guide_legend(title = "Legend")) +
  guides(col = FALSE) +
  facet_wrap(~region, nrow = 3) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1")

# Arrange graphs
plot8 <- grid.arrange(graph8)
plot9 <- grid.arrange(plot7, plot8, ncol = 2)

# Save NPV graph for all regions and ROIs
ggsave(filename = paste0("npv_ALL_", as.numeric(param.z), ".png"), 
       plot = plot8, 
       width = 22, 
       height = 26, 
       units = "cm",
       dpi = 600)

# Save combined PPV and NPV graph for all regions and ROIs
ggsave(filename = paste0("ppv_npv_ALL_", as.numeric(param.z), ".png"), 
       plot = plot9, 
       width = 22, 
       height = 26, 
       units = "cm",
       dpi = 600)


