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
##                    So we compare SCCs and SPM against or true known regions and compute 
##                    metrics if accuracy.
##
## Date Created: 2024-06-19
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
            "imager", "itsadug", "fields", "BPST", "Triangulation", "ImageSCC", 
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
rownames(dat) <- NULL; rm(x, y, dat)

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

#* Load data and Create DB ----
setwd("~/GitHub/SCCneuroimage/PETimg_masked for simulations")

number <- paste0("C", 1:25)
name <- paste0("masked_swww", number, "_tripleNormEsp_w00_rrec_OSEM3D_32_it1")
# Only 25 files are controls and they follow the above defined structure

database_CN <- data.frame(CN_number = integer(),z = integer(), x = integer(), y = integer(), pet = integer())

for (i in 1:length(name)) {
  
  temporal <- f.clean(name[i])
  CN_number <- rep(number[i], length.out = nrow(Z))
  temporal <- cbind(CN_number,temporal)
  database_CN <- rbind(database_CN,temporal)
}

nrow(database_CN[database_CN$pet < 0, ]) # No negative values
rm(temporal); rm(CN_number)


#* Create SCC Matrix ----

# Working on a functional data setup requires for the data to be in a concrete format which
# usually implies a long line of data points so that each row represents a function

SCC_CN <- matrix(nrow = length(name), ncol = nrow(Z))

for (i in 1:length(number)) {
  
  Y <- subset(database_CN, database_CN$CN_number == number[i] & database_CN$z == param.z) 
  Y <- Y[1:9919, 5] 
  Y <- as.matrix(Y)
  Y = t(Y) 
  Y[is.nan(Y)] <- 0
  SCC_CN[i, ] <- Y
  
}

# Sometimes R really doesn't want to remove NA so this might be necessary:

# na.zero <- function(x) {
#     x[is.na(x)] <- 0
# }
# SCC_CN <-  apply(SCC_CN, 2, na.zero)


setwd(paste0("~/GitHub/SCCneuroimage/z", as.character(param.z)))
save(SCC_CN, file = "SCC_CN.RData") # SCC matrix for Controls


