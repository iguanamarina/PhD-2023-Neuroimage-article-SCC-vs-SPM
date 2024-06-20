############################## ################### ############## ### 
##
## Script name: 1. MASTERSCRIPT FOR GROUP COMPARISONS.r
##
## Purpose of script: Use SCCs for the estimation of regions with different
##                    activity levels in PET images and also to compare the
##                    results with SPM.
##
## Date Created: 2024-06-19
##
## Author: Juan A. Arias (M.Sc.)
## Email: juanantonio.arias.lopez@usc.es
## Webpage: https://juan-arias.xyz
##   
############################## ################### ############## ### 

# 1) PREAMBLE  ---- 

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

  
# 2) CONTOURS OF NEURO-DATA ----------

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

