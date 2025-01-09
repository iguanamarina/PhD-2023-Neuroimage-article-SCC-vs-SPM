
### ########################################################## ###
#####             *SCRIPT YA PARA SCCs     *                  ####
### ########################################################## ###

library(gamair);library(oro.nifti);library(memisc);library(devtools);library(remotes);library(readr)
library(imager);library(itsadug);library(ggplot2);library(contoureR);library(fields);library(BPST)
library(Triangulation);library(ImageSCC)
load("contour35.RData")

# In order to be consistent we use common names Brain.V and Brain.Tr. From here onwards most of the names follow the ones provided by Wang et al (2019)

Brain.V <- VT[[1]]
Brain.Tr <- VT[[2]]

head(Brain.V);head(Brain.Tr)


V.est = as.matrix(Brain.V)
# Brain.v<-cbind(Brain.V[,2],Brain.V[,1]) # In case you need to transpose the data
Tr.est = as.matrix(Brain.Tr)

V.band = as.matrix(Brain.V)
Tr.band = as.matrix(Brain.Tr) 


# Response Variable:

load("SCC_CN.RData")
load("SCC_roiAD.RData")
load("SCC_w32.RData")
load("SCC_w79.RData")
load("SCC_w214.RData")
load("SCC_w271.RData")
load("SCC_w413.RData")

##########################################################
###           MEAN AVERAGE NORMALIZATION      (!!!!!)  ###  This can be reduced with functions (!!)
##########################################################

for (i in 1:nrow(SCC_CN)) {
  
  temp <- SCC_CN[i, ]
  mean <- mean(as.numeric(temp), na.rm = T)
  SCC_CN[i, ] <- (temp/mean)
}

for (i in 1:nrow(SCC_roiAD)) {
  
  temp <- SCC_roiAD[i, ]
  mean <- mean(as.numeric(temp), na.rm = T)
  SCC_roiAD[i, ] <- (temp/mean)
}

for (i in 1:nrow(SCC_w32)) {
  
  temp <- SCC_w32[i, ]
  mean <- mean(as.numeric(temp), na.rm = T)
  SCC_w32[i, ] <- (temp/mean)
}

for (i in 1:nrow(SCC_w79)) {
  
  temp <- SCC_w79[i, ]
  mean <- mean(as.numeric(temp), na.rm = T)
  SCC_w79[i, ] <- (temp/mean)
}

for (i in 1:nrow(SCC_w214)) {
  
  temp <- SCC_w214[i, ]
  mean <- mean(as.numeric(temp), na.rm = T)
  SCC_w214[i, ] <- (temp/mean)
}

for (i in 1:nrow(SCC_w271)) {
  
  temp <- SCC_w271[i, ]
  mean <- mean(as.numeric(temp), na.rm = T)
  SCC_w271[i, ] <- (temp/mean)
}

for (i in 1:nrow(SCC_w413)) {
  
  temp <- SCC_w413[i, ]
  mean <- mean(as.numeric(temp), na.rm = T)
  SCC_w413[i, ] <- (temp/mean)
}



### ########################################################## ###
#####        *OTHER PARAMETERS FOR SCC ESTIMATION*            ####
### ########################################################## ###

x <- rep(1:91, each = 109, length.out = 9919) 
y <- rep(1:109,length.out = 9919)
Z <- cbind(as.matrix(x),as.matrix(y)); Z

# Following Wang et al's recomendations:

d.est = 5 # degree of spline for mean function  5
d.band = 2 # degree of spline for SCC  2
r = 1 # smoothing parameter  1
lambda = 10^{seq(-6,3,0.5)} # penalty parameters
alpha.grid = c(0.10,0.05,0.01) # vector of confidence levels


### ########################################################## ###
#####               *CONSTRUCTION OF SCC'S*                   ####
### ########################################################## ###

setwd("~/Escritorio/SCCs_simulaciones_GAS/Results")

save(SCC_CN, file="SCC_CN.Rdata")
save(SCC_roiAD, file="SCC_roiAD.Rdata")
save(SCC_w32, file="SCC_w32.Rdata")
save(SCC_w79, file="SCC_w79.Rdata")
save(SCC_w214, file="SCC_w214.Rdata")
save(SCC_w271, file="SCC_w271.Rdata")
save(SCC_w413, file="SCC_w413.Rdata")

## Asignación previa:

Y_CN = SCC_CN
Y_AD = SCC_roiAD
Y_AD_32 = SCC_w32
Y_AD_79 = SCC_w79
Y_AD_214 = SCC_w214
Y_AD_271 = SCC_w271
Y_AD_413 = SCC_w413

################

setwd("~/Escritorio/SCCs_simulaciones_GAS/Results")

SCC_COMP_1 = scc.image(Ya = Y_AD, Yb = Y_CN, Z = Z, d.est = d.est, d.band = d.band, r = r,
                       V.est.a = V.est, Tr.est.a = Tr.est,
                       V.band.a = V.band, Tr.band.a = Tr.band,
                       penalty = TRUE, lambda = lambda, alpha.grid = alpha.grid,
                       adjust.sigma = TRUE)    

save(SCC_COMP_1, file = "SCC_COMP_1.RData")

################

SCC_COMP_32 = scc.image(Ya = Y_AD_32, Yb = Y_CN, Z = Z, d.est = d.est, d.band = d.band, r = r,
                       V.est.a = V.est, Tr.est.a = Tr.est,
                       V.band.a = V.band, Tr.band.a = Tr.band,
                       penalty = TRUE, lambda = lambda, alpha.grid = alpha.grid,
                       adjust.sigma = TRUE)    

save(SCC_COMP_32, file = "SCC_COMP_32.RData")

################

SCC_COMP_79 = scc.image(Ya = Y_AD_79, Yb = Y_CN, Z = Z, d.est = d.est, d.band = d.band, r = r,
                       V.est.a = V.est, Tr.est.a = Tr.est,
                       V.band.a = V.band, Tr.band.a = Tr.band,
                       penalty = TRUE, lambda = lambda, alpha.grid = alpha.grid,
                       adjust.sigma = TRUE)    

save(SCC_COMP_79, file = "SCC_COMP_79.RData")

################

SCC_COMP_214 = scc.image(Ya = Y_AD_214, Yb = Y_CN, Z = Z, d.est = d.est, d.band = d.band, r = r,
                       V.est.a = V.est, Tr.est.a = Tr.est,
                       V.band.a = V.band, Tr.band.a = Tr.band,
                       penalty = TRUE, lambda = lambda, alpha.grid = alpha.grid,
                       adjust.sigma = TRUE)    

save(SCC_COMP_214, file = "SCC_COMP_214.RData")

################

SCC_COMP_271 = scc.image(Ya = Y_AD_271, Yb = Y_CN, Z = Z, d.est = d.est, d.band = d.band, r = r,
                       V.est.a = V.est, Tr.est.a = Tr.est,
                       V.band.a = V.band, Tr.band.a = Tr.band,
                       penalty = TRUE, lambda = lambda, alpha.grid = alpha.grid,
                       adjust.sigma = TRUE)    

save(SCC_COMP_271, file = "SCC_COMP_271.RData")

################

SCC_COMP_413 = scc.image(Ya = Y_AD_413, Yb = Y_CN, Z = Z, d.est = d.est, d.band = d.band, r = r,
                       V.est.a = V.est, Tr.est.a = Tr.est,
                       V.band.a = V.band, Tr.band.a = Tr.band,
                       penalty = TRUE, lambda = lambda, alpha.grid = alpha.grid,
                       adjust.sigma = TRUE)    

save(SCC_COMP_413, file = "SCC_COMP_413.RData")

################

# With Yb included we now have a Two-group SCC for the difference between mean functions of Yb and Ya


    # Error in chol.default(lhs) : la submatriz de orden 1 no es definida positiva
    # The error you are seeing occurs when some of the eigenvectors of the matrix you are trying to operate on are not positive
    # Regularizing means may be the solution.
    # Alright, that was from a post in stackoverflow. Let's try a gran mean scaling in order to center our data.

## IMPORTANT NOTE: ESTIMATED MEAN FUNCTIONS APPEAR IN THE ORDER THEY APPEAR IN THE CODE (FIRST Ya THEN Yb)
## BUT SCC IS ACTUALLY CALCULATED FOR THE DIFFERENCE BETWEEN Yb MINUS Ya

    
    # FIRST WAY OF VISUALIZATION:

    plot(SCC_COMP_1,
         # breaks=c(0,2),
         # col="turquoise",
         breaks=seq(from=0,to=2,length.out = 65),
         xlab="Longitudinal (1-95)",
         ylab="Transversal (1-79)",
         sub="Difference between estimated mean functions: CNs - ADs",
         col.sub="red",
         family ="serif")

    # SECOND WAY: JUST ONE COLOR AND THEN WE OVERLAY A SERIES OF POINTS
    # RUN THIS NEXT PLOT CODE, STOP IN ONE OF THE ESTIMATED MEAN FUNCTIONS, THEN RUN "POINTS" TO OVERLAY THEM
    
    plot(SCC_COMP_1,
         breaks=c(0,3),
         col="turquoise",
         # breaks=seq(from=0,to=2,length.out = 65),
         xlab="Longitudinal (1-95)",
         ylab="Transversal (1-79)",
         sub="Difference between estimated mean functions: CNs - ADs",
         col.sub="red",
         family ="serif")


    
my_points <- function(aa){
  
  Z.band <- matrix(aa$Z.band,ncol=2) # Positions
  z1 <- unique(Z.band[,1]); z2 <- unique(Z.band[,2]); # Separated positions
  n1 <- length(z1); n2 <- length(z2) # Lengths of those positions
  scc <- matrix(NA,n1*n2,2) # Matrix with value of that SCC in that position
  ind.inside.band <- aa$ind.inside.cover # Keep only regions inside triangulation
  scc[ind.inside.band,] <- aa$scc[,,2] # Assign SCC to those areas
  scc.limit <- c(min(scc[,1],na.rm=TRUE), max(scc[,2],na.rm=TRUE)) # LIMITS: minimum of inferior, maximum of superior
  
  scc.l.mtx <- matrix(scc[,1],nrow=n2,ncol=n1) # Lower SCC for each location. 
  scc.u.mtx <- matrix(scc[,2],nrow=n2,ncol=n1) # Upper SCC for each location.
  scc.l.mtx[scc.l.mtx<0]=NA # The ones that work just fine are substituded by a NA as we don't want to represent them, only positive LowerSCCs are represented
  scc.u.mtx[scc.u.mtx>0]=NA # The ones that work just fine are substituded by a NA as we don't want to represent them, only positive UpperSCCs are represented
  
  points.P<-which(scc.l.mtx>0,arr.ind=TRUE) # Points where mean difference is positive (first image is stronger)
  points.N<-which(scc.u.mtx<0,arr.ind=TRUE) # Points where mean difference is negative (second image is stronger)
  
  pointers<-list(points.P,points.N)
  print(pointers)  
}

points_1 <- my_points(SCC_COMP_1) # Returns coordinates of points above or below estimated mean function (points.P,points.N; in that order)

## SO NOW IF YOU GO BACK, PLOT THE SCC YOU NEED AND THEN RUN THE FOLLOWING LINES UPON IT, IT WILL DRAW THE POINTS WHICH ARE UP OR DOWN SCC

  points(points_1[[1]],
         type="p",
         pch=20,
         col="navy")  
  
  points(points_1[[2]],
         type="p",
         pch=15,
         col="yellow") 

    # Type: p,l,b,c,o,s,S,h,n
    # pch = 0:18 =:46
    # col=
    # bg= background
    # lwd= line width
