########################################################################\
#   
#   THE OBJECTIVES OF THIS SCRIPT ARE:
#     
#     - TO IMPORT SIMULATED CHUS/USC NIFTI FILES INTO 'R' (NEW DIMENSIONS)
#     - TO TEST PREVIOUS SCRIPTS ON SCC WITH THESE IMAGES      
#     - IF EVERYTHING WORKS FINE, ESCALATE THIS COMPARISON
#
#######################################################################\


### ################################################## ###
#####                NEUROIMAGE DATA                  ####
### ################################################## ###


## PREVIOUS PROCESSING OF PET IMAGES (REALIGNEMENT, WRAPPING, CORREGISTRARION, NORMALIZATION, AND MASKING) DONE 


### ################################################## ###
#####                 *PREAMBLE*                      ####
####         Installation of necessary packgs.         ###
### ################################################## ###


install.packages(c("mgcv","gamair","oro.nifti","memsic"))
library(mgcv);library(gamair);library(oro.nifti);library(memisc)
 
# If necessary:

    # install.packages("unix") 
    # library(unix)
    # rlimit_as(1e12)  #increases to ~12GB
    # rlimit_all()


CNimg <- readNIfTI("w00", verbose = FALSE, warn = -1, reorient = TRUE, call = NULL, read_data = TRUE)

orthographic(CNimg,text="Plano Ortografico PET", col.y = hotmetal(),
             crosshairs = TRUE, col.crosshairs = "red",)

ADimg <- readNIfTI("w413_0_8", verbose = FALSE, warn = -1, reorient = TRUE, call = NULL, read_data = TRUE)

orthographic(ADimg,text="Plano Ortografico PET", col.y = hotmetal(),
             crosshairs = TRUE, col.crosshairs = "red",)


### ################################################## ###
#####            *IMPORT NIFTI FILES*                 ####
####            f.clean  (PPT,z,x,y,pet)               ###
### ################################################## ###


f.clean <- function(name) {
  
  ## Load Data
  
  file <- readNIfTI(fname = name, verbose = FALSE, warn = -1, reorient = TRUE, call = NULL, read_data = TRUE)
  namex <- as.character(name)
  n = img_data(file)
  n = to.data.frame(n)

  dataframe <- data.frame(z = integer(),x = integer(),y = integer(),pet = integer()) 
  
  # Loop for 91 slices of Z in the NiFtI image -> move to dataframe
  
  for (i in seq(1:91)) {
    
    n_lim = n[n$Var2==i,] # Select just one Z slice
    n_lim$Var1=NULL
    n_lim$Var2=NULL
    
    z <- rep(i, length.out = 9919)
    x <- rep(1:91, each = 109, length.out = 9919) 
    y <- rep(1:109, length.out = 9919)
    
    attach(n_lim)
    pet <- c(`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,`9`,`10`,`11`,`12`,`13`,`14`,
             `15`,`16`,`17`,`18`,`19`,`20`,`21`,`22`,`23`,`24`,`25`,`26`,
             `27`,`28`,`29`,`30`,`31`,`32`,`33`,`34`,`35`,`36`,`37`,`38`,
             `39`,`40`,`41`,`42`,`43`,`44`,`45`,`46`,`47`,`48`,`49`,`50`,
             `51`,`52`,`53`,`54`,`55`,`56`,`57`,`58`,`59`,`60`,`61`,`62`,
             `63`,`64`,`65`,`66`,`67`,`68`,`69`,`70`,`71`,`72`,`73`,`74`,
             `75`,`76`,`77`,`78`,`79`,`80`,`81`,`82`,`83`,`84`,`85`,`86`,
             `87`,`88`,`89`,`90`,`91`)
    detach(n_lim)
    
    temp0 = data.frame(z,x,y,pet) # temporal dataframe
    temp1 <- print(temp0) 
    dataframe <- rbind(dataframe,temp1)
  }
  
  # Demographics: PPT, group (AD/CN), sex, age.
  
  # demog <- demo[demo$PPT==namex,]
  # 
  # PPT <- rep(demog$PPT,length.out=7505)
  # group <-rep(demog$Group,length.out=7505)
  # sex <-rep(demog$Sex,length.out=7505)
  # age <-rep(demog$Age,length.out=7505)
  # temp2 <- data.frame(PPT,group,sex,age)
 
  temp2 <- rep("AD", length.out = 9919)
  dataframe <- cbind(temp2,dataframe)
  
  print(dataframe) # Necessary for assigning an object name

}

CN = f.clean("w00")
AD = f.clean("w413_0_8")
example_data = rbind(CN, AD)

database = example_data 

nrow(database[database$pet<0,]) # There are negative values (ilogical)

database$pet[database$pet<0]<- NaN  # Convert negatives to NaN 
database$pet[database$pet==0] <- NaN # Convert zeros to NaN 


### ########################################################## ###
#####          SIMULTANEOUS CONFIDENCE CORRIDORS              ####
### ########################################################## ###

# This part of the code needs some packages which to date (09/2020) are only available at GitHub

install.packages(c("devtools","remotes","readr","imager","itsadug","ggplot2","contoureR","fields"))

library(devtools);library(remotes);library(readr);library(imager);library(itsadug);library(ggplot2);library(contoureR);library(fields)

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE) # Forces installation to finish even with warning messages (EXTREMELY NECESSARY)
remotes::install_github("funstatpackages/BPST")
remotes::install_github("funstatpackages/Triangulation")
remotes::install_github("funstatpackages/ImageSCC")

library(BPST);library(Triangulation);library(ImageSCC)

Data<-database
attach(Data)
head(Data)
str(Data)

### ########################################################## ###
#####           *FUNCTIONS FOR SCCs DATA.FRAME*               ####
### ########################################################## ###


  Y <- subset(Data, Data$temp2=="CN" & Data$z==35) 
  Y <- Y[1:9919,5] 
  Y <- as.matrix(Y)
  Y = t(Y) 
  Y[is.nan(Y)] <- 0
  SCC_CN <- Y

  Y <- subset(Data, Data$temp2=="AD" & Data$z==35) 
  Y <- Y[1:9919,5] 
  Y <- as.matrix(Y)
  Y = t(Y) 
  Y[is.nan(Y)] <- 0
  SCC_AD <- Y
  
  # Use z=35 as alternative to z=30

### ########################################################## ###
#####                *CONTOURS OF NEURO-DATA*                 ####
### ########################################################## ###
        

# Z are the coordinates where data is measured:

x <- rep(1:91, each = 109, length.out = 9919) 
y <- rep(1:109,length.out = 9919)
Z <- cbind(as.matrix(x),as.matrix(y)); Z

dat <- cbind(Z,t(SCC_AD))

dat <- as.data.frame(dat)
dat[is.na(dat)] <- 0
sum(is.na(dat$pet)) # should be = 0
head(dat)

rownames(dat) <- NULL 

# df = getContourLines(dat[1:9919,], levels = c(0)) 
# No se ha hecho un enmascaramiento de los datos, por lo que el contour está
# mal. Voy a usar df_test temporalmente.

library(contoureR); library(ggplot2)

df_test = getContourLines(dat, levels = c(2))  

ggplot(df_test, aes(x, y, colour = z)) + geom_path() 

contour35 = df_test 
head(contour35);str(contour35)

contour <- function(x){
  
  aa <- contour35[contour35$GID == x,] # We keep GID==x (0,1,2,3...)
  a <- aa[, 5:6] # Then we keep just the coordinates 
  print(a) # and then print in order to loop and make a list
}

coord <- list()

for (i in 0:max(contour35$GID)) { #change contour30 to any other name previously assigned if necessary
   
  coord[[i + 1]] <- contour(i)
  rownames(coord[[ i + 1 ]]) <- NULL
}


### Test the results are coherent:

head(coord[[1]],10); plot(coord[[1]])   # external boundaries
head(coord[[2]],10); points(coord[[2]]) # first hole
head(coord[[3]],10); points(coord[[3]]) # second hole
# (...) 


### ########################################################## ###
#####             *TRIANGULATION PARAMETERS*                  ####
### ########################################################## ###

VT = TriMesh(coord[[1]], 8)

# VTtry = TriMesh(coord[[1]],8,list(as.matrix(coord[[2]])))
#              
# VT8 = TriMesh(coord[[1]], 8,
#               list(as.matrix(coord[[2]]),
#                    as.matrix(coord[[3]]), 
#                    as.matrix(coord[[4]]), 
#                    as.matrix(coord[[5]])))

head(VT8$V,10);head(VT8$Tr,10) 


# In order to be consistent we use common names Brain.V and Brain.Tr. From here onwards most of the names follow the ones provided by Wang et al (2019)

Brain.V <- VT[[1]]
Brain.Tr <- VT[[2]]

head(Brain.V);head(Brain.Tr)


V.est=as.matrix(Brain.V)
# Brain.v<-cbind(Brain.V[,2],Brain.V[,1]) # In case you need to transpose the data
Tr.est=as.matrix(Brain.Tr)

V.band=as.matrix(Brain.V)
Tr.band=as.matrix(Brain.Tr) 


# Response Variable:

Y_CN = SCC_CN
Y_AD = SCC_AD


### ########################################################## ###
#####        *OTHER PARAMETERS FOR SCC ESTIMATION*            ####
### ########################################################## ###


# Following Wang et al's recomendations:

d.est = 7 # degree of spline for mean function  5
d.band = 4 # degree of spline for SCC  2
r = 2 # smoothing parameter  1
lambda = 10^{seq(-6,3,0.5)} # penalty parameters
alpha.grid = c(0.10,0.05,0.01) # vector of confidence levels



### ########################################################## ###
#####               *CONSTRUCTION OF SCC'S*                   ####
### ########################################################## ###

Y_AD <- rbind(Y_AD,Y_AD)
Y_CN <- rbind(Y_CN,Y_CN)


# Run one sample SCC construction:

library(unix)
rlimit_as(1e12)  #increases to ~12GB
rlimit_all()

SCC_CN_1 = scc.image(Ya = Y_CN,Z=Z,d.est=d.est,d.band=d.band,
                     r=r,V.est.a=V.est,Tr.est.a=Tr.est,
                     V.band.a=V.band,Tr.band.a=Tr.band,
                     penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,
                     adjust.sigma=TRUE)

  plot(SCC_CN_1,
       breaks = seq(from = 0,to = 2,length.out = 65),
       col = ,
       xlab = "Longitudinal (1-95)",
       ylab = "Transversal (1-79)",
       sub = "Control Group",
       col.sub = "black",
       family = "serif")
  


SCC_AD_1=scc.image(Ya=Y_AD,Z=Z,d.est=d.est,d.band=d.band,r=r,
                   V.est.a=V.est,Tr.est.a=Tr.est,
                   V.band.a=V.band,Tr.band.a=Tr.band,
                   penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,
                   adjust.sigma=TRUE)


  plot(SCC_AD_1,
       breaks=seq(from=0,to=2,length.out = 65),
       col=,
       xlab="Longitudinal (1-95)",
       ylab="Transversal (1-79)",
       sub="Alzheimer Group",
       col.sub="black",
       family ="serif")

   # Same here for the other group (AD). These images are orientative, our main goal is obtaining SCCs for the difference
   # between groups' mean functions.
  
   ## IN CASE YOU WANT TO EXPORT PLOTS THIS IS EASIER AFTER (...) :  
  
      par(las=1,
          col.main="white",
          col.lab="white",
          col.sub="white")
      
      plot(SCC_AD_1,
           axes=FALSE,
           breaks=seq(from=0,to=2,length.out = 65),
           col=,
           family ="serif",
           horizontal=F)
      dev.off()


      
SCC_COMP_1=scc.image(Ya=Y_AD,Yb=Y_CN,Z=Z,d.est=d.est,d.band=d.band,r=r,
                     V.est.a=V.est,Tr.est.a=Tr.est,
                     V.band.a=V.band,Tr.band.a=Tr.band,
                     penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,
                     adjust.sigma=TRUE)    


# With Yb included we now have a Two-group SCC for the difference between mean functions of Yb and Ya

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
         breaks=c(0,2),
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


  # NOT NECESSARY TO RUN THESE TWO LINES
  image.plot(z2,z1,scc.l.mtx, zlim = scc.limit) # Regions where the difference between mean functions is positive (falls above 0). That is, activity in Image1 is stronger than in Image 2 in that area
  image.plot(z2,z1,scc.u.mtx, zlim = scc.limit) # Regions where the difference between mean functions is negative (falls below 0). That is, activity in Image2 is stronger than in Image 1 in that area 
  #

points_1 <- my_points(SCC_COMP_1) # Returns coordinates of points above or below estimated mean function (points.P,points.N; in that order)

## SO NOW IF YOU GO BACK, PLOT THE SCC YOU NEED AND THEN RUN THE FOLLOWING LINES UPON IT, IT WILL DRAW THE POINTS WHICH ARE UP OR DOWN SCC

  points(points_1[[1]],
         type="p",
         pch=15,
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


##############################//


# Percentage of locations where the true difference between mean functions value fall within the SCC

mfit$Yhat <- rbind(SCC_COMP_1$mu.hat.a,SCC_COMP_1$mu.hat.b)

xx <- Z[,1]
uu <- unique(xx)
n1 <- length(uu)
yy <- Z[,2]
vv <- unique(yy)
n2 <- length(vv)
npix <- (as.numeric(Z[nrow(Z),1]))*(as.numeric(Z[nrow(Z),2]));

# Generate true mean function
if(mu.func==1){
  beta.true <- 20*((xx-0.5)^2+(yy-0.5)^2)
}else if(mu.func==2){
  beta.true <- 5*exp(-15*((xx-0.5)^2+(yy-0.5)^2))+0.5
}else if(mu.func==3){
  beta.true <- 3.2*(-xx^3+yy^3)+2.4
}else if(mu.func==4){
  # beta.true <- -10*sin(7*pi*(xx+0.09))+10*sin(7*pi*(yy-0.14))+2.5
  beta.true <- -10*sin(5*pi*(xx+0.22))+10*sin(5*pi*(yy-0.18))+2.8
}
ind.outside <- setdiff(1:npix,ind.inside)
beta.true[ind.outside] <- NA
beta.true <- as.matrix(beta.true)
beta.diff=beta.true[,2]-beta.true[,1]


apply(SCC_COMP_2$scc,3,FUN=function(scc){
  sum((scc[,1]<beta.diff[ind.inside]) & (scc[,2]>beta.diff[ind.inside]))/length(ind.inside)
})


# Check if 0 fall within the SCC everywhere
apply(SCC_COMP_2$scc,3,FUN=function(scc){
  any((scc[,1]>0) | (scc[,2]<0))
})


# Simulation 1000 times to obtain empirical coverage rate and empirical power
nsim=100
coverage=lapply(1:nsim,FUN=function(iter){
  set.seed(iter)
  cat('Iteration No.',iter,'\t')
  ptm0=Sys.time()
  dat=data2g.image(na=na,nb=nb,Z=Z,ind.inside=ind.inside,mu1.func=mu1.func,
                   noise.type='Func',lam1=lam1,lam2=lam2,iter=iter,delta=delta)
  Ya=dat$Ya
  Yb=dat$Yb
  beta.diff=dat$beta.true[,2]-dat$beta.true[,1]
  
  out=scc.image(Ya=Ya,Yb=Yb,Z=Z,V.est.a=V.est.a,Tr.est.a=Tr.est.a,V.band.a=V.band.a,Tr.band.a=Tr.band.a,
                V.est.b=V.est.b,Tr.est.b=Tr.est.b,V.band.b=V.band.b,Tr.band.b=Tr.band.b,
                d.est=d.est,d.band=d.band,r=r,
                penalty=TRUE,lambda=lambda,alpha.grid=alpha.grid,adjust.sigma=TRUE)
  cover.true=apply(out$scc,3,FUN=function(scc){
    sum((scc[,1]<beta.diff[ind.inside]) & (scc[,2]>beta.diff[ind.inside]))/length(ind.inside)
  })
  reject.null=apply(out$scc,3,FUN=function(scc){
    any((scc[,1]>0) | (scc[,2]<0))
  })
  ptm1=Sys.time()
  cat('Time: ',ptm1-ptm0,'\n')
  return(list(bw=out$bw,cover.true=cover.true,reject.null=reject.null))
})

# Empirical test power:
reject.all=sapply(coverage,'[[',3)
apply(reject.all,1,mean)

# Empirical coverage rate
cover.all=sapply(coverage,'[[',2)
apply(cover.all==1,1,mean)

# Mean bandwidth
bw.all=sapply(coverage,'[[',1)
apply(bw.all,1,mean)