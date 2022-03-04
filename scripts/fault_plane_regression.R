# fault plane regression
# Emerson Lynch
# July 1, 2020
# last updated June 22, 2021

# INTRODUCTION --------------------
# This script regresses a fault plane through scarp midpoints surveyed in the field, 
# but the quality of that fault plane is dependent on how many midpoints were surveyed, 
# whether they actually fall along one plane, and the range in elevations surveyed. 
# If there is a small range in elevations, fault dip will not be well constrained.


# RSTUDIO SET-UP: run at beginning of each session --------------- 
# set working directory to folder where script & profiles are kept
setwd("~/file-path")

# install packages if needed. on macs, xquartz must also be installed separately for rgl to work
# install.packages("rgl") #rgl is the 3D visualization package
# install.packages("matlib")
#install.packages("todor")

library(rgl) 
library(matlib) 

# SET-UP: run at beginning of each fault/strand-----------------
rm(list=ls()) #clear workspace
Mids<-read.csv(file="midpoints.txt",header=TRUE,sep=",") #'[read in text file of strand midpoints]

#optional, sets E,N origin to 0
Mids$Easting<-(Mids$Easting - min(Mids[,4], na.rm=T))
Mids$Northing<-(Mids$Northing - min(Mids[,3], na.rm=T))


#plot midpoints 
# open3d() #open a new xquartz window if you already have one open and don't want to overwrite its contents
plot3d(x=Mids$Easting, y=Mids$Northing,z=Mids$HAE,aspect="iso", 
         xlab="Easting (m)",ylab="Northing (m)",zlab="Elevation (m)")

# FAULT PLANE REGRESSION -------------

MidsElm=lm(Easting~HAE+Northing,data=Mids)

coefsE <- coef(MidsElm)
a <- -1
b <- coefsE["Northing"]
c <- coefsE["HAE"]
d <- coefsE["(Intercept)"]
planes3d(a, b, c, d, alpha = 0.25)

# CALCULATE STRIKE ---------------
# let's choose the median Z as our shared elevation between the two points
Z = median(Mids$HAE)
N1 = median(Mids$Northing)
N2 = median(Mids$Northing) + 5 #choose a distance that scales with your data

#calculate the easting of these points based on our plane equation
E1 = (b*N1 + c*Z + d)
points3d(x = E1, y = N1, z = Z, col="blue")

E2 = b*N2 + c*Z + d
points3d(x = E2, y = N2, z = Z, col="red")

#calculate strike from these points (NB: visually check plot to see which of N1 and N2 is larger, same for E)
rise = N2 - N1 #should be the distance you chose above
run = E1 - E2

alpha_rad = atan(rise/run)
strike_rad = alpha_rad + (pi/2) # add 90° to get angle from N=000
strike_deg<-strike_rad*(180/pi)

#'[visually check RHR on the plot! if not...]
strike_deg=strike_deg+180; strike_rad = strike_rad + pi
if (strike_deg>360){
  strike_deg = strike_deg-360
}


# CALCULATE DIP ---------------------
#dip direction
dd_deg=strike_deg+90
if (dd_deg>360){
  dd_deg = dd_deg-360
}

dd_rad=strike_rad + (pi/2)
if (dd_rad > (2*pi)){
  dd_rad = dd_rad-(2*pi)
}

#convert to radians
dd_rad_1=dd_deg*(pi/180)

#calculate dip 
#generate a random N value within your plane
npt1 <-runif(1,min(Mids[,3],na.rm=T),max(Mids[,3],na.rm=T))
pt1=data.frame(Northing=npt1)

#pick the max Z value so you can go down dip within your plotted plane
zpt1 <-max(Mids[,5],na.rm=T)
pt1$HAE=zpt1

#calculate E from N, Z using fault plane eqn
ept1 <- b*pt1$Northing + c*pt1$HAE + d
pt1$Easting=ept1

#plot our trial point
points3d(x=pt1$Easting,y=pt1$Northing,z=pt1$HAE,color="violet")

#let's say our unit vector is 2m
run = 2
#generate E, N coords for a point 2m away along new dip direction
ept2 = ept1 + (sin(dd_rad))*run
pt2=data.frame(Easting=ept2)
npt2 = npt1 + (cos(dd_rad))*run 
pt2$Northing=npt2

#calculate Z of that point using plane equation
zpt2 = (a*pt2$Easting + b*pt2$Northing + d)/(-c)
pt2$HAE=zpt2

# plot our second point
points3d(x=pt2$Easting,y=pt2$Northing,z=pt2$HAE,color="green")

# calculate dip
rise = pt1$HAE - pt2$HAE
dip_rad=atan(rise/run)
dip_deg=dip_rad*(180/pi)

if (dip_rad < (0)){
  dip_rad = dip_rad*(-1)
}
if (dip_deg < (0)){
  dip_deg = dip_deg*(-1)
}

{} #'[blank line, no need to execute - don't want to auto-expand next section

# CREATING THE FAULT DATAFRAME ---------------

# first time running script create data frame for fault regressions
strand<-c('strand') #'[name of first fault strand run]
faults<-data.frame(a,b,c,strike_rad,strike_deg,dd_rad,dd_deg,dip_rad,dip_deg)
row.names(faults)<-strand
names(faults)[names(faults) == "a"] <- "Easting"
names(faults)[names(faults) == "b"] <- "Northing"
names(faults)[names(faults) == "c"] <- "HAE"

#save data frame to be able to load later in 3D regression
save(faults,file="faults.Rda") #'[decide what you want to name your data frame]


# ADDIING TO THE FAULT DATAFRAME ------------------

# subsequent runs (subsequent fault strands) append those data
strand<-c('strand') #'[name of new strand to be appended]
load("~/Documents/faults.Rda") #file path to the data frame you created above
new_fault<-data.frame(a,b,c,strike_rad,strike_deg,dd_rad,dd_deg,dip_rad,dip_deg)
row.names(new_fault)<-strand
names(new_fault)[names(new_fault) == "a"] <- "Easting"
names(new_fault)[names(new_fault) == "b"] <- "Northing"
names(new_fault)[names(new_fault) == "c"] <- "HAE"
faults<-rbind(faults,new_fault)

#save data frame to be able to load later in 3D regression
save(faults,file="faults.Rda") #same data frame name as above

