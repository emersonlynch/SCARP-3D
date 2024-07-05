# SCARP-3D
# Script Calculating displAcements of lineaR Profiles in 3D
# Emerson Lynch
# Last updated March 29, 2024
# cleaned up comments
# added point uncertainty to linear regression

# Known issues:
# RL vs LL currently determined by direction of delta easting - ONLY FOR FAULTS WITH SIMILAR ORIENTATIONS AS BRF
# Currently defines downhill-up as positive dip slip, downhill-down as negative dip slip

##### Introduction --------
# this script takes a total station profile along a geomorphic piercing point
# and calculates a linear regression along each profile segment. Using a fault
# plane with a specified strike and dip with a surveyed midpoint as the origin,
# we calculate the intersections of the profile segments with the fault plane.
# We then calculate heave, throw, dip slip, strike slip, and the true displacement
# vector.

# This is repeated for a range of S/D values within uncertainty and for multiple
# selections of the upper and lower profiles.

# commands to plot at each step are included as comments to avoid generating hundreds of plots with each run

#####  R Set-up Section -------------
# set working directory to folder where script & profiles are kept
setwd("~/file-path")

# install packages if needed. on macs, xquartz must also be installed separately for rgl to work
# install.packages("rgl")#rgl is the 3D visualization package, only needed to produce plots while checking/debugging
# install.packages("matlib")

library(rgl) 
library(matlib) 

rm(list=ls()) #clear workspace

##### Model inputs -------
####'[EDIT THIS SECTION WITH YOUR CHOSEN VALUES]
rm(list=ls()) #clear workspace if you have already run script in this session

# prepared survey csv or txt file (with extension)
og_survey <- read.csv(file="./profile.txt", header=TRUE, sep=",") 

profile = "profile" # name of profile

og_strike_deg <- 273 #'[enter RHR strike

s_sd <- 5 # enter standard deviation to use for strike selection
og_dip_deg <- 72 #'[enter dip

og_dip_deg <- og_dip_deg + 5 #add 5 degrees to correct for underestimating dip from surveyed midpoints
d_sd <- 5 # enter standard deviation to use for dip selection


n <- 9 #column  with the final profile selection
i <- c(6:n) # run for each column of individual choices

unit_length <- 10 # choose a length for the unit vector used in fault plane

all_runs=data.frame() #create a blank data frame to store your outputs

##### Model run ------

for (val in i){

  j = 0
  repeat{
    j = j+1 
    if (j == 101){
      break
    }
    strike_deg <- rnorm(1,og_strike_deg,sd=s_sd) # choose 1 random strike from 
    # a normal distribution with mean strike_deg, standard deviation 5 degrees
    
    dip_deg <- rnorm(1, og_dip_deg, sd=d_sd) # choose 1 random dip from a normal distribution 
    # with mean dip_deg, standard deviation 5 degrees
      if (dip_deg>90){
        dip_diff=dip_deg-90
        dip_deg=90-dip_diff
        
        strike_deg = strike_deg +180
          if (strike_deg>360){
            strike_deg = strike_deg-360
          }
      }
     # surveyed midpoint for plotting profile location later
     og_midpt = subset(og_survey, og_survey[, 5]=="M")
     # midpoint selection by different people
     init_midpt = subset(og_survey, og_survey[, val] == "M") 
     
     survey=og_survey
     survey$Easting<-(og_survey$Easting - init_midpt$Easting)
     survey$Northing<-(og_survey$Northing - init_midpt$Northing)
     survey$HAE<-(og_survey$HAE - init_midpt$HAE)
     
     #redefine midpoint as origin
     midpt = subset(survey, survey[, 5] == "M")

     ## how to add in uncertainty of points? world's clunkiest monte carlo EML 29 Mar 2024
     new_survey = survey
     new_survey$Easting <-rnorm(nrow(survey),survey$Easting,survey$E.error)
     new_survey$Northing <-rnorm(nrow(survey),survey$Northing,survey$N.error)
     
    # #plot the profile in 3D space, with option to click and drag to spin the plot
    # plot3d(x = survey$Easting, y = survey$Northing, z = survey$HAE,aspect="iso",
    #        xlab="Easting (m)",ylab="Northing (m)",zlab="Elevation (m)")
    # #aspect="iso" means isometric - E/N/Z are all to scale - no vertical exaggeration

    ##### PROFILE REGRESSION -----
    # define sections of the profile separated by fault strands
    sect1=subset(survey, survey[,val] == "U") # uphill section
    sect2=subset(survey, survey[,val] == "L") # downhill section
    
    # create XY (N,E) line to use for HAE prediction
    # predict easting as a function of northing
    sect1xy=lm(Easting ~ Northing, data = sect1)
    sect2xy=lm(Easting ~ Northing, data = sect2)
    
    # predict plane of HAE (Height Above Ellipsoid, aka elevation) based on Northing and Easting
    sect1lm=lm(HAE ~ Easting + Northing, data = sect1)
    sect2lm=lm(HAE ~ Easting + Northing, data = sect2)
    
    # sequence 100 new northing points to use to calculate easting, 
    # within the northing bounds of your data
    N_new1=seq(min(survey$Northing),max(survey$Northing),length.out = 100)
    new1=data.frame(Northing=N_new1)
    N_new2=seq(min(survey$Northing),max(survey$Northing),length.out = 100)
    new2=data.frame(Northing=N_new2)
    
    # predict easting values for each new northing point, using your xy linear model for each section
    E_new1 = predict(sect1xy,newdata = new1)
    new1$Easting=E_new1
    E_new2 = predict(sect2xy,newdata = new2)
    new2$Easting=E_new2
    
    # predict elevation values for each new E,N point, using your HAE regression plane
    Z_new1=predict(sect1lm,newdata = new1)
    new1$HAE=Z_new1
    Z_new2=predict(sect2lm,newdata = new2)
    new2$HAE=Z_new2
    
    # # plot these newly generated lines
    # lines3d(new1$Easting,new1$Northing,new1$HAE, col="blue")
    # lines3d(new2$Easting,new2$Northing,new2$HAE, col="red")

    ##### FAULT SECTION ------
    # unit conversions
    strike_rad = strike_deg * pi/180;  # convert to radians
    dip_rad = dip_deg * pi/180;
    
    dd_deg = strike_deg+90; 
    if (dd_deg > 360){
      dd_deg = dd_deg-360;
    }
    dd_rad = dd_deg*pi/180;
    
    # Convert a line of strike to x,y,x data points 
    
    # by definition the z values on a line of strike are the same. We have 
    # to assume a starting elevation. This value is given the name init_strike_z.
    # In this example I assume they are zero. For a scarp profile you could
    # assume this is the XYZ coordinates of the midpoint. 
    # (For my profiles, I define the origin as the midpoint so still zero)
    
    init_strike_z = 0;
    
    #  We also have to assume the starting x,y coordinates of the strike line.
    #  These variable names are called init_strike_x, and init_strike_y.
    #  In this example I assume the initial x,y, coordinates are (0,0)
    
    init_strike_x = 0;
    init_strike_y = 0;
    
    #  Therefore the coordinates of Point 1 (P1) are
    P1 = data.frame(Easting=init_strike_x, Northing=init_strike_y, HAE=init_strike_z)
    
    # For a plane of known strike, if you start at position (init_strike_x,
    # init_strike_y), and move a distance (unit_length), along the line of
    # strike, by trig, you will have changed your distance in x and y by:
    delta_x = sin(strike_rad)*unit_length;
    delta_y = cos(strike_rad)*unit_length;
    
    # now you can add delta_x and _y to init_strike_x and _y to obtain the new
    # coordinates of the end of the line.  
    
    ## this is Point 2 (P2) 
    P2 = data.frame(Easting=(init_strike_x + delta_x), Northing=(init_strike_y + delta_y), 
                    HAE=init_strike_z)
    
    # Calculate points down-dip 
    # Calculate X and Y coordinates of point P3 by moving out a distance 
    # unit_length along a the dip direction from P1
    dipdir_delta_x = sin(dd_rad)*unit_length;
    dipdir_delta_y = cos(dd_rad)*unit_length;
    
    P3X = P1[1]+dipdir_delta_x;
    P3Y = P1[2]+dipdir_delta_y;
    
    # calculate the z value that corresponds to the point on the plane 
    # below P3X,P3Y, which we can calculate given the fault dip
    delta_z = tan(dip_rad) * unit_length;
    P3Z = init_strike_z - delta_z;
    P3 = data.frame(Easting=P3X, Northing=P3Y, HAE=P3Z)
    
    # calculate plane 
    fault_points=data.frame(P1)
    fault_points<-rbind(fault_points,P2,P3)
    
    fault.lm=lm(Easting~HAE+Northing,data=fault_points)
    coefs <- coef(fault.lm)
    
    # set up to plot fault plane 
    a <- -1
    b <- coefs["Northing"]
    c <- coefs["HAE"]
    d <- coefs["(Intercept)"]

#     planes3d(a, b, c, d, alpha = 0.25)
#     grid3d(c("x","y","z"))

    ##### INTERSECTIONS ------
    coefs1 <- coef(sect1xy)
    coefs2 <- coef(sect1lm)
    coefs3 <- setNames(c(a,b,c,d), c("Easting", "Northing","HAE","(Intercept)"))
    
    A<-matrix(c(1,-coefs1["Northing"],0,
                -coefs2["Easting"],-coefs2["Northing"],1,
                -coefs3["Easting"],-coefs3["Northing"],-coefs3["HAE"]),3,3,byrow=TRUE)
    colnames(A)<-c("x","y","z")
    B<-matrix(c(coefs1["(Intercept)"],coefs2["(Intercept)"],coefs3["(Intercept)"]))
    PtB<-solve(A, B)
    
    #intersection of downhill segment with fault plane
    coefs1 <- coef(sect2xy)
    coefs2 <- coef(sect2lm)
    
    A<-matrix(c(1,-coefs1["Northing"],0,
                -coefs2["Easting"],-coefs2["Northing"],1,
                -coefs3["Easting"],-coefs3["Northing"],-coefs3["HAE"]),3,3,byrow=TRUE)
    colnames(A)<-c("x","y","z")
    B<-matrix(c(coefs1["(Intercept)"],coefs2["(Intercept)"],coefs3["(Intercept)"]))
    PtA<-solve(A, B)
    
    # #plot intersections!
    # points3d(t(PtB),col="blue",size=6,label="B")
    # points3d(t(PtA),col="red", size=6,label="A")

    ##### CALCULATIONS -----
    # Calculate throw
    throw = PtA[3] - PtB[3]
    
    # Calculate dip slip
    dip_slip = throw/sin(dip_rad)
    
    #Calculate heave
    heave = throw/tan(dip_rad)
    heave2 = sqrt(dip_slip^2 - throw^2)
    
    #Calculate Point C
    #Starting with the upper intersection point (Int2 for uphill facing scarps), 
    #plot a point that is dip_slip away down dip
    #generate E, N coords for a point dip_slip away along dip direction
    
    deltaE=sin(dd_rad)*heave
    deltaN=cos(dd_rad)*heave
    
    
    ePtC = PtA[1] + deltaE
    nPtC = PtA[2] + deltaN
    
    PtC=data.frame(Easting=ePtC)
    PtC$Northing=nPtC
    #calculate Z of that point using plane equation
    
    zPtC = (-a*PtC$Easting - b*PtC$Northing - d)/(c)
    PtC$HAE=zPtC
    
    # # plot point C
    # points3d(x=PtC$Easting,y=PtC$Northing,z=PtC$HAE,size=6,color="green")
    
    ## Calculate strike slip
    # calculate delta E
    delta_E <- PtB[1]-PtC$Easting
    #calculate delta N
    delta_N <-PtB[2]-PtC$Northing
    
    s_slip <- sqrt(delta_E^2 + delta_N^2)
      
      ## calculate direction of strike slip
    dE = PtB[1] - PtA[1]
    dN = PtB[2] - PtA[2]
    dZ = PtB[3] - PtA[3]
    alpha_rad = atan(dE/dN)
    alpha_deg = alpha_rad * (180/pi)
    if (dN > 0){
      if (dE > 0){
        trend = alpha_deg
      }
      if (dE < 0){
        trend = 360 - abs(alpha_deg)
      }
    }
    if (dN < 0){
      if (dE > 0){
        trend = 180 - abs(alpha_deg)
      }
      if (dE < 0){
        trend = 180 + abs(alpha_deg)
      }
    }
    
    ###FOR determining RL (+) vs LL (-) for BRF ONLY###
    if (dE < 0) {
      s_slip = -1 * s_slip
    }
    
    ###FIX THIS LATER TO USE WITH OTHER FAULT SYSTEMS###
    # # determine sinistral vs dextral. defining convention of (+) for RL and (-) for LL
    # # this method can return false LL when planes dip SW
    # if (trend > dd_deg){
    #   theta = trend - dd_deg
    #   if (90 < theta){
    #     if (theta < 180){
    #       s_slip = -1 * s_slip
    #     }
    #   }
    #   if (270 < theta){
    #     if (theta < 360){
    #       s_slip = -1 * s_slip
    #     }
    #   }
    # }
    # 
    # if (trend < dd_deg){
    #   theta = dd_deg - trend
    #   if (0 < theta){
    #     if (theta < 90){
    #       s_slip = -1 * s_slip
    #     }
    #   }
    #   if (180 < theta){
    #     if (theta < 270){
    #       s_slip = -1 * s_slip
    #     }
    #   }
    # }
    
    ## Calculate true displacement
    true_disp <- sqrt(s_slip^2 + dip_slip^2)
    
    # calculate plunge, for completion
    dist = sqrt(dE^2 + dN^2)
    plunge_rad = atan(-dZ/dist) # minus dZ because plunge is measured down from horizontal
    plunge_deg = plunge_rad * (180/pi)
   
    
    ## Calculate the ratio of dip slip to strike slip
    ratio <- s_slip/dip_slip
    
    ##### SAVE VALUES TO DATA FRAME -----
    name = colnames(survey[val])
    
    modelrun <-data.frame(strike_deg,dip_deg,heave,throw,dip_slip,
                          s_slip,true_disp,trend,plunge_deg,ratio,dE,dN)
    row.names(modelrun)<-name
    
    all_runs <- rbind(all_runs,modelrun)
  }
 }

##### SUMMARY STATS - RUN ONCE AT END OF FOR LOOP -----
####'[EDIT THIS SECTION WITH YOUR CHOSEN VALUES]
profile <- all_runs #'[enter the name of the profile]

save(profile,file="./profile.Rda") #'[re-enter the name of the profile as 
  #'[the first arg, specify file location/name in second]

###

heave_m = mean(all_runs$heave)
heave_sd = sd(all_runs$heave)

throw_m = mean(all_runs$throw)
throw_sd = sd(all_runs$throw)

dip_slip_m = mean(all_runs$dip_slip)
dip_slip_sd = sd(all_runs$dip_slip)

s_slip_m = mean(all_runs$s_slip)
s_slip_sd = sd(all_runs$s_slip)

true_disp_m = mean(all_runs$true_disp)
true_disp_sd = sd(all_runs$true_disp)

ratio_m = mean(all_runs$ratio)
ratio_sd = sd(all_runs$ratio)

trend_m = mean(all_runs$trend)
trend_sd = sd(all_runs$trend)

plunge_m = mean(all_runs$plunge_deg)
plunge_sd = sd(all_runs$plunge_deg)


##### plot profile & fault to check offsets if needed #####
# these are the commands throughout the script condensed into one place here
# #plot the profile in 3D space, with option to click and drag to spin the plot
# plot3d(x = survey$Easting, y = survey$Northing, z = survey$HAE,aspect="iso",
#        xlab="Easting (m)",ylab="Northing (m)",zlab="Elevation (m)")
# #aspect="iso" means isometric - E/N/Z are all to scale - no vertical exaggeration
# # plot these newly generated lines
# lines3d(new1$Easting,new1$Northing,new1$HAE, col="blue")
# lines3d(new2$Easting,new2$Northing,new2$HAE, col="red")
#     planes3d(a, b, c, d, alpha = 0.25)
#     grid3d(c("x","y","z"))
# #plot intersections!
# points3d(t(PtB),col="blue",size=6,label="B")
# points3d(t(PtA),col="red", size=6,label="A")
# # plot point C
# points3d(x=PtC$Easting,y=PtC$Northing,z=PtC$HAE,size=6,color="green")

#
##### first time running script create data frame for fault regressions-----
offsets<-data.frame(dip_slip_m,dip_slip_sd,heave_m,heave_sd,throw_m,throw_sd,s_slip_m,
                    s_slip_sd,true_disp_m,true_disp_sd,ratio_m,ratio_sd,trend_m,trend_sd,plunge_m,plunge_sd)
row.names(offsets)<-profile

#save data frame to be able to load later in 3D regression
save(offsets,file="./offsets.Rda") #'[give file a useful name]

##### save offsets on subsequent runs -------
load("./offsets.Rda")
new_offsets<-data.frame(dip_slip_m,dip_slip_sd,heave_m,heave_sd,
                        throw_m,throw_sd,s_slip_m,s_slip_sd,true_disp_m,
                        true_disp_sd,ratio_m,ratio_sd,trend_m,trend_sd,plunge_m,plunge_sd)
row.names(new_offsets)<-profile
offsets<-rbind(offsets,new_offsets)

save(offsets,file="./offsets.Rda")

# if you want to save your file to csv to open in another program
write.csv(offsets,file="./offsets.csv")
