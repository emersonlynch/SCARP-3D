README for 3D offset calculations
By Emerson Lynch
Updated March 4, 2022

This repository contains everything needed to recreate the 3D offset calculations of Lynch et al. (in prep)

In order to calculate the 3D offset of your own geomorphic piercing points, you need
1) topographic profiles as XYZ points (can be surveyed in the field or extracted from a DEM)
2) strike and dip of your fault plane(s)
3) selections of points to include in regressions

I have included a script to calculate the strike and dip of a fault plane from surveyed scarp midpoints: fault-plane-regression.R

To account for the uncertainty in profile selection, I had other members of my lab select the points they would include in each regression. The script to share 3D plots of your profiles is here: profile-viewer.Rmd
>> knit to html to share with others
>> instructions are included for others to install the required tools for viewing the 3D plots

The script to calculate 3D offset is here: 3d_regression.R

I have included the raw profiles and subgroups of fault midpoints included in my paper as txt files in the [topo-profiles] and [fault-midpoints] folders