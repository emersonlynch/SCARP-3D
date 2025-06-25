Three-dimensional offsets of geomorphic piercing lines displaced by the Quaternary-active Beaufort Range fault, northern Cascadia forearc, BC, Canada

Dryad data repository: https://doi.org/10.5061/dryad.nvx0k6dvc


GENERAL INFORMATION
	1.	Title of Dataset: Three-dimensional offsets of geomorphic piercing lines displaced by the Quaternary-active Beaufort Range fault, northern Cascadia forearc, BC, Canada
	2.	Author Information: Principal Investigator Contact Information Name: Emerson M. Lynch; Institution: Northern Arizona University, Flagstaff, AZ, USA; Email: eml365@nau.edu, emerson.lynch@gmail.com, now at USGS Geologic Hazards Science Center, Golden, CO, USA, elynch@usgs.gov. Authors and affiliations: Emerson M. Lynch, Northern Arizona University; Christine Regalla, Northern Arizona University; Kristin D. Morell, UC Santa Barbara; Nicolas Harrichhausen, University of Alaska Anchorage; Lucinda J. Leonard, University of Victoria
	3.	Date of data collection (06/2018-10/2018)
	4.	Geographic location of data collection: Port Albert, British Columbia, Canada
	5.	Information about funding sources that supported the collection of the data: Funding for this project was supported by NSF EAR Tectonics 1756834/ 2004684 to C. Regalla and 1756943 to K. Morell, by the Boston University Department of Earth and Environmental Sciences, and the Northern Arizona University Duebendorfer-Barnes Structure Endowment. 

SHARING/ACCESS INFORMATION
	1.	Licenses/restrictions placed on the data: This work is licensed under a CC0 1.0 Universal (CC0 1.0) Public Domain Dedication license. Files fault-plane-regression.R and SCARP-3D.R are copyright under a GNU GPL v3 or later license.
	2.	Links to publications or works that cite or use the data:
	    ◦	Lynch, E.M., Regalla, C., Morell, K.D., Harrichhausen, N., Leonard, L.J., (in revision for Seismica). Evidence for an active transtensional Beaufort Range fault in the northern Cascadia forearc
	3.	Recommended citation for this dataset:
	    ◦	Lynch, E.M., Regalla, C., Morell, K.D., Harrichhausen, N., Leonard, L.J., (2025) Three-dimensional offsets of geomorphic piercing lines displaced by the Quaternary-active Beaufort Range fault, northern Cascadia forearc, BC, Canada. Dryad, Dataset.

DESCRIPTION OF THE DATA AND FILE STRUCTURE
This dataset contains raw total station survey data, R scripts, and a markdown file.
	1.	Total station survey data: We include two zip files of total station surveys: fault-midpoints.zip and topo-profiles.zip. The survey text files contain the fields Point (integer identifier), Easting (m), Northing (m), HAE (height above ellipsoid, or elevation, m), and Code (see below). Easting and Northing coordinates are in WGS 84 UTM Zone 10.
	    ⁃	fault-midpoints.zip contains surveyed fault midpoints or inflection points used to calculate strike and dip of the fault plane. Each fault segment is in a separate .txt file. The Code field includes an abbreviation of what was surveyed and can be cross-referenced with Table S2 in the accompanying article (e.g., TH1A1 refers to a point collected in the thalweg 1 profile at strand A1).
	    ⁃	topo-profiles.zip contains surveys of geomorphic piercing lines, with each profile in a separate .txt file with names corresponding to Table S2 in the accompanying article. The Code field consists of point classification: P = profile; M = midpoint; S = scarp.
	2.	R scripts:
	    ⁃	fault-plane-regression.R: this script calculates the strike and dip of a fault plane
	    ⁃	required inputs: fault scarp midpoints as xyz points.
	    ⁃	SCARP-3D.R: Script Calculating displAcements of lineaR Profiles in 3D. This script calculates 3D offset of geomorphic piercing points and is hosted on a linked Zenodo repository.
	    ⁃	required inputs: topographic profiles as xyz points (can be surveyed in the field or extracted from a digital elevation model), strike and dip of the fault plane, selected points to include in the regressions
	3.	Markdown file: profile-viewer.Rmd
	    ⁃	this markdown file can be used to share profiles with others to facilitate multiple choices of points to include in linear regressions.
	    ⁃	required inputs: topographic profiles as xyz points with code field indicating scarp midpoints.
	    ⁃	instructions are included for recipients to install required tools for viewing 3D plots in the resulting html file.

METHODOLOGICAL INFORMATION
	1.	Description of methods used for collection/generation of data:
	    ⁃	Total station survey data were collected in summer and fall of 2018 over the course of two field seasons. These data were collected using Spectra Precision Focus 6 and Nikon XS total stations.
	2.	Instrument- or software-specific information needed to interpret the data:
	    ⁃	The programming language R is required, and our scripts were written in the GUI RStudio. We include commands at the beginning of each script to install all required packages. RStudio is available at: https://posit.co/download/rstudio-desktop/
	    ⁃	XQuartz is required for viewing and manipulating 3D plots generated by the R scripts. XQuartz is available at: https://www.xquartz.org/
	    ⁃	webgl is required for viewing 3D plots in html files. webgl is available at: https://get.webgl.org/
	3.	People involved with data collection, processing, analysis, and/or submission:
	    ⁃	Lynch, E.M., Regalla, C., Harrichhausen, N., and Leonard, L.J. participated in total station surveying, assisted by Johanna Fischi (Boston University) and Tatum MacLeod (Western Washington University).
	    ⁃	Lynch, E.M. wrote the R scripts and markdown file.
