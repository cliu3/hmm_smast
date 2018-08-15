#install.packages("gdata")
#install.packages("R.matlab")
#install.packages("zoo")
library(gdata)
library(R.matlab)
library(zoo)
source("matlab_time.R")
source("datenum.R")

# raw tag directory
tagdir <- "raw_tags"
outdir <- "processed_tags"

df <- read.xls("Halibut_inventory.xlsx",stringsAsFactors=FALSE, fileEncoding="latin1")
#df[] <- lapply(df, as.character)

ptags <- c(6);

for(i in ptags)
  #i <- 1
{
  i = which(df[,"Fish.number"]==i)
  tag <- list()
  
  
  
  source("metadata.R")
  
  
  ####
  tag[["fish_id"]]   <- df[i,"Fish.number"]
  tag[["tag_id"]]    <- toString(df[i,"PTT.ID"])
  tag[["type"]]      <- df[i,"Tag.Type"]
  tag[["length"]]    <- df[i,"fork.length.flat..cm."]
  #tag[["sex"]]       <- df[i,"SEX"]
  #tag[["maturity"]]  <- df[i,"MATURITY"]
  tag[["release_dnum"]] <- datenum(format(as.POSIXct(paste(df[i,"release.date."],df[i,"Release.time."]), tz="America/New_York"),tz="UTC"))
  tag[["release_lon"]] <- as.numeric(df[i,"Longitude..decimal.D."])
  tag[["release_lat"]] <- as.numeric(df[i,"Latitude..decimal.D."])
  tag[["recapture_lon"]] <- as.numeric(df[i,"Argos.Location.Longitude"])
  tag[["recapture_lat"]] <- as.numeric(df[i,"Argos.Location.Latitude"])
  tag[["recapture_dnum"]] <- datenum(df[i,"Earliest.Argos.Location.time..UTC."])
  tag[["recap_uncertainty_km"]] <- as.numeric(df[i,"Error.radius..m."])
  
  tag[["datafile"]]  <- file.path(tagdir, paste(tag[["tag_id"]],'/', tag[["tag_id"]], '-Series.csv', sep=""))
  
  #set time-dependent fields
  if(tag[["type"]]=="DST"){
    source("process_staroddi.R")
    tag = process_staroddi(tag)
  } else if(tag[["type"]]=="PSAT"){
    if(i==14) {
      source("process_wildlife_archived.R")
      tag = process_wildlife_archived(tag)
    }else {
    source("process_wildlife.R")
    tag = process_wildlife(tag)
  } } else {
    print(paste("type is", tag[["type"]]))
    print("not setup to read the type")
  }
  
  matfilename <- file.path(outdir, paste(tag$fish_id,"_raw.mat",sep=""))
  writeMat(matfilename, tag=tag)
  
}