best <- read.csv("https://salmonetics.github.io/data-Spawning/estimates/combined_abundance_geomeans_3-04-2019.csv",stringsAsFactors=F)
pops <- read.csv("https://salmonetics.github.io/CAX/data/ca-data-all 02-08-2019 18 55_Populations.csv",stringsAsFactors=F)
pops <- select(pops,NMFS_POPID,POPID,ESU_DPS,MAJORPOPGROUP)
best <- merge(pops,best)

best <- filter(best, ESU_DPS != "Salmon, Chinook (Lower Columbia River ESU)")
best <- filter(best, ESU_DPS != "Salmon, Chinook (Upper Willamette River ESU)")
best <- filter(best, ESU_DPS != "Salmon, chum (Columbia River ESU)")
best <- filter(best, ESU_DPS != "Salmon, chum (Columbia River ESU)")
best <- filter(best, ESU_DPS != "Salmon, coho (Lower Columbia River ESU)")
best <- filter(best, ESU_DPS != "Salmon, chum (Columbia River ESU)")
best <- filter(best, ESU_DPS != "Salmon, coho (Oregon Coast ESU)")
best <- filter(best, ESU_DPS != "Steelhead (Lower Columbia River DPS)")
best <- filter(best, ESU_DPS != "Steelhead (Upper Willamette River DPS)")
best <- filter(best, ESU_DPS != "")
best <- filter(best, SPAWNINGYEAR >= 1980)

best$NOSAIJ <- ifelse(is.na(best$NOSAIJ),0,best$NOSAIJ)
best$NOSAEJ <- ifelse(is.na(best$NOSAEJ),0,best$NOSAEJ)
best$TSAIJ <- ifelse(is.na(best$TSAIJ),0,best$TSAIJ)
best$TSAEJ <- ifelse(is.na(best$TSAEJ),0,best$TSAEJ)
best$count <- ifelse(best$NOSAIJ>best$NOSAEJ,best$NOSAIJ,best$NOSAEJ)
best$tcount <- ifelse(best$TSAIJ>best$TSAEJ,best$TSAIJ,best$TSAEJ)
best$hcount <- round(best$tcount,0) - best$count

best2 <- summarise(group_by(best,ESU_DPS,MAJORPOPGROUP, SPAWNINGYEAR),COUNT=sum(count),HCOUNT=sum(hcount))
best2 <- filter(best2,!(ESU_DPS =="Steelhead (Middle Columbia River DPS)" && MAJORPOPGROUP != "Yakima"))
write.csv(best2,"CAX_compdata_03-08-2019.csv",row.names=FALSE)
########################
#get dam counts
getSRSpringCounts <- function(proj){
  url1 <- 'http://www.cbr.washington.edu/dart/cs/php/rpt/adult_annual.php?sc=1&outputFormat=csv&proj='
  url2 <- '&startdate=3%2F1&enddate=8%2F17&run='
  url <- paste(url1,proj,sep="")
  url <- paste(url,url2,sep="")
  temp <- read.csv(url)
  temp <- filter(temp,!is.na(Year))
  return(temp)
}

getUCSpringCounts <- function(proj){
  url1 <- 'http://www.cbr.washington.edu/dart/cs/php/rpt/adult_annual.php?sc=1&outputFormat=csv&proj='
  url2 <- '&startdate=4%2F14&enddate=8%2F18&run='
  url <- paste(url1,proj,sep="")
  url <- paste(url,url2,sep="")
  temp <- read.csv(url)
  temp <- filter(temp,!is.na(Year))
  return(temp)
}

getSRFallCounts <- function(proj){
  url1 <- 'http://www.cbr.washington.edu/dart/cs/php/rpt/adult_annual.php?sc=1&outputFormat=csv&proj='
  url2 <- '&startdate=8%2F18&enddate=12%2F15&run='
  url <- paste(url1,proj,sep="")
  url <- paste(url,url2,sep="")
  temp <- read.csv(url)
  temp <- filter(temp,!is.na(Year))
  return(temp)
}

getAnnualDamCounts <- function(proj){
  url1 <- 'http://www.cbr.washington.edu/dart/cs/php/rpt/adult_annual.php?sc=1&outputFormat=csv&proj='
  url2 <- '&startdate=1%2F1&enddate=12%2F31&run='
  url <- paste(url1,proj,sep="")
  url <- paste(url,url2,sep="")
  temp <- read.csv(url)
  temp <- filter(temp,!is.na(Year))
  return(temp)
}



SRSpring <- getSRSpringCounts("LWG")
SRSpring <- filter(data.frame(group="srchinook", year=SRSpring$Year, count=SRSpring$Chinook),count>0)
SRFall <- getSRFallCounts("LWG")
SRFall <- filter(data.frame(group="fallchinook", year=SRFall$Year, count=SRFall$Chinook),count>0)
SRSthd <- getAnnualDamCounts("LWG")
SRSthd <- filter(data.frame(group="srsteelhead", year=SRSthd$Year, count=SRSthd$Wild.Steelhead) ,count>0) 
UCSthd <- getAnnualDamCounts("PRD")
UCSthd <- filter(data.frame(group="ucsteelhead", year=UCSthd$Year, count=UCSthd$Wild.Steelhead),count>0)
UCSpring <- getUCSpringCounts("PRD")
UCSpring <- filter(data.frame(group="ucchinook", year=UCSpring$Year, count=UCSpring$Chinook),count>0)
Sock <- getAnnualDamCounts("LWG")
Sock <- filter(data.frame(group="sockeye", year=Sock$Year, count=Sock$Sockeye),count>0)
Yak <- getAnnualDamCounts("PRO")
Yak <- filter(data.frame(group="yakima", year=Yak$Year, count=Yak$Wild.Steelhead),count>0)

data1 <- rbind(SRSpring,SRFall,SRSthd,UCSthd,UCSpring,Sock,Yak)
data1 <- filter(data1,year>=1980,year<=2019)

write.csv(data1,"TAC_Dam_Counts-03-08-2018.csv",row.names=FALSE)


##########################################
#get dam flow data

getMeanDamData <- function(year){
  df <- getMeanDamDataFrame(year)
  lst <- list()
  lst$srchinook <- getSRChinookDamData(df)
  lst$fallchinook <- getFallChinookDamData(df)
  lst$ucchinook <- getUCChinookDamData(df)
  lst$srsteelhead <- getSRSteelheadDamData(df)
  lst$ucsteelhead <- getUCSteelheadDamData(df)
  lst$sockeye <- getSRSockeyeDamData(df)
  return(lst)
}

getMeanDamDataFrame <- function(year){
  url1 <- "http://www.cbr.washington.edu/dart/cs/php/rpt/mg.php?sc=1&mgconfig=adult&outputFormat=csvSingle&year%5B%5D=" 
  url2 <- "&loc%5B%5D=LWG&loc%5B%5D=PRD&ftype%5B%5D=fc&ftype%5B%5D=fb&ftype%5B%5D=fs&data%5B%5D=&data%5B%5D=Outflow&data%5B%5D=Spill&data%5B%5D=Temp+%28WQM%29&startdate=1%2F1&enddate=12%2F31&avgyear=0&sumAttribute=none&consolidate=1&zeros=1&grid=1&y1min=0&y1max=&y2min=&y2max=&size=medium"
  url <- paste(url1,year,url2,sep="")
  df <- read.csv(url,stringsAsFactors=FALSE)
  return(df)
}

mergeDamData <- function(df.counts,df.flow,df.spill,df.temp){
  df.counts <- data.frame(day=df.counts$mm.dd,counts=df.counts$value)
  df.flow <- data.frame(day=df.flow$mm.dd,flow=df.flow$value)
  df.spill <- data.frame(day=df.spill$mm.dd,spill=df.spill$value)
  df.temp <- data.frame(day=df.temp$mm.dd,temp=df.temp$value)
  df <- merge(df.counts,merge(df.flow,merge(df.spill,df.temp,all.x=TRUE),all.x=TRUE),all.x=TRUE)
}

getSRChinookDamData <- function(df){
  df.counts <- filter(df,location=="LWG", datatype=="Adult Passage", parameter=="Chin")
  df.flow <- filter(df,location=="LWG", datatype=="Outflow")
  df.spill <- filter(df,location=="LWG", datatype=="Spill")
  df.temp <- filter(df,location=="LWG", datatype=="Temp (WQM)")
  df <- mergeDamData(df.counts,df.flow,df.spill,df.temp)
  df <- df[60:229,]
  return(getMeanValues(df))
}

getFallChinookDamData <- function(df){
  df.counts <- filter(df,location=="LWG", datatype=="Adult Passage", parameter=="Chin")
  df.flow <- filter(df,location=="LWG", datatype=="Outflow")
  df.spill <- filter(df,location=="LWG", datatype=="Spill")
  df.temp <- filter(df,location=="LWG", datatype=="Temp (WQM)")
  df <- mergeDamData(df.counts,df.flow,df.spill,df.temp)
  df <- df[230:365,]
  return(getMeanValues(df))
}  

getUCChinookDamData <- function(df){
  df.counts <- filter(df,location=="PRD", datatype=="Adult Passage", parameter=="Chin")
  df.flow <- filter(df,location=="PRD", datatype=="Outflow")
  df.spill <- filter(df,location=="PRD", datatype=="Spill")
  df.temp <- filter(df,location=="PRD", datatype=="Temp (WQM)")
  df <- mergeDamData(df.counts,df.flow,df.spill,df.temp)
  df <- df[60:225,]
  return(getMeanValues(df))
}  

getSRSteelheadDamData <- function(df){
  df.counts <- filter(df,location=="LWG", datatype=="Adult Passage", parameter=="Stlhd")
  df.flow <- filter(df,location=="LWG", datatype=="Outflow")
  df.spill <- filter(df,location=="LWG", datatype=="Spill")
  df.temp <- filter(df,location=="LWG", datatype=="Temp (WQM)")
  df <- mergeDamData(df.counts,df.flow,df.spill,df.temp)
  df <- df[1:365,]
  return(getMeanValues(df))
}  

getUCSteelheadDamData <- function(df){
  df.counts <- filter(df,location=="PRD", datatype=="Adult Passage", parameter=="Stlhd")
  df.flow <- filter(df,location=="PRD", datatype=="Outflow")
  df.spill <- filter(df,location=="PRD", datatype=="Spill")
  df.temp <- filter(df,location=="PRD", datatype=="Temp (WQM)")
  df <- mergeDamData(df.counts,df.flow,df.spill,df.temp)
  df <- df[1:365,]
  return(getMeanValues(df))
}

getSRSockeyeDamData <- function(df){
  df.counts <- filter(df,location=="LWG", datatype=="Adult Passage", parameter=="Sock")
  df.flow <- filter(df,location=="LWG", datatype=="Outflow")
  df.spill <- filter(df,location=="LWG", datatype=="Spill")
  df.temp <- filter(df,location=="LWG", datatype=="Temp (WQM)")
  df <- mergeDamData(df.counts,df.flow,df.spill,df.temp)
  df <- df[1:365,]
  return(getMeanValues(df))
}  

getMeanValues <- function(df){ 
  df$counts <- as.numeric(df$counts)
  df$flow <- as.numeric(df$flow)
  df$spill <- as.numeric(df$spill)
  df$temp <- as.numeric(df$temp)
  df <- filter(df,counts>0)
  sum.counts <- sum(df$counts)
  df$counts <- df$counts/sum.counts
  lst <- list()
  lst$flow <- getSum(select(df,counts,flow))
  lst$spill <- getSum(select(df,counts,spill))
  lst$temp <- getSum(select(df,counts,temp))
  return(lst)
}

getSum <- function(df){
  names(df) <- c("counts","value")
  df <- filter(df,value>0)
  return(sum(df$counts*df$value))
}

getAllMeans <- function(start,stop){
  options(warn=-1)
  df <- data.frame()
  for(year in start:stop){
    print(year)
    df1 <- data.frame()
    lst <- getMeanDamData(year)
    df1 <- rbind(df1,data.frame(year=year,group="srchinook",flow=lst$srchinook$flow,spill=lst$srchinook$spill,temp=lst$srchinook$temp))
    df1 <- rbind(df1,data.frame(year=year,group="fallchinook",flow=lst$fallchinook$flow,spill=lst$fallchinook$spill,temp=lst$fallchinook$temp))
    df1 <- rbind(df1,data.frame(year=year,group="srsteelhead",flow=lst$srsteelhead$flow,spill=lst$srsteelhead$spill,temp=lst$srsteelhead$temp))
    df1 <- rbind(df1,data.frame(year=year,group="ucchinook",flow=lst$ucchinook$flow,spill=lst$ucchinook$spill,temp=lst$ucchinook$temp))
    df1 <- rbind(df1,data.frame(year=year,group="ucsteelhead",flow=lst$ucsteelhead$flow,spill=lst$ucsteelhead$spill,temp=lst$ucsteelhead$temp))
    df1 <- rbind(df1,data.frame(year=year,group="sockeye",flow=lst$sockeye$flow,spill=lst$sockeye$spill,temp=lst$sockeye$temp))
    df <- rbind(df,df1)
  }
  options(warn=0)
  return(df)
}

allmeans <- getAllMeans(1980,2018)
write.csv(filter(allmeans,year<2019),"WeightedMeanDamData.csv",row.names=FALSE)

