#' load WU and join to create a DF
#'
#' Loads text files scraped from the Weather Undeground and wrangles them into a spatial data frame.
#' @filesDir location of the files scraped from Weather Underground
#' @pwsLocations shapefile of point locations of weather stations. 
#' @datesList list of datetimes for the found in the list of files
#' @keywords PWS, Weather Underground munging
#' @export
#' @examples
#' loadAndJoinWU()

#take a directory of CSVs and data frame of observations joined with coordinates
loadAndJoinWU <- function(filesDir,pwsLocations,datesList) {
  #create a df
  df= data.frame(t(rep(NA,length(datesList))))
  names(df) <- c(as.character(datesList))
  df <- df[-1,]


  #add a column for the pws id / name
  name = data.frame(0)
  name <- name[-1,]

  #for each point try and load a csv of observations from the safe dir above
  for(i in 1:nrow(pwsLocations)) {
    p <- pwsLocations[i,]
    pwsId = p$name

    fileNameToLoad = paste(filesDir,pwsId,"_hourly.csv",sep="")
    #load the table
    if (file.exists(fileNameToLoad ) == TRUE) {
      myTable = read.csv(fileNameToLoad)


      myTableAsRow = t(myTable[,2]) #transpose the precip column to row
      myTableDates = t(myTable[,1]) # transpsose the dates of the to row
      colnaDf = colnames(df)
      colnames(myTableAsRow) = colnaDf #copy the field names so they match

      name = rbind(name,as.character(pwsId))

      df = rbind(df,myTableAsRow) #copy the row into the dataframe
     }
  }

  #bind the loaded data with an id field
  df2= cbind(df,name)

  #joint the obs data with geo file so we have coords for the interpolation
  merged = merge(df2,pwsLocations,by = "name",all = FALSE)

  return(merged)

} # end of load function
