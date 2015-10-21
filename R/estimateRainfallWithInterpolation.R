#' Rainfall estimation
#'
#' Takes a dataframe of locations and observed variables and computes a single estimates at each timestep for a given area i.e. the watershed.
#' @merged data frame of merged PWS data observation time-series
#' @datesList list of datetimes for the interpolation
#' @keywords PWS interpolation
#' @export
#' @examples
#' estimateRainfallWithInterpolation()

#This needs some restructuring to remove the watershed input (this should be a second analysis step)
estimateRainfallWithInterpolation<- function(merged,datesList,watershed) {
  #got to do this to make sure spatstat / gstat clash is avoided
  detach(package:gstat)
  library(gstat)

  library(plyr)
  library(automap)

  saveToTempFile = TRUE

  returnThiessen = TRUE #if false then use IDW estimate
  computeIDW = TRUE
  computeAutomapKrige = FALSE

  plotDensity = TRUE
  plotThiessen = FALSE

  #ditch this package so no confusion with gstat idw function
  #detach("package:spatstat",unload=TRUE)

  #create an id column to help joining dfs later
  idCol = 1:nrow(merged)
  #create df for holding cross-validation residuals
  residualsFinalDf = data.frame(matrix(ncol = length(datesList), nrow = nrow(merged)))
 # columnNamesDf = merged[1,1:length(residualsFinalDf)]
#  columnNamesDf =
 # colnames(residualsFinalDf) = colnames(columnNamesDf)
  colnames(residualsFinalDf) = datesList

cbind(name = merged$name, residualsFinalDf)
  #create a grid for the interp surface based on bbox of pws data
  #assuming projection set up to metres
  bbox = pwsLocations@bbox
  bboxEdge = 100 #give a margin
  x.range <- as.numeric(c(bbox[1,1]-bboxEdge, bbox[1,2]+bboxEdge))  # min/max longitude of the interpolation area
  y.range <- as.numeric(c(bbox[2,1]-bboxEdge, bbox[2,2]+bboxEdge))  # min/max latitude of the interpolation area

  cell = 1000 # metre cell size
  grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = cell), y = seq(from = y.range[1],  to = y.range[2], by = cell))  # expand points to grid
  coordinates(grd) <- ~x + y

  #x11()

  gridded(grd) <- TRUE

  #now declare coordinate fields
  merged$x <- merged$coords.x1
  merged$y <- merged$coords.x2
  coordinates(merged) <- c("x", "y")



  merged@data$name <- NULL
  merged@data$coords.x1 <- NULL
  merged@data$coords.x2 <- NULL

  #colnames(merged)[colnames(merged)=="x2013-01-01"] <- "new_name"
  #dw <- idw(X2013.01.30 ~ 1, merged, grd)  # apply idw model for the data
  #idw <- idw(X20130130 ~ 1, merged, grd)  # apply idw model for the data
  #idw <- idw(formula = X20130101 ~ 1, merged, grd)  # apply idw model for the data

  precipValues = numeric(0) #vector to hold area stats
  dateValues = character(0)


  #loop over the dates, predict the surface, compute the area stats
  for (j in 1:length(datesList)) {
  #for (j in 1:10) {

    print(paste0("Iteration: ",j ," of ",length(datesList)))

    areaPropIDW = NULL
    areaPropThiessen = NULL
    areaPropKrige = NA

    naVals = !is.na(merged@data[which(as.POSIXlt(names(merged@data),tz = "UTC") == datesList[j])]) #clunky logical vector of na values
    naValsDf = data.frame(cbind(naVals, 1:nrow(naVals)))

    numb=  sum( naVals == TRUE )
    if (numb == 0){
      areaPropThiessen = NA
      areaPropIDW = NA
      areaPropKrige = NA

      #cleaned = merged[,j]
    } else{

      idCol = 1:nrow(merged)
      mergedCopy = merged
      mergedCopy@data = cbind(merged@data,idCol)

      cleaned = mergedCopy[naVals[,1],] # select out the non-na vals

      listOfColNames = names(cleaned@data)
      listOfColNames= listOfColNames[1:length(listOfColNames)-1]
      colna = as.POSIXlt(listOfColNames,tz = "UTC")
      colToRename = which(colna==datesList[j])
      names(cleaned@data)[colToRename] <- "z"
      if (computeIDW == TRUE) {
        #rename function vector as z, seems like a hack - cant seem to pass in dependent variable (field name) in otherwise

        #do the interpolation
        idw <- idw( z ~ 1, cleaned, grd)
        #check we have proj
        idw@proj4string =  pwsLocations@proj4string

        #compute statistics
        r <- raster(idw)
        # extract mean values for polygons
        areaPropIDW <- extract(r, watershed, fun=mean)



        if (sum(cleaned@data$z) == 0) {
          print("Not IDW cross-val as all values are 0")

          rowNamesNotNa= rownames(cleaned@data)

          idColDf = data.frame(idCol)
          zeroResidualsDf = data.frame(cbind(as.numeric(rowNamesNotNa), cleaned@data$z))
          colnames(zeroResidualsDf) = c("idCol","residual")
          idColDf = data.frame(idCol)

          #join the non-na values to the right ids for going back in the data frame
          residualsForFinalDf= join(idColDf,zeroResidualsDf,by="idCol")

          #insert this residuals liit into the output residuals frame
          #residualsFinalDf[which(colna==datesList[j])] = residualsForFinalDf$residual
          residualsFinalDf[j] = residualsForFinalDf$residual
        } else {
          #use the gstat idw.cv to get some residuals
          idw.cv = krige.cv(z~1, cleaned,  nfold=100)
          rowNamesNotNa= rownames(cleaned@data)
          idColDf = data.frame(idCol)
          rowNames=rownames(idw.cv@data)
          residualsSetDf = data.frame(cbind(as.numeric(rowNames), idw.cv@data$residual))

          colnames(residualsSetDf) = c("idCol","residual")

          #join the non-na values to the right ids for going back in the data frame
          residualsForFinalDf = join(idColDf,residualsSetDf,by="idCol")

          #insert this residuals liit into the output residuals frame
          #residualsFinalDf[which(colna==datesList[j])] = residualsForFinalDf$residual
          residualsFinalDf[j] = residualsForFinalDf$residual
        }
      }

      if (computeAutomapKrige == TRUE) {
        #use automap package to estimate variogram
        if (sum(cleaned@data$z) == 0) {
          print("Not kriging as all values are 0")

          areaPropKrige = 0
          rowNamesNotNa= rownames(cleaned@data)

          idColDf = data.frame(idCol)
          zeroResidualsDf = data.frame(cbind(as.numeric(rowNamesNotNa), cleaned@data$z))
          colnames(zeroResidualsDf) = c("idCol","residual")
          idColDf = data.frame(idCol)

          #join the non-na values to the right ids for going back in the data frame
          residualsForFinalDf= join(idColDf,zeroResidualsDf,by="idCol")

          #insert this residuals liit into the output residuals frame
          #residualsFinalDf[which(colna==datesList[j])] = residualsForFinalDf$residual
          residualsFinalDf[j] = residualsForFinalDf$residual

        }else {
          print(paste0("Starting kriging interp..."))
          kriging_result = autoKrige(z~1, cleaned, grd)

          plot(kriging_result)


          #Do cross-validation using the interpolation routine
          nFolds = nrow(cleaned@data)
          #nFolds = 10
          print(paste0("Using ", nFolds, " folds"))

          #do the cross-val krige
          kriging_result_cv = autoKrige.cv(z~1, cleaned, nfold = nFolds)

          print(paste0("print krige CV summary: "))
                summary(kriging_result_cv)


          #get kriging result df table
          krigResultDf = kriging_result_cv$krige.cv_output@data
          #get the rownames for joining
          rowNames= rownames(kriging_result_cv$krige.cv_output@data)
          #get the residuals
          residualsSetDf = data.frame(cbind(as.numeric(rowNames), krigResultDf$residual))
          colnames(residualsSetDf) = c("idCol","residual")
          idColDf = data.frame(idCol)

          #join the non-na values to the right ids for going back in the data frame
          residualsForFinalDf = join(idColDf,residualsSetDf,by="idCol")

          #insert this residuals liit into the output residuals frame
          #residualsFinalDf[which(colna==datesList[j])] = residualsForFinalDf$residual
          residualsFinalDf[j] = residualsForFinalDf$residual

          krigeR = raster(kriging_result$krige_output)
          areaPropKrige <- extract(krigeR, watershed, fun=mean)

          Sys.sleep(0)
        }
      }

      cleanedDf = cleaned@data
      cleanedDf2 = cbind(cleaned$x,cleaned$y, cleanedDf[colToRename])
      newPoints = SpatialPointsDataFrame(cbind(cleaned$x,cleaned$y),cleanedDf2)
      newPoints = SpatialPointsDataFrame(cbind(cleaned$x,cleaned$y),cleanedDf[colToRename],proj4string = CRS("+init=epsg:27700"))

      #check there is enough obs for thiessen estimation
      if (nrow(newPoints@data)==1) {
        areaPropThiessen =  as.numeric(newPoints@data[1])
      } else{
        voronPoly = voronoipolygons(newPoints)
        voronPoly@data <- cleanedDf2
        areaPropThiessen = arealInterpolation(voronPoly,watershed)
      }

      print(paste0("Thiesssen precip ", areaPropThiessen))
      print(paste0("IDW precip ", areaPropIDW))
      print(paste0("Automap Krige precip ", areaPropKrige))

      #station coords for plotting
      cleanedCoords  = data.frame(cbind(cleaned$x,cleaned$y))

      idw.output=as.data.frame(idw)
      #names(idw.output)[1:3]<-c("x","y","var1.pred")

      if (plotDensity == TRUE) {

        Sys.sleep(0)
        print(paste0("Precip Surface ",j))

        print(ggplot() + geom_tile(data = idw.output, aes(x = x, y = y, fill = var1.pred)) +
                geom_path(data = watershed, aes(long, lat, group = group), colour = "grey"))

        Sys.sleep(0)
      }

      #
      #     voronoi <- deldir(merged$x, merged$y)
      #
      #     #Now we can make a plot
      #     ggplot(data=merged@data, aes(x=long,y=lat)) +
      #       #Plot the voronoi lines
      #       geom_segment(
      #         aes(x = x1, y = y1, xend = x2, yend = y2),
      #         size = 2,
      #         data = voronoi$dirsgs,
      #         linetype = 1,
      #         color= "#FFB958") +
      #       #Plot the points
      #       geom_point(
      #         fill=rgb(70,130,180,255,maxColorValue=255),
      #         pch=21,
      #         size = 4,
      #         color="#333333")
      #

      if (plotThiessen == TRUE) {
        print( ggplot()+
                 geom_path(data = watershed, aes(long, lat, group = group), colour = "grey") +
                 geom_path(data = voronPoly, aes(long, lat, group = group), colour = "blue") +
                 geom_point(data=cleanedCoords, mapping=aes(x=X1, y=X2)) +
                 coord_fixed())
        Sys.sleep(0)
      }

    }#end of if to check there is some values to interp or nearest neighbour

    # assign these to a new variable
    precipValues = rbind(precipValues,t(c(areaPropThiessen,areaPropIDW,areaPropKrige)))
    dateValues = rbind(dateValues,datesList[j])

    #save results to temp file
    if (saveToTempFile == TRUE){
      dfPrecipTemp = data.frame(cbind(dateValues,precipValues))
      write.table(dfPrecipTemp, file = "temp.csv", sep = ",")
    }

  } #end for



  dfPrecip = data.frame(cbind(dateValues,precipValues))
  names(dfPrecip) = c("datetime","thiessen","idw","autoKrige")


  #return(dfPrecip)
  MyList<- list("a"=dfPrecip, "b"=residualsFinalDf)
  return(MyList)
  #return(dfPrecip,residualsFinalDf)
}#end of interp rainfall
