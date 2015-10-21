#Takes two polygon layers, intersects them, computes the weighted average of the intersecting area (where weighting is based on area proportion)
arealInterpolation<- function(voronPoly,watershedBdy){
  source("voronoipolygons.R",echo=T)
  
  
  pi <- intersect(voronPoly, watershedBdy)
    
  totalArea = gArea(watershedBdy)
  aveValues = numeric(0)
  weightVector = numeric(0)
  polyValues = pi@data['z']

  polyValues = as.vector(as.matrix(polyValues))
  
  #lazy for loop
  for(k in 1:length(pi@polygons)) {
    smallArea = gArea(SpatialPolygons(pi@polygons[k]))
    weight= smallArea/totalArea
    #weight = 0.3
    weightVector = c(weightVector,weight)
    #print(paste0("poly Weight: ", k , " ",weight))
    contrib = weight * pi@data[k,'z']
    #print(paste0("poly: ", k , " ",pi@data[k,'z']))
    aveValues = c(aveValues,contrib)
  }

  
  print("polyValues: ")
  print(polyValues)
  
  print("weights: ")
  print(weightVector)
  
  weightedMean = weighted.mean(polyValues,weightVector)
  print("weighted Mean: ")
  print(weightedMean)
  
  areaProp = sum(aveValues)
  
  return(areaProp)
}
