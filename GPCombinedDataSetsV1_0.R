#####################################################
#
# Functions to combine a given set of data
#
# Used by for some integrated data sets where
# correlation is calculated over a combined set of spatial data
# 
# 
#       generateXCoordsGivenGenData   - retrieve the coordinates from generated data sets  
# Functions to generate combined set of spatial distance data
#       generateCombinedDis           - given two sets of coords
#       generateCombinedDisGenData    - given two set of simulated spatial data
#       generateCombinedDisReal       - given two set of Real data assuming contain XYCoords
#       generateCombinedDisMixed      - given two mixed sets of data: 
#                                       1 set of presence-only and 1 set of abundance  data
#       generateCombinedDisAny        - given two sets of data in a range of data format

################################
# generateXCoordsGivenGenData
#
# returns xcoords given spatially generated data
# data returned in format xcoords (lat, long)
# 
# spg_gen_data    - set of spatially generated data
# 
generateXCoordsGivenGenData <- function(spg_gen_data) {
  lat = (1:spg_data$grid_unit_y)*spg_data$scale_y
  long = (1:spg_data$grid_unit_x)*spg_data$scale_x
  return(xcoords = cbind(lat,long))
}


################################
# generateCombinedDis
#
# Combine distances from two spatial datasets given xycoords for both
# returns spatial distance matrix
#
# coords1      - xycoords for dataset 1
# coords2      - xycoords for dataset 2
#
generateCombinedDis <- function(coords1, coords2){
  xcoords <- as.numeric(rbind(coords1,coords2) )
  dim_1 = dim(coords1)[1] + dim(coords2)[1]
  dim_2 = dim(coords1)[2]
  dim(xcoords) <- c(dim_1,dim_2)
  dis <- (fields::rdist(xcoords, xcoords))
  return(list(dis=dis, comb_coords = xcoords))
}


################################
# generateCombinedDisGenData
#
# Combine distances from two spatial datasets
# for simulated data sets
# returns spatial distance matrix
#
# spg_gen_data1     - dataset 1 in generated data format
# spg_gen_data1     - dataset 2 in generated data format
#
generateCombinedDisGenData <- function(spg_gen_data1, spg_gen_data2){
  xcoords1 <- generateXCoordsGivenGenData(spg_gen_data1)
  xcoords2 <- generateXCoordsGivenGenData(spg_gen_data2)
  return(generateCombinedDis(xcoords1,xcoords2))
}


################################
# generateCombinedDisGenData
#
# Combine distances from two spatial datasets
# for real data sets
# returns spatial distance matrix
#
# real_data1        - real dataset 1 with xy coordinates fields coordXY
# spg_gen_data1     - real dataset 2 with xy coordinates fields coordXY
#
generateCombinedDisReal <- function(real_data1, real_data2) {
  xcoords1 <- real_data1$coordXY
  xcoords2 <- real_data2$coordXY
  return(generateCombinedDis(xcoords1,xcoords2))
}


################################
# generateCombinedDisGenData
#
# Combine distances from two spatial datasets
# for mixed presence-only and abundance data sets
# returns spatial distance matrix
#
# mix_data1   - presence-only data (mgp_data_PO)
# mix_data2   - abundance data (mgp_data_Abund)
#
generateCombinedDisMixed <- function(mix_data1, mix_data2) {
  if ( is.null(mix_data1) || is.null(mix_data2)) {cat("generateCombinedDisMixed ","Must  have non-null data","\n");return(NULL)}
  if ( !is.null(mix_data1$coordXY) ) {
    xcoords1 <- mix_data1$coordXY
  } else if (!(is.null(mix_data1$grid_unit_y)) ) {
    xcoords1 <- generateXCoordsGivenGenData(mix_data1)
  } else if (is.matrix(mix_data1) && dim(mix_data1)[2]>=2 ) {
    xcoords1 <- mix_data1
  } else { cat("generateCombinedDisMixed: ","data 1 type unknown error"); return(NULL)}
  if ( !is.null(mix_data2$coordXY) ) {
    xcoords2 <- mix_data2$coordXY
  } else if (!(is.null(mix_data2$grid_unit_y)) ) {
    xcoords2 <- generateXCoordsGivenGenData(mix_data2)
  } else if (is.matrix(mix_data2) && dim(mix_data2)[2]>=2 ) {
    xcoords2 <- mix_data2
  } else { cat("generateCombinedDisMixed: ","data 2 type unknown error"); return(NULL)}
  return(generateCombinedDis(xcoords1,xcoords2))
}



################################
# generateCombinedDisAny
#
# generateds combined distances for large number of different
# coordinate formats
# returns spatial distance matrix
#
# my_data1     - dataset 1 in real or generated data format
# my_data2     - dataset 2 in real or generated data format
#
generateCombinedDisAny <- function(my_data1, my_data2) {
  if ( is.null(my_data1) || is.null(my_data2)) {cat("generateCombinedCoordsAny ","Must  have non-null data","\n");return(NULL)}
  if ( !is.null(my_data1$coordXY)  && !is.null(my_data2$coordXY)) {
    return(generateCombinedDisReal(my_data1, my_data2))
  } else if ( !(is.null(my_data1$grid_unit_y)) && !(is.null(my_data2$grid_unit_y))  ) {
    return(generateCombinedDisGenData(my_data1, my_data2))
  } else if ( is.matrix(my_data1)  && is.matrix(my_data2) && dim(my_data1)[2]>=2 && dim(my_data2)[2]>=2) {
    # assume data is in coord form with first col lat and second col long
    # any other fields ignored
    return(generateCombinedDis(my_data1,my_data1))
  }
  return(generateCombinedDisMixed(my_data1, my_data2))
}



