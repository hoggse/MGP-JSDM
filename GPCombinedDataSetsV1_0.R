################################
# generateXCoordsGivenGenData
#
# returns xcoords (lat, long) from spatial data
generateXCoordsGivenGenData <- function(spg_gen_data) {
  lat = (1:spg_data$grid_unit_y)*spg_data$scale_y
  long = (1:spg_data$grid_unit_x)*spg_data$scale_x
  return(xcoords = cbind(lat,long))
}


################################
# generateCombinedDis
#
# Combine distances from two spatial datasets
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
# mix_data1   - presence-only data (mgp_data_PO)
# mix_data2   - abundance data (mgp_data_Abund)
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