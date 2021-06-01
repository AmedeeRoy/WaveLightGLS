getSSTPoint <- function(path, coord, time){
  # Download the base map
  data.sst.sat <- ncdf4::nc_open(path)
  
  ## LECTURE TEMPERATURE
  lat <- ncdf4::ncvar_get(data.sst.sat, varid="lat")
  lon <- ncdf4::ncvar_get(data.sst.sat, varid="lon")
  t <-ncdf4::ncvar_get(data.sst.sat, varid="time")
  t <- as.POSIXct("1981-01-01- 00:00:00", tz = "GMT") + seconds(t)
  
  sst <- sapply(1:nrow(coord), function(i){
    idx_x = which.min(abs(lon-coord[i,1]))
    idx_y = which.min(abs(lat-coord[i,2]))
    
    z_date <- as.Date(time[i])
    z = as.numeric(z_date-as.Date("1950-01-01"), units = "hours")
    idx_t = which.min(abs(t - time[i]))
    
    if(length(idx_y) == 0){
      s = NA
    } else {
      if(idx_t>1 & idx_t<length(t)& idx_x>1 & idx_x<length(lon) & idx_y>1 & idx_y<length(lat)){
        s <- ncdf4::ncvar_get(data.sst.sat, varid = "analysed_sst", start = c(idx_x, idx_y, idx_t),
                              count = c(1,1,1))
      } else {
        s = NA
      }
    }
    
    return(s)
  })
  
  nc_close( data.sst.sat )
  return(sst)
}

getSSTArea <- function(path, time){
  # Download the base map
  data.sst.sat <- ncdf4::nc_open(path)
  
  ## LECTURE TEMPERATURE
  lat <- ncdf4::ncvar_get(data.sst.sat, varid="lat")
  lon <- ncdf4::ncvar_get(data.sst.sat, varid="lon")
  t <-ncdf4::ncvar_get(data.sst.sat, varid="time")
  t <- as.POSIXct("1981-01-01- 00:00:00", tz = "GMT") + seconds(t)
  
  width = 2
  
  sst <- sapply(1:length(time), function(i){
    z_date <- as.Date(time[i])
    z = as.numeric(z_date-as.Date("1950-01-01"), units = "hours")
    idx_t = which.min(abs(t - time[i]))
    
    s <- max(ncdf4::ncvar_get(data.sst.sat, varid = "analysed_sst", start = c(111-width, 85-width, idx_t),
                              count = c(5,5,24)), na.rm = TRUE)
    # s <- mean(ncdf4::ncvar_get(data.sst.sat, varid = "analysed_sst", start = c(552, 425, idx_t),
    #                            count = c(1,1,4)), na.rm = TRUE)
    
    return(s)
  })
  
  nc_close( data.sst.sat )
  return(sst)
}
