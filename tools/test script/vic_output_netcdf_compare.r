
library(ncdf4)
library(ncdf4.helpers)

ERROR <- 1
SUCCESS <- 0

assert <- function(expression, description) {
    if ( identical(expression, FALSE) ) {
        message(description)
        quit(status=ERROR)
    }
}

checkValues <- function(value1, value2) {
    if (is.na(value1) && value2 == 0) {
        return(TRUE)
    } else if (value1 == 0 && is.na(value2)) {
        return(TRUE)
    } else {
        return(value1 == value2)
    }
}

compare.variable <- function(nc1, nc2, var, var2, start=NA, count=NA) {
          data1 <- ncvar_get(nc1, var$name, start, count)
          data2 <- ncvar_get(nc2, var2$name, start, count)
          data1[is.na(data1)] <- 0
          data2[is.na(data2)] <- 0
  
          #write(data1, file=paste("internalData/internal", var$name, "Data1", sep=""))
          #write(data2, file=paste("internalData/internal", var$name, "Data2", sep=""))
   
          assert(var$nvars == var2$nvars, paste("The data must have the same size", var$nvars, var2$nvars))
          assert(var$ndims == var2$ndims, paste("The data must have the same number of dimensions", var$ndims, var2$ndims))
          assert(length(data1) == length(data2), paste("The data must have the same size", length(data1), length(data2)))
  
          diffs <- abs(data1 - data2)
          #diffs[is.na(diffs)] <- 0
          s <- sum(diffs)
          differentEntries <- diffs[diffs > 0.001]
#          message("dimensions: ", str(dim(diffs)))
#          message("maximum difference: ", str(max(diffs)))
#          message("accumulated differences: ", s)
#          message("number of strictly different entries: ", length(diffs[diffs!=0]))
#          message("number of unacceptable entries (diff > 0.001) : ", length(differentEntries))

          assert(length(differentEntries) == 0, paste("Difference indices ", str(differentEntries), "diffs: ", str(diffs[differentEntries]),"data1:",str(data1[differentEntries]), "data2:",str(data2[differentEntries]) ))
}

compare.netCDF <- function(file1Name = NA, file2Name = NA) {
  # args comes from invoking th
  args <- commandArgs(TRUE)

  if (is.na(file1Name)) {
    file1Name <- args[1]
  }  
  if (is.na(file2Name)) {
    file2Name <- args[2]
  }
     
  assert(is.na(file1Name) == FALSE && is.na(file2Name) == FALSE, "Must supply two netcdf files to compare either as function parameters, or Rscript arguments")
  
  #file1Name <- "../out/automated_4.1.2_netcdf_2013_10_30__14_03_57/results.nc"
  #file2Name <- "outputs/vic5fluxes.nc"
  
  message("comparing files: ", file1Name, " and ", file2Name)
  
  nc1 <- nc_open(file1Name)
  nc2 <- nc_open(file2Name)
  
  assert(nc1$nvars == nc2$nvars, "Files must contain the same number of variables")
  
  #for each variable
  for( i in 1:nc1$nvars) {
      var <- nc1$var[[i]]
      var2 <- nc2$var[[i]]
      message("checking variable: ", var$name, " and ", var2$name)

      if (var$ndims == 4) {
        for (d in 1:var$varsize[[1]]) {
          start <- c(d,1,1,1)
          count <- c(1,-1,-1,-1)
          compare.variable(nc1, nc2, var, var2, start, count) 
        }
      } else {
        compare.variable(nc1, nc2, var, var2)
      }
  }
  
  nc_close(nc1)
  nc_close(nc2)
  
  message("The data in these files is identical")
  quit(status=SUCCESS)
  
}

compare.netCDF()

