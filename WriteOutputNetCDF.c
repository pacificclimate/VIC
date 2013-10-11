#include "WriteOutputNetCDF.h"
#include "netcdfcpp.h"

const char* WriteOutputNetCDF::getDescriptionOfOutputType() {
  return "NETCDF";
}

void WriteOutputNetCDF::openFile() {
  //TODO: all the things
}

void WriteOutputNetCDF::write_data(out_data_struct* out_data, const dmy_struct* dmy, int dt,
    const ProgramState* state) {

  // Loop over output files
  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {

    NcFile ncFile(dataFiles[file_idx]->filename, NcFile::Replace, NULL, 0, NcFile::Classic);

    if (!ncFile.is_valid()) {
      fprintf(stderr, "Error: could not open netCDF file for output, skipping... (%s)\n", dataFiles[file_idx]->filename);
      continue;
    }

    // Add global attributes.
    ncFile.add_att("institution", "Pacific Climate Impacts Consortium (PCIC), Victoria, BC, www.pacificclimate.org");


    NcDim* timeDim = ncFile.add_dim("time"); //time is unlimited.
    NcDim* latDim = ncFile.add_dim("lat");
    NcDim* lonDim = ncFile.add_dim("lon");


#if !OUTPUT_FORCE
    // Write the date
    if (dt < 24) {
      // Write year, month, day, and hour
      fprintf(dataFiles[file_idx]->fh, "%04i\t%02i\t%02i\t%02i\t",
              dmy->year, dmy->month, dmy->day, dmy->hour);
    }
    else {
      // Only write year, month, and day
      fprintf(dataFiles[file_idx]->fh, "%04i\t%02i\t%02i\t",
              dmy->year, dmy->month, dmy->day);
    }
#endif
    // Loop over this output file's data variables
    for (int var_idx = 0; var_idx < dataFiles[file_idx]->nvars; var_idx++) {
      // Loop over this variable's elements
      for (int elem_idx = 0; elem_idx < out_data[dataFiles[file_idx]->varid[var_idx]].nelem; elem_idx++) {
        if (!(var_idx == 0 && elem_idx == 0)) {
          fprintf(dataFiles[file_idx]->fh, "\t ");
        }
        fprintf(dataFiles[file_idx]->fh, out_data[dataFiles[file_idx]->varid[var_idx]].format, out_data[dataFiles[file_idx]->varid[var_idx]].aggdata[elem_idx]);
      }
    }
    fprintf(dataFiles[file_idx]->fh, "\n");
  }
}

void WriteOutputNetCDF::compressFiles() {
  //TODO: implement this
}

void WriteOutputNetCDF::write_header(out_data_struct* out_data, const dmy_struct* dmy,
    const ProgramState* state) {
}


