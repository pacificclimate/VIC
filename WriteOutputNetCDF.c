#include "WriteOutputNetCDF.h"

FILE* WriteOutputNetCDF::openFile(const char* string) {
  return NULL;//TODO: implement this
}

void WriteOutputNetCDF::write_data(out_data_file_struct* out_data_files,
    out_data_struct* out_data, const dmy_struct* dmy, int dt,
    const ProgramState* state) {
}

const char* WriteOutputNetCDF::getDescriptionOfOutputType() {
  return "NETCDF";
}

void WriteOutputNetCDF::write_header(out_data_file_struct* out_data_files,
    out_data_struct* out_data, const dmy_struct* dmy,
    const ProgramState* state) {
}


