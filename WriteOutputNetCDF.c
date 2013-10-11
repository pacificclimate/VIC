#include "WriteOutputNetCDF.h"
#include "netcdfcpp.h"

WriteOutputNetCDF::WriteOutputNetCDF(const ProgramState* state) : WriteOutputFormat(state), netCDF(NULL) {
  netCDFOutputFileName = state->options.NETCDF_FULL_FILE_PATH;
}

WriteOutputNetCDF::~WriteOutputNetCDF() {
  fprintf(stderr, "NetCDF output file closed from destructor\n");
  if (netCDF != NULL) {
    delete netCDF;
  }
}

const char* WriteOutputNetCDF::getDescriptionOfOutputType() {
  return "NETCDF";
}

// Called once for each grid cell. This could potentially support parallel writes if a library supports it.
void WriteOutputNetCDF::openFile() {
  netCDF = new NcFile(netCDFOutputFileName.c_str(), NcFile::Write, NULL, 0, NcFile::Netcdf4);
}

// Called automatically only once at the beginning of the program.
void WriteOutputNetCDF::initializeFile() {
  NcFile ncFile(netCDFOutputFileName.c_str(), NcFile::Replace, NULL, 0, NcFile::Netcdf4);

  if (!ncFile.is_valid()) {
    char errorString [MAXSTRING];
    sprintf(errorString, "Fatal Error: could not initialise netCDF file for output. File name: \"%s\"\n", netCDFOutputFileName.c_str());
    vicerror(errorString);
  }

  static const int NDIMS = 2;
  static const int NX = 6;
  static const int NY = 12;
  // This is the data array we will write. It will just be filled
  // with a progression of numbers for this example.
  int dataOut[NX][NY];

  // Create some pretend data. If this wasn't an example program, we
  // would have some real data to write, for example, model output.
  for(int i = 0; i < NX; i++)
     for(int j = 0; j < NY; j++)
  dataOut[i][j] = i * NY + j;

  // Add global attributes.
  ncFile.add_att("institution", "Pacific Climate Impacts Consortium (PCIC), Victoria, BC, www.pacificclimate.org");

  // When we create netCDF dimensions, we get back a pointer to an
  // NcDim for each one.
  NcDim* xDim = ncFile.add_dim("x", NX);
  NcDim* yDim = ncFile.add_dim("y", NY);

  // Define a netCDF variable. The type of the variable in this case
  // is ncInt (32-bit integer).
  NcVar *data = ncFile.add_var("data", ncInt, xDim, yDim);

  // Write the pretend data to the file. Although netCDF supports
  // reading and writing subsets of data, in this case we write all
  // the data in one operation.
  data->put(&dataOut[0][0], NX, NY);


/*
  // Set up the dimensions and variables.
  NcDim* timeDim = ncFile.add_dim("time");    // Time is unlimited.
  NcDim* latDim = ncFile.add_dim("lat", 500); // TODO: where do these dimension sizes come from?
  NcDim* lonDim = ncFile.add_dim("lon", 500);*/

  ncFile.sync();
  ncFile.close();
  // The file will be automatically close when the NcFile object goes
  // out of scope. This frees up any internal netCDF resources
  // associated with the file, and flushes any buffers.
}

// This is called per time record.
void WriteOutputNetCDF::write_data(out_data_struct* out_data, const dmy_struct* dmy, int dt,
    const ProgramState* state) {

  if (!netCDF->is_valid()) {
    fprintf(stderr, "Warning: could not write to netCDF file. Skipping data. Record %04i\t%02i\t%02i\t%02i\t Lat: %f, Lon: %f. File: \"%s\".",
        dmy->year, dmy->month, dmy->day, dmy->hour, lat, lon, netCDFOutputFileName.c_str());
    return;
  }
/*
  // Loop over output files
  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {

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
  */
}

void WriteOutputNetCDF::compressFiles() {
  //TODO: implement this
}

void WriteOutputNetCDF::write_header(out_data_struct* out_data, const dmy_struct* dmy,
    const ProgramState* state) {
  // This is not applicable for netCDF output format. Intentionally empty.
}


