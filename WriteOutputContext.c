#include "WriteOutputContext.h"
#include "WriteOutputAscii.h"
#include "WriteOutputBinary.h"
#include "WriteOutputNetCDF.h"

WriteOutputContext::WriteOutputContext(const ProgramState* state) {
  if (state->options.OUTPUT_FORMAT == OutputFormat::BINARY_FORMAT) {
    outputFormat = new WriteOutputBinary(state);
  } else if (state->options.OUTPUT_FORMAT == OutputFormat::NETCDF_FORMAT) {
    outputFormat = new WriteOutputNetCDF(state);
  } else {
    outputFormat = new WriteOutputAscii(state);
  }
}

WriteOutputContext::~WriteOutputContext() {
  delete outputFormat;
}

