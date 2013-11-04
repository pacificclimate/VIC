#include "WriteOutputContext.h"
#include "WriteOutputAscii.h"
#include "WriteOutputBinary.h"
#include "WriteOutputNetCDF.h"

WriteOutputContext::WriteOutputContext(const ProgramState* state) {
  if (state->options.OUTPUT_FORMAT == OutputFormat::BINARY_FORMAT) {
    outputFormat = new WriteOutputBinary(state);
#if NETCDF_OUTPUT_AVAILABLE
  } else if (state->options.OUTPUT_FORMAT == OutputFormat::NETCDF_FORMAT) {
    outputFormat = new WriteOutputNetCDF(state);
#endif
  } else {
    outputFormat = new WriteOutputAscii(state);
  }
}

WriteOutputContext::~WriteOutputContext() {
  delete outputFormat;
}

