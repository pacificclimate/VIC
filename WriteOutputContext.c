#include "WriteOutputContext.h"
#include "WriteOutputAscii.h"
#include "WriteOutputBinary.h"
#include "WriteOutputNetCDF.h"

WriteOutputContext::WriteOutputContext(OutputFormat::Type type) {
  if (type == OutputFormat::BINARY_FORMAT) {
    outputFormat = new WriteOutputBinary();
  } else if (type == OutputFormat::NETCDF_FORMAT) {
    outputFormat = new WriteOutputNetCDF();
  } else {
    outputFormat = new WriteOutputAscii();
  }
}

WriteOutputContext::~WriteOutputContext() {
  delete outputFormat;
}

