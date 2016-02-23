#include "WriteOutputContext.h"
#include "WriteOutputAscii.h"
#include "WriteOutputBinary.h"
#include "WriteOutputNetCDF.h"

WriteOutputContext::WriteOutputContext(const ProgramState* state) {
    outputFormat = new WriteOutputNetCDF(state);
}

WriteOutputContext::~WriteOutputContext() {
  delete outputFormat;
}

