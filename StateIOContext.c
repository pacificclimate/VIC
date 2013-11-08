#include "StateIOContext.h"

#include "StateIOASCII.h"
#include "StateIOBinary.h"
#include "StateIONetCDF.h"

StateIOContext::StateIOContext(std::string filename, const ProgramState* state) : stream(NULL) {
  if (state->options.STATE_FORMAT == StateOutputFormat::BINARY_STATEFILE) {
    stream = new StateIOBinary(filename, state);
  } else if (state->options.STATE_FORMAT == StateOutputFormat::NETCDF_STATEFILE) {
    stream = new StateIONetCDF(filename, state);
  } else {
    stream = new StateIOASCII(filename, state);
  }
}

StateIOContext::~StateIOContext() {
  if (stream != NULL) {
    delete stream;
  }
}

