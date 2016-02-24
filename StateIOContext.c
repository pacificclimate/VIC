#include "StateIOContext.h"
#include "user_def.h"
#include "StateIOASCII.h"
#include "StateIOBinary.h"
#include "StateIONetCDF.h"


StateIOContext::StateIOContext(std::string filename, StateIO::IOType ioType, const ProgramState* state) : stream(NULL) {
	if (state->options.STATE_FORMAT == StateOutputFormat::NETCDF_STATEFILE) {
		stream = new StateIONetCDF(filename, ioType, state);
	}
}

StateIOContext::~StateIOContext() {
	if (stream != NULL) {
		delete stream;
	}
}

