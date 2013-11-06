#include "StateIONetCDF.h"

StateIONetCDF::StateIONetCDF(const ProgramState* state) : StateIO(state) {
  // TODO Auto-generated constructor stub

}

StateIONetCDF::~StateIONetCDF() {
  // TODO Auto-generated destructor stub
}

int StateIONetCDF::write(const int* data, int numValues,
    const StateVariableMetaData* meta) {
  return 0;
}

int StateIONetCDF::write(const double* data, int numValues,
    const StateVariableMetaData* meta) {
  return 0;
}

int StateIONetCDF::write(const char* data, int numValues,
    const StateVariableMetaData* meta) {
  return 0;
}

int StateIONetCDF::read(int* data, int numValues) {
  return 0;
}



int StateIONetCDF::read(double* data, int numValues) {
  return 0;
}


