#include "StateIONetCDF.h"

StateIONetCDF::StateIONetCDF(std::string filename, IOType ioType, const ProgramState* state) : StateIO(filename, ioType, state) {
  // TODO Auto-generated constructor stub

}

StateIONetCDF::~StateIONetCDF() {
  // TODO Auto-generated destructor stub
}

void StateIONetCDF::initializeOutput() {
  //TODO: initialize global attributes, setup variable dimensions and structure.
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

int StateIONetCDF::read(int* data, int numValues, const StateVariableMetaData* meta) {
  return 0;
}

int StateIONetCDF::read(double* data, int numValues, const StateVariableMetaData* meta) {
  return 0;
}

int StateIONetCDF::read(char* data, int numValues, const StateVariableMetaData* meta) {
  return 0;
}

StateHeader StateIONetCDF::readHeader() {
  return StateHeader(-1, -1, -1, -1, -1);
}

int StateIONetCDF::seekToCell(int cellid, int* nVeg, int* nBand) {
  return 0;
}

void StateIONetCDF::flush() {
  //TODO: implement if applicable
}

void StateIONetCDF::rewindFile() {

}


