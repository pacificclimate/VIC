#include "StateIOBinary.h"

#include <cstdio>

StateIOBinary::StateIOBinary(FILE* file, const ProgramState* state) : StateIO(state), file(file) {
}

StateIOBinary::~StateIOBinary() {
}

int StateIOBinary::write(const int* data, int numValues, const StateVariableMetaData* meta) {
  return fwrite(data, sizeof(int), numValues, file);
}

int StateIOBinary::write(const double* data, int numValues, const StateVariableMetaData* meta) {
  return fwrite(data, sizeof(double), numValues, file);
}

int StateIOBinary::write(const char* data, int numValues, const StateVariableMetaData* meta) {
  return fwrite(data, sizeof(char), numValues, file);
}

int StateIOBinary::read(int* data, int numValues) {
  return 0;
}



int StateIOBinary::read(double* data, int numValues) {
  return 0;
}


