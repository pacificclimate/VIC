#include "StateIOASCII.h"

#include <cstdio>

StateIOASCII::StateIOASCII(FILE* file, const ProgramState* state) :  StateIO(state), file(file), firstValueOnLine(true) {
  // TODO Auto-generated constructor stub

}

StateIOASCII::~StateIOASCII() {
  // TODO Auto-generated destructor stub
}

int StateIOASCII::write(const int* data, int numValues, const StateVariableMetaData* meta) {
  int errorCode = 0;
  for (int i = 0; i < numValues; i++) {
    if (!firstValueOnLine) {
      errorCode += fprintf(file, " ");
    }
    firstValueOnLine = false;
    errorCode += fprintf(file, "%i", data[i]);
  }
  return errorCode;
}

int StateIOASCII::write(const double* data, int numValues, const StateVariableMetaData* meta) {
  int errorCode = 0;
  for (int i = 0; i < numValues; i++) {
    if (!firstValueOnLine) {
      errorCode += fprintf(file, " ");
    }
    firstValueOnLine = false;
    errorCode += fprintf(file, "%f", data[i]);
  }
  return errorCode;
}

int StateIOASCII::write(const char* data, int numValues, const StateVariableMetaData* meta) {
  int errorCode = 0;
  for (int i = 0; i < numValues; i++) {
    if (!firstValueOnLine) {
      errorCode += fprintf(file, " ");
    }
    firstValueOnLine = false;
    errorCode += fprintf(file, "%i", (int) data[i]);
  }
  return errorCode;
}

int StateIOASCII::writeNewline() {
  int errorCode = fprintf(file, "\n");
  firstValueOnLine = true;
  return errorCode;
}

int StateIOASCII::read(int* data, int numValues) {
  return 0;
}


int StateIOASCII::read(double* data, int numValues) {
  return 0;
}

void StateIOASCII::flush() {
  fflush(file);
}



