#include "StateIOASCII.h"

#include <stdio.h>
#include <cstdio>
#include <sstream>

#include "vicNl.h"

StateIOASCII::StateIOASCII(std::string filename, IOType ioType, const ProgramState* state) :  StateIO(filename, ioType, state), firstValueOnLine(true) {
  std::string openType = "r";
  if (ioType == StateIO::Writer) {
    openType = "a";
  }
  file = open_file(filename.c_str(), openType.c_str());
}

StateIOASCII::~StateIOASCII() {
  if (file != NULL) {
    fclose(file);
  }
}

void StateIOASCII::initializeOutput() {
  // Wipe out any existing file of the same name if it exists.
  fclose(file);
  file = open_file(filename.c_str(), "w");

  /* Write save state date information */
  write(&state->global_param.stateyear, 1, StateVariables::NONE);
  write(&state->global_param.statemonth, 1, StateVariables::NONE);
  write(&state->global_param.stateday, 1, StateVariables::NONE);
  processNewline();
  /* Write simulation flags */
  write(&state->options.Nlayer, 1, StateVariables::NONE);
  write(&state->options.Nnode, 1, StateVariables::NONE);
  processNewline();
}

int StateIOASCII::write(const int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
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

int StateIOASCII::write(const double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  int errorCode = 0;
  for (int i = 0; i < numValues; i++) {
    if (!firstValueOnLine) {
      errorCode += fprintf(file, " ");
    }
    firstValueOnLine = false;
    errorCode += fprintf(file, "%.18e", data[i]);
  }
  return errorCode;
}

int StateIOASCII::write(const float* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) { //new
  int errorCode = 0;
  for (int i = 0; i < numValues; i++) {
    if (!firstValueOnLine) {
      errorCode += fprintf(file, " ");
    }
    firstValueOnLine = false;
    errorCode += fprintf(file, "%.18e", data[i]);
  }
  return errorCode;
}

int StateIOASCII::write(const bool* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) { //new
  int errorCode = 0;
  for (int i = 0; i < numValues; i++) {
    if (!firstValueOnLine) {
      errorCode += fprintf(file, " ");
    }
    firstValueOnLine = false;
    errorCode += fprintf(file, "%d", data[i]);
  }
  return errorCode;
}

int StateIOASCII::write(const char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
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

int StateIOASCII::processNewline() {
  if (ioType == StateIO::Writer) {
    int errorCode = fprintf(file, "\n");
    firstValueOnLine = true;
    return errorCode;
  }
  return 0;
}


int StateIOASCII::read(int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  int numRead = 0;
  for (int i = 0; i < numValues; i++) {
    numRead += fscanf(file, "%d", &data[i]);
  }
  if (feof(file)) {
    throw VICException("End of model state file found unexpectedly");
  }
  return numRead;
}

int StateIOASCII::read(double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  int numRead = 0;
  for (int i = 0; i < numValues; i++) {
    numRead += fscanf(file, "%lf", &data[i]);  //TODO: possibly read " %lf"
  }
  if (feof(file)) {
    throw VICException("End of model state file found unexpectedly");
  }
  return numRead;
}

int StateIOASCII::read(float* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) { //new
  int numRead = 0;
  for (int i = 0; i < numValues; i++) {
    numRead += fscanf(file, "%f", &data[i]);
  }
  if (feof(file)) {
    throw VICException("End of model state file found unexpectedly");
  }
  return numRead;
}

int StateIOASCII::read(bool* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) { //new
  int numRead = 0;
  for (int i = 0; i < numValues; i++) {
    numRead += fscanf(file, "%d", &data[i]);
  }
  if (feof(file)) {
    throw VICException("End of model state file found unexpectedly");
  }
  return numRead;
}

int StateIOASCII::read(char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  int numRead = 0;
  for (int i = 0; i < numValues; i++) {
    int temp;
    numRead += fscanf(file, "%d", &temp);
    data[i] = (char)temp;
  }
  if (feof(file)) {
    throw VICException("End of model state file found unexpectedly");
  }
  return numRead;
}

StateHeader StateIOASCII::readHeader() {
  int year, month, day, nLayer, nNode;
  fscanf(file, "%d %d %d", &year, &month, &day);
  fscanf(file, "%d %d", &nLayer, &nNode);
  return StateHeader(year, month, day, nLayer, nNode);
}

int StateIOASCII::seekToCell(int cellid, int* nVeg, int* nBand) {
  int tmpCellNum, tmpNveg, tmpNband;
  char buffer[MAXSTRING];
  /* read cell information */
  fscanf(file, "%d %d %d", &tmpCellNum, &tmpNveg, &tmpNband);
  // Skip over unused cell information
  while (tmpCellNum != cellid && !feof(file)) {
    // skip rest of current cells info
    fgets(buffer, MAXSTRING, file); // skip rest of general cell info
#if EXCESS_ICE
        fgets(buffer, MAXSTRING, file); //excess ice info
#endif
    for (int veg = 0; veg <= tmpNveg; veg++) {
      fgets(buffer, MAXSTRING, file); // skip dist precip info
      for (int band = 0; band < tmpNband; band++)
        fgets(buffer, MAXSTRING, file); // skip snowband info
    }
    if (state->options.LAKES) {
      fgets(buffer, MAXSTRING, file); // skip lake info
    }
    // read info for next cell
    fscanf(file, "%d %d %d", &tmpCellNum, &tmpNveg, &tmpNband);
  } //end while

  *nVeg = tmpNveg;
  *nBand = tmpNband;

  if (feof(file)) {
    return -1;
  }

  return 0;
}

void StateIOASCII::flush() {
  if (ioType == StateIO::Writer) {
    fflush(file);
  }
}

void StateIOASCII::rewindFile() {
  rewind(file);
}




