#include "StateIOBinary.h"

#include <cstdio>

#include "vicNl.h"

StateIOBinary::StateIOBinary(std::string filename, const ProgramState* state) : StateIO(filename, state) {
  file = open_file(filename.c_str(), "ab");
}

StateIOBinary::~StateIOBinary() {
  if (file != NULL) {
    fclose(file);
  }
}

void StateIOBinary::initializeOutput() {
  // Wipe out any existing file of the same name if it exists.
  fclose(file);
  file = open_file(filename.c_str(), "wb");

  // The write functions are not used here because the automatic NBytes field is not applicable for the file header.

  /* Write save state date information */
  fwrite(&state->global_param.stateyear, sizeof(int), 1, file);
  fwrite(&state->global_param.statemonth, sizeof(int), 1, file);
  fwrite(&state->global_param.stateday, sizeof(int), 1, file);

  /* Write simulation flags */
  fwrite(&state->options.Nlayer, sizeof(int), 1, file);
  fwrite(&state->options.Nnode, sizeof(int), 1, file);

  fflush(file);
}

int StateIOBinary::write(const int* data, int numValues, const StateVariableMetaData* meta) {
  //return fwrite(data, sizeof(int), numValues, file);
  int dataLength = numValues * sizeof(int);
  dataToWrite.append((const char*)data, dataLength);
  return dataLength;
}

int StateIOBinary::write(const double* data, int numValues, const StateVariableMetaData* meta) {
  //return fwrite(data, sizeof(double), numValues, file);
  int dataLength = numValues * sizeof(double);
  dataToWrite.append((const char*)data, dataLength);
  return dataLength;
}

int StateIOBinary::write(const char* data, int numValues, const StateVariableMetaData* meta) {
  //return fwrite(data, sizeof(char), numValues, file);
  int dataLength = numValues * sizeof(char);
  dataToWrite.append((const char*)data, dataLength);
  return dataLength;
}

int StateIOBinary::read(int* data, int numValues) {
  return 0;
}

int StateIOBinary::read(double* data, int numValues) {
  return 0;
}

/*
 * The StateIOBinary takes care of finding out how much data is being written per cell internally (here).
 * The expected format for binary statefile per cell is: cellNum, numVeg, numBands, numBytes, data...
 * where numBytes is the number of bytes from after this value to the end of the line.
 * Instead of counting the bytes of each value that will be inserted to the file manually and writing this
 * before everything else, the Binary Writer just appends all values to an internal string and counts everything
 * when the data is flushed. This is much better for maintenance since the list of values added doesn't need to be updated
 * manually by hand every time a new value is added to the state file.
 */
void StateIOBinary::flush() {
  int headerLength = 3 * sizeof(int); // Three integers at the beginning of the line are not counted.
  int numBytes = (dataToWrite.size() - headerLength);
  dataToWrite.insert(headerLength, (const char *)&numBytes, sizeof(int));
  fwrite(dataToWrite.c_str(), sizeof(char), dataToWrite.size() ,file);
  fflush(file);
  dataToWrite = "";
}



