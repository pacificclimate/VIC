#include "StateIOBinary.h"

#include <cstdio>

StateIOBinary::StateIOBinary(FILE* file, const ProgramState* state) : StateIO(state), file(file) {
}

StateIOBinary::~StateIOBinary() {
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
  fputs(dataToWrite.c_str(),file);
  fflush(file);
  dataToWrite = "";
}



