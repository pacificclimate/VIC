#include "StateIOBinary.h"

#include <cstdio>
#include <stdio.h>

#include "vicNl.h"

StateIOBinary::StateIOBinary(std::string filename, IOType ioType, const ProgramState* state) : StateIO(filename, ioType, state) {
  std::string openType = "rb";
  if (ioType == StateIO::Writer) {
    openType = "ab";
  }
  file = open_file(filename.c_str(), openType.c_str());
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

int StateIOBinary::write(const int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  //return fwrite(data, sizeof(int), numValues, file);
  int dataLength = numValues * sizeof(int);
  dataToWrite.append((const char*)data, dataLength);
  return dataLength;
}

int StateIOBinary::write(const double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  //return fwrite(data, sizeof(double), numValues, file);
  int dataLength = numValues * sizeof(double);
  dataToWrite.append((const char*)data, dataLength);
  return dataLength;
}

int StateIOBinary::write(const char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  //return fwrite(data, sizeof(char), numValues, file);
  int dataLength = numValues * sizeof(char);
  dataToWrite.append((const char*)data, dataLength);
  return dataLength;
}

int StateIOBinary::read(int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  int numRead = fread(data, sizeof(int), numValues, file);
  if (feof(file)) {
    throw new VICException("End of model state file found unexpectedly");
  }
  return numRead;
}

int StateIOBinary::read(double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  int numRead = fread(data, sizeof(double), numValues, file);
  if (feof(file)) {
    throw new VICException("End of model state file found unexpectedly");
  }
  return numRead;
}

int StateIOBinary::read(char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  int numRead = fread(data, sizeof(char), numValues, file);
  if (feof(file)) {
    throw new VICException("End of model state file found unexpectedly");
  }
  return numRead;
}

StateHeader StateIOBinary::readHeader() {
  int year, month, day, nLayer, nNode;
  fread(&year, sizeof(int), 1, file);
  fread(&month, sizeof(int), 1, file);
  fread(&day, sizeof(int), 1, file);
  fread(&nLayer, sizeof(int), 1, file);
  fread(&nNode, sizeof(int), 1, file);
  return StateHeader(year, month, day, nLayer, nNode);
}

int StateIOBinary::seekToCell(int cellid, int* nVeg, int* nBand) {
  int tmpCellNum, tmpNVeg, tmpNBand, tmpNBytes;
  char tmpByte;
  /* read cell information */
  fread(&tmpCellNum, sizeof(int), 1, file);
  fread(&tmpNVeg, sizeof(int), 1, file);
  fread(&tmpNBand, sizeof(int), 1, file);
  fread(&tmpNBytes, sizeof(int), 1, file);
  // Skip over unused cell information
  while (tmpCellNum != cellid && !feof(file)) {
    // skip rest of current cells info
    for (int byte = 0; byte < tmpNBytes; byte++)
      fread(&tmpByte, 1, 1, file);
    // read info for next cell
    fread(&tmpCellNum, sizeof(int), 1, file);
    fread(&tmpNVeg, sizeof(int), 1, file);
    fread(&tmpNBand, sizeof(int), 1, file);
    fread(&tmpNBytes, sizeof(int), 1, file);
  } //end while

  *nVeg = tmpNVeg;
  *nBand = tmpNBand;

  if (feof(file)) {
    return -1;
  }

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
  if (ioType == StateIO::Writer) {
    int headerLength = 3 * sizeof(int); // Three integers at the beginning of the line are not counted.
    int numBytes = (dataToWrite.size() - headerLength);
    dataToWrite.insert(headerLength, (const char *) &numBytes, sizeof(int));
    fwrite(dataToWrite.c_str(), sizeof(char), dataToWrite.size(), file);
    fflush(file);
    dataToWrite = "";
  }
}

void StateIOBinary::rewindFile() {
  rewind(file);
}


