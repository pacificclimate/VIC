#ifndef STATEIOBINARY_H_
#define STATEIOBINARY_H_

#include <string>

#include "StateIO.h"

class StateIOBinary: public StateIO {
public:
  StateIOBinary(std::string filename, IOType ioType, const ProgramState* state);
  virtual ~StateIOBinary();
  void initializeOutput();
  int write(const int* data, int numValues, const StateVariableMetaData* meta);
  int write(const double* data, int numValues, const StateVariableMetaData* meta);
  int write(const char* data, int numValues, const StateVariableMetaData* meta);
  int read(int* data, int numValues, const StateVariableMetaData* meta);
  int read(double* data, int numValues, const StateVariableMetaData* meta);
  int read(char* data, int numValues, const StateVariableMetaData* meta);
  StateHeader readHeader();
  int seekToCell(int cellid, int* nVeg, int* nBand);
  void flush();
  void rewindFile();
private:
  FILE* file;
  std::string dataToWrite;
};

#endif /* STATEIOBINARY_H_ */
