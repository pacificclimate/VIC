#ifndef STATEIOASCII_H_
#define STATEIOASCII_H_

#include <string>

#include "StateIO.h"

class StateIOASCII: public StateIO {
public:
  StateIOASCII(std::string filename, IOType ioType, const ProgramState* state);
  virtual ~StateIOASCII();
  void initializeOutput();
  int write(const int* data, int numValues, const StateVariableMetaData* meta);
  int write(const double* data, int numValues, const StateVariableMetaData* meta);
  int write(const char* data, int numValues, const StateVariableMetaData* meta);
  int writeNewline();
  StateHeader readHeader();
  int seekToCell(int cellid, int* nVeg, int* nBand);
  int read(int* data, int numValues, const StateVariableMetaData* meta);
  int read(double* data, int numValues, const StateVariableMetaData* meta);
  int read(char* data, int numValues, const StateVariableMetaData* meta);
  void flush();
  void rewindFile();
private:
  FILE* file;
  bool firstValueOnLine;
};

#endif /* STATEIOASCII_H_ */
