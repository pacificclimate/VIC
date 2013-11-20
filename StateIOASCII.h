#ifndef STATEIOASCII_H_
#define STATEIOASCII_H_

#include <string>

#include "StateIO.h"

class StateIOASCII: public StateIO {
public:
  StateIOASCII(std::string filename, IOType ioType, const ProgramState* state);
  virtual ~StateIOASCII();
  void initializeOutput();
  int write(const int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int write(const double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int write(const char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int processNewline();
  StateHeader readHeader();
  int seekToCell(int cellid, int* nVeg, int* nBand);
  int read(int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int read(double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int read(char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  void flush();
  void rewindFile();
private:
  FILE* file;
  bool firstValueOnLine;
};

#endif /* STATEIOASCII_H_ */
