#ifndef STATEIOASCII_H_
#define STATEIOASCII_H_

#include <string>

#include "StateIO.h"

class StateIOASCII: public StateIO {
public:
  StateIOASCII(std::string filename, const ProgramState* state);
  virtual ~StateIOASCII();
  void initializeOutput();
  int write(const int* data, int numValues, const StateVariableMetaData* meta);
  int write(const double* data, int numValues, const StateVariableMetaData* meta);
  int write(const char* data, int numValues, const StateVariableMetaData* meta);
  int writeNewline();
  int read(int* data, int numValues);
  int read(double* data, int numValues);
  void flush();
private:
  FILE* file;
  bool firstValueOnLine;
};

#endif /* STATEIOASCII_H_ */
