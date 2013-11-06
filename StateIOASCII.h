#ifndef STATEIOASCII_H_
#define STATEIOASCII_H_

#include "StateIO.h"

class StateIOASCII: public StateIO {
public:
  StateIOASCII(FILE* file, const ProgramState* state);
  virtual ~StateIOASCII();
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
