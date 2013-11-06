#ifndef STATEIOBINARY_H_
#define STATEIOBINARY_H_

#include "StateIO.h"

class StateIOBinary: public StateIO {
public:
  StateIOBinary(FILE* file, const ProgramState* state);
  virtual ~StateIOBinary();
  int write(const int* data, int numValues, const StateVariableMetaData* meta);
  int write(const double* data, int numValues, const StateVariableMetaData* meta);
  int write(const char* data, int numValues, const StateVariableMetaData* meta);
  int read(int* data, int numValues);
  int read(double* data, int numValues);
  void flush();
private:
  FILE* file;
  std::string dataToWrite;
};

#endif /* STATEIOBINARY_H_ */
