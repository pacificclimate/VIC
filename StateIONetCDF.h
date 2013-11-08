#ifndef STATEIONETCDF_H_
#define STATEIONETCDF_H_

#include "StateIO.h"

class StateIONetCDF: public StateIO {
public:
  StateIONetCDF(const ProgramState* state);
  virtual ~StateIONetCDF();
  void initializeOutput(FILE** f, const char* filename, const ProgramState* state);
  int write(const int* data, int numValues, const StateVariableMetaData* meta);
  int write(const double* data, int numValues, const StateVariableMetaData* meta);
  int write(const char* data, int numValues, const StateVariableMetaData* meta);
  int read(int* data, int numValues);
  int read(double* data, int numValues);
  void flush();
};

#endif /* STATEIONETCDF_H_ */
