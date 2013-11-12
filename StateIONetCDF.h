#ifndef STATEIONETCDF_H_
#define STATEIONETCDF_H_

#include <string>

#include "StateIO.h"

class StateIONetCDF: public StateIO {
public:
  StateIONetCDF(std::string filename, IOType ioType, const ProgramState* state);
  virtual ~StateIONetCDF();
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
};

#endif /* STATEIONETCDF_H_ */
