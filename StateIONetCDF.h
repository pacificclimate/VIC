#ifndef STATEIONETCDF_H_
#define STATEIONETCDF_H_

#include <map>
#include <string>

#include "StateIO.h"

namespace netCDF {
  class NcFile;
}

class StateIONetCDF: public StateIO {
public:
  StateIONetCDF(std::string filename, IOType ioType, const ProgramState* state);
  virtual ~StateIONetCDF();
  void initializeOutput();
  int write(const int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int write(const double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int write(const char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int read(int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int read(double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int read(char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  StateHeader readHeader();
  int seekToCell(int cellid, int* nVeg, int* nBand);
  void flush();
  void rewindFile();

private:
  void populateMetaData();
  std::map<int, StateVariableMetaData> metaData;
  void checkAndRead(std::string name, int * storage, int length);

  void openFile();
  void closeFile();

  netCDF::NcFile* netCDF;
  std::string filename;
};

#endif /* STATEIONETCDF_H_ */
