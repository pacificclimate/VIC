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
  void notifyCellLocation(float lat, float lng);
  void notifyDimensionUpdate(StateVariables::StateVariableLastDimension dimension, int value = -1);
  int seekToCell(int cellid, int* nVeg, int* nBand);
  void flush();
  void rewindFile();

private:
  template<typename T> int generalWrite(const T* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  template<typename T> int generalRead(T* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  void populateMetaData();
  void populateMetaDimensions();
  void initializeDimensionIndices();
  std::map<StateVariables::StateMetaDataVariableIndices, StateVariableMetaData> metaData;
  std::map<StateVariables::StateVariableLastDimension, StateVariableDimension> metaDimensions;
  void checkAndRead(std::string name, int * storage, int length);

  void openFile();
  void closeFile();

  netCDF::NcFile* netCDF;
  std::string filename;
  int curLatIndex;
  int curLonIndex;
  std::map<StateVariables::StateVariableLastDimension, int> curDimensionIndices;
};

#endif /* STATEIONETCDF_H_ */
