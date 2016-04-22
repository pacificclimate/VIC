#ifndef STATEIONETCDF_H_
#define STATEIONETCDF_H_
#include "user_def.h"   // Need to know the NETCDF_OUTPUT_AVAILABLE option here.

#if NETCDF_OUTPUT_AVAILABLE

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
  int write(const float* data, int numValues, const StateVariables::StateMetaDataVariableIndices id); //new
  int write(const bool* data, int numValues, const StateVariables::StateMetaDataVariableIndices id); //new
  int write(const char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int read(int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int read(double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int read(float* data, int numValues, const StateVariables::StateMetaDataVariableIndices id); //new
  int read(bool* data, int numValues, const StateVariables::StateMetaDataVariableIndices id); //new
  int read(char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  StateHeader readHeader();
  void notifyDimensionUpdate(StateVariables::StateVariableDimensionId dimension, int value = -1);
  void initializeDimensionIndices();
  int getCurrentDimensionIndex(StateVariables::StateVariableDimensionId dimension);
  int seekToCell(int cellid, int* nVeg, int* nBand);
  void flush();
  void rewindFile();

private:
  template<typename T> int generalWrite(const T* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  template<typename T> int generalRead(T* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  void populateMetaData(const ProgramState* state);
  void populateMetaDimensions();
  void initializeLatLonDims();
  std::map<StateVariables::StateMetaDataVariableIndices, StateVariableMetaData> metaData;
  std::map<StateVariables::StateVariableDimensionId, StateVariableDimension> metaDimensions;
  void checkAndRead(std::string name, int * storage, int length);

  void openFile();
  void closeFile();

  netCDF::NcFile* netCDF;
  std::map<StateVariables::StateVariableDimensionId, int> curDimensionIndices;
};

#endif // NETCDF_OUTPUT_AVAILABLE

#endif /* STATEIONETCDF_H_ */
