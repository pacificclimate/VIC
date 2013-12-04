#ifndef STATEIO_H_
#define STATEIO_H_

#include "vicNl_def.h"
#include <string>
#include <vector>
#if NETCDF_OUTPUT_AVAILABLE
#include <netcdf>
#endif // NETCDF_OUTPUT_AVAILABLE

// For backwards compatibility, only add values to the end of the enum,
// and never remove values (just stop using them if necessary).
namespace StateVariables {
enum StateMetaDataVariableIndices {
  NONE = 0,
  GRID_CELL,
  VEG_TYPE_NUM,
  NUM_BANDS,
  SOIL_DZ_NODE,
  SOIL_ZSUM_NODE,
  SOIL_DEPTH,
  SOIL_EFFECTIVE_POROSITY,
  SOIL_DP,
  PRCP_MU,
  INIT_STILL_STORM,
  INIT_DRY_TIME,
  HRU_VEG_INDEX,
  HRU_BAND_INDEX,
  LAYER_MOIST,
  LAYER_SOIL_ICE,
  LAYER_ICE_CONTENT,
  HRU_VEG_VAR_WDEW,
  SNOW_LAST_SNOW,
  SNOW_MELTING,
  SNOW_COVERAGE,
  SNOW_SWQ,
  SNOW_SURF_TEMP,
  SNOW_SURF_WATER,
  SNOW_PACK_TEMP,
  SNOW_PACK_WATER,
  SNOW_DENSITY,
  SNOW_COLD_CONTENT,
  SNOW_CANOPY,
  ENERGY_T,
  LAKE_LAYER_MOIST,
  LAKE_LAYER_SOIL_ICE,
  LAKE_LAYER_ICE_CONTENT,
  LAKE_SNOW_LAST_SNOW,
  LAKE_SNOW_MELTING,
  LAKE_SNOW_COVERAGE,
  LAKE_SNOW_SWQ,
  LAKE_SNOW_SURF_TEMP,
  LAKE_SNOW_SURF_WATER,
  LAKE_SNOW_PACK_TEMP,
  LAKE_SNOW_PACK_WATER,
  LAKE_SNOW_DENSITY,
  LAKE_SNOW_COLD_CONTENT,
  LAKE_SNOW_CANOPY,
  LAKE_ENERGY_T,
  LAKE_ACTIVENOD,
  LAKE_DZ,
  LAKE_SURFDZ,
  LAKE_LDEPTH,
  LAKE_SURFACE,
  LAKE_SAREA,
  LAKE_VOLUME,
  LAKE_TEMP,
  LAKE_TEMPAVG,
  LAKE_AREAI,
  LAKE_NEW_ICE_AREA,
  LAKE_ICE_WATER_EQ,
  LAKE_HICE,
  LAKE_TEMPI,
  LAKE_SWE,
  LAKE_SURF_TEMP,
  LAKE_PACK_TEMP,
  LAKE_SALBEDO,
  LAKE_SDEPTH
};

enum StateVariableDimensionId {
  NO_DIM = 0,
  LAT_DIM,
  LON_DIM,
  BNDS_DIM,
  LAYERS_DIM,
  NODES_DIM,
  LAKE_NODES_DIM,
  FROST_LAYER_AREAS_DIM,
  FROST_AREAS_DIM,
  HRU_DIM,
  VEG_DIM,
  DIST_DIM,
};
}

#if NETCDF_OUTPUT_AVAILABLE

using std::string;
class StateVariableDimension {
public:
  StateVariableDimension() : name("invalid"), size(-1) {}
  StateVariableDimension(string name, int size) : name(name), size(size) {}
  string name;
  int size;
};

using StateVariables::StateVariableDimensionId;
using StateVariables::NO_DIM;
// This is just a wrapper class for now in case more attributes are needed in netCDF state files.
class StateVariableMetaData {
public:
  StateVariableMetaData() : name("invalid"), type(netCDF::NcType::nc_DOUBLE) {}
  StateVariableMetaData(string name, StateVariableDimensionId d1 = NO_DIM, StateVariableDimensionId d2 = NO_DIM,
      StateVariableDimensionId d3 = NO_DIM, StateVariableDimensionId d4 = NO_DIM) : name(name), type(netCDF::NcType::nc_DOUBLE) {
    dimensions.push_back(StateVariables::LAT_DIM);  // All variables are index by lat/long.
    dimensions.push_back(StateVariables::LON_DIM);
    dimensions.push_back(d1);
    dimensions.push_back(d2);
    dimensions.push_back(d3);
    dimensions.push_back(d4);
  }
  string name;
  netCDF::NcType::ncType type;
  std::vector<StateVariables::StateVariableDimensionId> dimensions;
};

#endif // NETCDF_OUTPUT_AVAILABLE

class StateHeader {
public:
  StateHeader(int year, int month, int day, int nLayer, int nNode) : year(year), month(month), day(day), nLayer(nLayer), nNode(nNode) {}
  int year;
  int month;
  int day;
  int nLayer;
  int nNode;
  bool isValid() { return (year > 0 && month > 0 && day > 0); }
};

class StateIO {
public:
  enum IOType { Reader, Writer };
  StateIO(std::string filename, IOType type, const ProgramState* state);
  virtual ~StateIO();
  int process(int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(const int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(const double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(const char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);

  virtual void initializeOutput() = 0;
  virtual int write(const int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int write(const double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int write(const char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int processNewline();
  virtual StateHeader readHeader() = 0;
  virtual void notifyDimensionUpdate(StateVariables::StateVariableDimensionId dimension, int value = -1) {}
  virtual void initializeDimensionIndices() {}
  virtual int getCurrentDimensionIndex(StateVariables::StateVariableDimensionId dimension) { return -1; }
  virtual int seekToCell(int cellid, int* nVeg, int* nBand) = 0;
  virtual int read(int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int read(double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int read(char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual void flush() = 0;
  virtual void rewindFile() = 0;
  IOType getType() { return ioType; }
protected:
  std::string filename;
  const ProgramState* state;
  const IOType ioType;
};

#endif /* STATEIO_H_ */
