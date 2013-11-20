#ifndef STATEIO_H_
#define STATEIO_H_

#include "vicNl_def.h"
#include <string>

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
}

using std::string;
// This is just a wrapper class for now in case more attributes are needed in netCDF state files.
class StateVariableMetaData {
public:
  StateVariableMetaData() {}
  StateVariableMetaData(string name) : name(name) {}
  string name;
};

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
