#include "StateIONetCDF.h"

#include <netcdf>
#include <sstream>

#include "vicNl_def.h"

using netCDF::NcFile;
using netCDF::NcDim;
using netCDF::NcVar;
using netCDF::ncFloat;
using netCDF::ncInt;

// The following functions are defined in the WriteOutputNetCDF.c file.
std::string getCompilationMachineInfo();
std::string getDateOfCompilation();
std::string getRuntimeMachineInfo();
std::string getSourceVersion();
int latitudeToIndex(double lat, const ProgramState* state);
int longitudeToIndex(double lon, const ProgramState* state);
void verifyGlobalAttributes(const NcFile& file);

extern const char* version; // Defined in global.h

const std::string stateYear = "state_year";
const std::string stateMonth = "state_month";
const std::string stateDay = "state_day";
const std::string stateNLayer = "state_nlayer";
const std::string stateNNode = "state_nnode";

StateIONetCDF::StateIONetCDF(std::string filename, IOType ioType, const ProgramState* state) : StateIO(filename, ioType, state), netCDF(NULL), filename(filename) {
  populateMetaData();
  populateMetaDimensions();
  initializeDimensionIndices();
  openFile();
}

StateIONetCDF::~StateIONetCDF() {
  closeFile();
}

void StateIONetCDF::closeFile() {
  // The file will be automatically close when the NcFile object has
  // its destructor called. This frees up any internal netCDF resources
  // associated with the file, and flushes any buffers.
  if (netCDF != NULL) {
    delete netCDF;
    netCDF = NULL;
  }
}

void StateIONetCDF::openFile() {

  closeFile();  // Safety check against memory leaks.

  if (ioType == StateIO::Reader) {
    // Do not specify the type of file that will be read here (NcFile::nc4). The library will throw an exception
    // if it is provided for open types read or write (since the file was already created with a certain format).
    try {
      netCDF = new NcFile(filename, NcFile::read);
    } catch (netCDF::exceptions::NcException& e) {
      fprintf(stderr, "Error: could not open input netCDF state file \"%s\" check that it exists and is actually a netCDF formatted file\n", filename.c_str());
      throw;
    }
  } else {
    try {
      netCDF = new NcFile(filename, NcFile::write);
    } catch (netCDF::exceptions::NcNotNCF& e) { // File doesn't exist or is not netCDF format.
      // Ignored. This is expected the first time an instance of this class is created without calling initializeOutput.
    }
  }

}

  // Initialize global attributes, setup variable dimensions and structure.
void StateIONetCDF::initializeOutput() {

  closeFile();  // Close the file that was opened in the constructor so we can recreate it properly.

  netCDF = new NcFile(filename.c_str(), NcFile::replace, NcFile::nc4);

  // Add global attributes here. (These could potentially be overwritten by inputs from the global file)
  netCDF->putAtt("title", "VIC output.");
  std::string source = std::string("VIC ") + version + ". ";
  source += "(Built from source: " + getSourceVersion() + ").";
  std::string history = "Created by " + source + ".\n";
  history += " Compiled on: " + getDateOfCompilation() + ".\n";
  history += " Compiled by machine: " + getCompilationMachineInfo() + ".\n";
  history += " Model run by machine: " + getRuntimeMachineInfo();

  // Add global attributes specified in the global input file.
  // Special cases for "source" and "history" attributes: they can't be overwritten, just appended to.
  for (std::vector<std::pair<std::string, std::string> >::const_iterator it =
      state->global_param.netCDFGlobalAttributes.begin();
      it != state->global_param.netCDFGlobalAttributes.end(); ++it) {
    if (it->first.compare("source") == 0) {
      source += it->second;
    } else if (it->first.compare("history") == 0) {
      history += it->second;
    } else {
      netCDF->putAtt(it->first, it->second);
    }
  }

  netCDF->putAtt("source", source.c_str());
  netCDF->putAtt("history", history.c_str());
  netCDF->putAtt("frequency", state->global_param.out_dt < 24 ? "hour" : "day");
  netCDF->putAtt("Conventions", "CF-1.6");

  /* Write save state date information */
  netCDF->putAtt(stateYear, netCDF::ncInt, state->global_param.stateyear);
  netCDF->putAtt(stateMonth, netCDF::ncInt, state->global_param.statemonth);
  netCDF->putAtt(stateDay, netCDF::ncInt, state->global_param.stateday);
  netCDF->putAtt(stateNLayer, netCDF::ncInt, state->options.Nlayer);
  netCDF->putAtt(stateNNode, netCDF::ncInt, state->options.Nnode);

  verifyGlobalAttributes(*netCDF);

  // Set up the dimensions and variables.

  fprintf(stderr, "Setting up grid dimensions, lat size: %ld, lon size: %ld\n",
      (size_t) state->global_param.gridNumLatDivisions,
      (size_t) state->global_param.gridNumLonDivisions);

  NcDim latDim = netCDF->addDim("lat", (size_t) state->global_param.gridNumLatDivisions);
  NcDim lonDim = netCDF->addDim("lon", (size_t) state->global_param.gridNumLonDivisions);
  NcDim boundsDim = netCDF->addDim("bnds", 2);

  std::map<StateVariables::StateVariableLastDimension, NcDim> allDynamicDimensions;
  for (std::map<StateVariables::StateVariableLastDimension, StateVariableDimension>::iterator it = metaDimensions.begin();
      it != metaDimensions.end(); ++it) {
    NcDim tempDim = netCDF->addDim(it->second.name, it->second.size);
    allDynamicDimensions[it->first] = tempDim;
  }

  // Define the coordinate variables.
  NcVar latVar = netCDF->addVar("lat", ncFloat, latDim);
  NcVar lonVar = netCDF->addVar("lon", ncFloat, lonDim);

  latVar.putAtt("axis", "Y");
  latVar.putAtt("units", "degrees_north");
  latVar.putAtt("standard name", "latitude");
  latVar.putAtt("long name", "latitude");
  latVar.putAtt("bounds", "lat_bnds");

  lonVar.putAtt("axis", "X");
  lonVar.putAtt("units", "degrees_east");
  lonVar.putAtt("standard name", "longitude");
  lonVar.putAtt("long name", "longitude");
  lonVar.putAtt("bounds", "lon_bnds");

  using namespace StateVariables;
  for (std::map<StateVariables::StateMetaDataVariableIndices, StateVariableMetaData>::iterator it = metaData.begin();
      it != metaData.end(); ++it) {

    const std::string varName = it->second.name;
    int id = it->first;

    std::vector<NcDim> dimensions;
    dimensions.push_back(latDim);
    dimensions.push_back(lonDim);

    for (std::vector<StateVariableLastDimension>::iterator dimIt =
        it->second.dimensions.begin(); dimIt != it->second.dimensions.end();
        ++dimIt) {
      if (*dimIt != StateVariables::NO_DIM) {
        dimensions.push_back(allDynamicDimensions[*dimIt]);
      }
    }

    try {
      NcVar data = netCDF->addVar(varName, netCDF::NcType(it->second.type), dimensions);
      data.putAtt("internal_id", ncInt, id); // This is basically just for reference, it might change between versions.
      if (state->options.COMPRESS) {
        data.setCompression(false, true, 1); // Some reasonable compression level - not too intensive.
      }
    } catch (const netCDF::exceptions::NcException& except) {
      fprintf(stderr, "Error adding variable: %s with id: %d Internal netCDF exception\n", varName.c_str(), id);
      throw;
    }
  }
}

// This template function means that we don't have to copy and paste the same code just for different types.
// The only type sensitive code in this function is the variable.putVar() call and this is overloaded for each type.
template<typename T> int StateIONetCDF::generalWrite(const T* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  std::vector<size_t> start;
  std::vector<size_t> count;
  start.push_back(curLatIndex);
  start.push_back(curLonIndex);
  count.push_back(1);
  count.push_back(1);
  StateVariables::StateVariableLastDimension lastDimensionId = StateVariables::NO_DIM;
  try {
    for (std::vector<StateVariables::StateVariableLastDimension>::iterator it = metaData[id].dimensions.begin();
        it != metaData[id].dimensions.end(); ++it) {
      if (*it != StateVariables::NO_DIM) {
        lastDimensionId = *it;
        start.push_back(curDimensionIndices[*it]);
        count.push_back(1);
      }
    }
    count[count.size() - 1] = numValues;  // Assumes that the list of values will go to the last dimension.

    NcVar variable = netCDF->getVar(metaData[id].name);
    variable.putVar(start, count, data);

  } catch (std::exception& e) {
    fprintf(stderr, "Error writing variable: %s, at latIndex: %d, lonIndex: %d. numValues = %d, last dimensionId = %d, last dimension length = %d\n",
        metaData[id].name.c_str(), curLatIndex, curLonIndex, numValues,
        lastDimensionId, metaDimensions[lastDimensionId].size);
    throw;
  }
  return numValues;
}

// This template function means that we don't have to copy and paste the same code just for different types.
// The only type sensitive code in this function is the variable.getVar() call and this is overloaded for each type.
template<typename T> int StateIONetCDF::generalRead(T* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  std::vector<size_t> start;
  std::vector<size_t> count;
  start.push_back(curLatIndex);
  start.push_back(curLonIndex);
  count.push_back(1);
  count.push_back(1);
  StateVariables::StateVariableLastDimension lastDimensionId = StateVariables::NO_DIM;
  try {
    for (std::vector<StateVariables::StateVariableLastDimension>::iterator it = metaData[id].dimensions.begin();
        it != metaData[id].dimensions.end(); ++it) {
      if (*it != StateVariables::NO_DIM) {
        lastDimensionId = *it;
        start.push_back(curDimensionIndices[*it]);
        count.push_back(1);
      }
    }
    count[count.size() - 1] = numValues;  // Assumes that the list of values will go to the last dimension.

    NcVar variable = netCDF->getVar(metaData[id].name);
    variable.getVar(start, count, data);

  } catch (std::exception& e) {
    fprintf(stderr, "Error writing variable: %s, at latIndex: %d, lonIndex: %d. numValues = %d, last dimensionId = %d, last dimension length = %d\n",
        metaData[id].name.c_str(), curLatIndex, curLonIndex, numValues,
        lastDimensionId, metaDimensions[lastDimensionId].size);
    throw;
  }
  return numValues;
}

int StateIONetCDF::write(const int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  return generalWrite(data, numValues, id);
}

int StateIONetCDF::write(const double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  return generalWrite(data, numValues, id);
}

int StateIONetCDF::write(const char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  return generalWrite(data, numValues, id);
}

int StateIONetCDF::read(int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  return generalRead(data, numValues, id);
}

int StateIONetCDF::read(double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  return generalRead(data, numValues, id);
}

int StateIONetCDF::read(char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  return generalRead(data, numValues, id);
}

void StateIONetCDF::checkAndRead(std::string name, int * storage, int length) {
  netCDF::NcGroupAtt att = netCDF->getAtt(name);
  if (att.getAttLength() == (unsigned int)length) {
    att.getValues(storage);
  } else {
    std::stringstream ss;
    ss << "Error: mismatch in variable sizes when reading " << name << " expected=" << length << " actual=" << att.getAttLength();
    throw new VICException(ss.str());
  }
}

StateHeader StateIONetCDF::readHeader() {
  int year, month, day, nLayer, nNode;
  checkAndRead(stateYear, &year, 1);
  checkAndRead(stateMonth, &month, 1);
  checkAndRead(stateDay, &day, 1);
  checkAndRead(stateNLayer, &nLayer, 1);
  checkAndRead(stateNNode, &nNode, 1);
  return StateHeader(year, month, day, nLayer, nNode);
}

void StateIONetCDF::notifyCellLocation(float lat, float lng) {
  curLatIndex = latitudeToIndex(lat, state);
  curLonIndex = longitudeToIndex(lng, state);
  initializeDimensionIndices(); // This is reset for each cell.
}

void StateIONetCDF::notifyDimensionUpdate(StateVariables::StateVariableLastDimension dimension, int value) {
  if (curDimensionIndices.find(dimension) == curDimensionIndices.end() || metaDimensions.find(dimension) == metaDimensions.end()) {
    std::stringstream ss;
    ss << "Error: could not find dimension " << dimension << ". Make sure it has an entry in both maps: curDimensionIndices and metaDimensions.";
    throw new VICException(ss.str());
  }
  int newValue = value;
  if (newValue < 0) {
    newValue = curDimensionIndices[dimension] + 1;
  }
  if (newValue < 0 || newValue >= metaDimensions[dimension].size) {
    std::stringstream ss;
    ss << "Error: the dimension " << dimension << " cannot have an index value of " << newValue;
    throw new VICException(ss.str());
  }
  curDimensionIndices[dimension] = newValue;
}

int StateIONetCDF::seekToCell(int cellid, int* nVeg, int* nBand) {
  NcDim latDim = netCDF->getDim("lat");
  NcDim lonDim = netCDF->getDim("lon");
  int latSize = latDim.getSize();
  int lonSize = lonDim.getSize();
  NcVar cellIds = netCDF->getVar("GRID_CELL");    //TODO: reuse string
  for (int i = 0; i < latSize; i++) {
    for (int j = 0; j < lonSize; j++) {
      double cellIdRead = -1;
      std::vector<size_t> start;
      start.push_back((size_t)i);
      start.push_back((size_t)j);
      cellIds.getVar(start, &cellIdRead);
      if ((int)cellIdRead == cellid) {
        NcVar veg = netCDF->getVar("VEG_TYPE_NUM");
        NcVar band = netCDF->getVar("NUM_BANDS"); //TODO: string reuse
        veg.getVar(start, nVeg);
        band.getVar(start, &nBand);
        return 0;
      }
    }
  }
  return -1;
}

void StateIONetCDF::flush() {
  //TODO: implement if applicable
}

void StateIONetCDF::rewindFile() {

}

// Reset all dimension indices to zero.
void StateIONetCDF::initializeDimensionIndices() {
  for (std::map<StateVariables::StateVariableLastDimension, StateVariableDimension>::iterator it = metaDimensions.begin();
      it != metaDimensions.end(); ++it) {
    curDimensionIndices[it->first] = 0;  // This will change to 0 on the first call to notifyDimensionUpdate.
  }
}

void StateIONetCDF::populateMetaDimensions() {
  using namespace StateVariables;
  metaDimensions[NO_DIM] = StateVariableDimension("no_dim", 1);
  metaDimensions[LAYERS_DIM] = StateVariableDimension("Nlayers", state->options.Nlayer);
  metaDimensions[NODES_DIM] = StateVariableDimension("Nnodes", state->options.Nnode);
  metaDimensions[LAKE_NODES_DIM] = StateVariableDimension("lake_active_nodes", MAX_LAKE_NODES+1);
  metaDimensions[FROST_LAYER_AREAS_DIM] = StateVariableDimension("frost_layer_subareas", state->options.Nlayer * FROST_SUBAREAS);
  metaDimensions[FROST_AREAS_DIM] = StateVariableDimension("frost_subareas", FROST_SUBAREAS);
  metaDimensions[HRU_DIM] = StateVariableDimension("hru", state->veg_lib->NVegLibTypes * MAX_BANDS);
  metaDimensions[VEG_DIM] = StateVariableDimension("vegetation", state->veg_lib->NVegLibTypes);
  metaDimensions[DIST_DIM] = StateVariableDimension("dist", 2); // Wet and dry.
}

void StateIONetCDF::populateMetaData() {
  using namespace StateVariables;
  metaData[NONE] =                    StateVariableMetaData("NONE");
  metaData[GRID_CELL] =               StateVariableMetaData("GRID_CELL");
  metaData[VEG_TYPE_NUM] =            StateVariableMetaData("VEG_TYPE_NUM");
  metaData[NUM_BANDS] =               StateVariableMetaData("NUM_BANDS");
  metaData[SOIL_DZ_NODE] =            StateVariableMetaData("SOIL_DZ_NODE", NODES_DIM);
  metaData[SOIL_ZSUM_NODE] =          StateVariableMetaData("SOIL_ZSUM_NODE", NODES_DIM);
  metaData[SOIL_DEPTH] =              StateVariableMetaData("SOIL_DEPTH", LAYERS_DIM);
  metaData[SOIL_EFFECTIVE_POROSITY] = StateVariableMetaData("SOIL_EFFECTIVE_POROSITY", LAYERS_DIM);
  metaData[SOIL_DP] =                 StateVariableMetaData("SOIL_DP");
  metaData[PRCP_MU] =                 StateVariableMetaData("PRCP_MU", VEG_DIM);
  metaData[INIT_STILL_STORM] =        StateVariableMetaData("INIT_STILL_STORM", VEG_DIM);
  metaData[INIT_DRY_TIME] =           StateVariableMetaData("INIT_DRY_TIME", VEG_DIM);
  metaData[HRU_VEG_INDEX] =           StateVariableMetaData("HRU_VEG_INDEX", HRU_DIM);
  metaData[HRU_BAND_INDEX] =          StateVariableMetaData("HRU_BAND_INDEX", HRU_DIM);
  metaData[LAYER_MOIST] =             StateVariableMetaData("LAYER_MOIST", HRU_DIM, DIST_DIM, LAYERS_DIM);
  metaData[LAYER_SOIL_ICE] =          StateVariableMetaData("LAYER_SOIL_ICE", HRU_DIM, DIST_DIM, FROST_LAYER_AREAS_DIM);
  metaData[LAYER_ICE_CONTENT] =       StateVariableMetaData("LAYER_ICE_CONTENT", HRU_DIM, DIST_DIM, LAYERS_DIM);
  metaData[HRU_VEG_VAR_WDEW] =        StateVariableMetaData("HRU_VEG_VAR_WDEW", HRU_DIM, DIST_DIM);
  metaData[SNOW_LAST_SNOW] =          StateVariableMetaData("SNOW_LAST_SNOW", HRU_DIM);
  metaData[SNOW_MELTING] =            StateVariableMetaData("SNOW_MELTING", HRU_DIM);
  metaData[SNOW_COVERAGE] =           StateVariableMetaData("SNOW_COVERAGE", HRU_DIM);
  metaData[SNOW_SWQ] =                StateVariableMetaData("SNOW_SWQ", HRU_DIM);
  metaData[SNOW_SURF_TEMP] =          StateVariableMetaData("SNOW_SURF_TEMP", HRU_DIM);
  metaData[SNOW_SURF_WATER] =         StateVariableMetaData("SNOW_SURF_WATER", HRU_DIM);
  metaData[SNOW_PACK_TEMP] =          StateVariableMetaData("SNOW_PACK_TEMP", HRU_DIM);
  metaData[SNOW_PACK_WATER] =         StateVariableMetaData("SNOW_PACK_WATER", HRU_DIM);
  metaData[SNOW_DENSITY] =            StateVariableMetaData("SNOW_DENSITY"), HRU_DIM;
  metaData[SNOW_COLD_CONTENT] =       StateVariableMetaData("SNOW_COLD_CONTENT", HRU_DIM);
  metaData[SNOW_CANOPY] =             StateVariableMetaData("SNOW_CANOPY", HRU_DIM);
  metaData[ENERGY_T] =                StateVariableMetaData("ENERGY_T", HRU_DIM, NODES_DIM);
  metaData[LAKE_LAYER_MOIST] =        StateVariableMetaData("LAKE_LAYER_MOIST", DIST_DIM, LAYERS_DIM);
  metaData[LAKE_LAYER_SOIL_ICE] =     StateVariableMetaData("LAKE_LAYER_SOIL_ICE", DIST_DIM, FROST_LAYER_AREAS_DIM);
  metaData[LAKE_LAYER_ICE_CONTENT] =  StateVariableMetaData("LAKE_LAYER_ICE_CONTENT", DIST_DIM, LAYERS_DIM);
  metaData[LAKE_SNOW_LAST_SNOW] =     StateVariableMetaData("LAKE_SNOW_LAST_SNOW");
  metaData[LAKE_SNOW_MELTING] =       StateVariableMetaData("LAKE_SNOW_MELTING");
  metaData[LAKE_SNOW_COVERAGE] =      StateVariableMetaData("LAKE_SNOW_COVERAGE");
  metaData[LAKE_SNOW_SWQ] =           StateVariableMetaData("LAKE_SNOW_SWQ");
  metaData[LAKE_SNOW_SURF_TEMP] =     StateVariableMetaData("LAKE_SNOW_SURF_TEMP");
  metaData[LAKE_SNOW_SURF_WATER] =    StateVariableMetaData("LAKE_SNOW_SURF_WATER");
  metaData[LAKE_SNOW_PACK_TEMP] =     StateVariableMetaData("LAKE_SNOW_PACK_TEMP");
  metaData[LAKE_SNOW_PACK_WATER] =    StateVariableMetaData("LAKE_SNOW_PACK_WATER");
  metaData[LAKE_SNOW_DENSITY] =       StateVariableMetaData("LAKE_SNOW_DENSITY");
  metaData[LAKE_SNOW_COLD_CONTENT] =  StateVariableMetaData("LAKE_SNOW_COLD_CONTENT");
  metaData[LAKE_SNOW_CANOPY] =        StateVariableMetaData("LAKE_SNOW_CANOPY");
  metaData[LAKE_ENERGY_T] =           StateVariableMetaData("LAKE_ENERGY_T", LAKE_NODES_DIM);
  metaData[LAKE_ACTIVENOD] =          StateVariableMetaData("LAKE_ACTIVENOD");
  metaData[LAKE_DZ] =                 StateVariableMetaData("LAKE_DZ");
  metaData[LAKE_SURFDZ] =             StateVariableMetaData("LAKE_SURFDZ");
  metaData[LAKE_LDEPTH] =             StateVariableMetaData("LAKE_LDEPTH");
  metaData[LAKE_SURFACE] =            StateVariableMetaData("LAKE_SURFACE", LAKE_NODES_DIM);
  metaData[LAKE_SAREA] =              StateVariableMetaData("LAKE_SAREA");
  metaData[LAKE_VOLUME] =             StateVariableMetaData("LAKE_VOLUME");
  metaData[LAKE_TEMP] =               StateVariableMetaData("LAKE_TEMP", LAKE_NODES_DIM);
  metaData[LAKE_TEMPAVG] =            StateVariableMetaData("LAKE_TEMPAVG");
  metaData[LAKE_AREAI] =              StateVariableMetaData("LAKE_AREAI");
  metaData[LAKE_NEW_ICE_AREA] =       StateVariableMetaData("LAKE_NEW_ICE_AREA");
  metaData[LAKE_ICE_WATER_EQ] =       StateVariableMetaData("LAKE_ICE_WATER_EQ");
  metaData[LAKE_HICE] =               StateVariableMetaData("LAKE_HICE");
  metaData[LAKE_TEMPI] =              StateVariableMetaData("LAKE_TEMPI");
  metaData[LAKE_SWE] =                StateVariableMetaData("LAKE_SWE");
  metaData[LAKE_SURF_TEMP] =          StateVariableMetaData("LAKE_SURF_TEMP");
  metaData[LAKE_PACK_TEMP] =          StateVariableMetaData("LAKE_PACK_TEMP");
  metaData[LAKE_SALBEDO] =            StateVariableMetaData("LAKE_SALBEDO");
  metaData[LAKE_SDEPTH] =             StateVariableMetaData("LAKE_SDEPTH");

  // Make type adjustments if required.
  metaData[INIT_STILL_STORM].type = netCDF::NcType::nc_CHAR;
  metaData[SNOW_MELTING].type = netCDF::NcType::nc_CHAR;
  metaData[LAKE_SNOW_MELTING].type = netCDF::NcType::nc_CHAR;
}
