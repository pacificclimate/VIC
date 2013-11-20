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
void verifyGlobalAttributes(const NcFile& file);

extern const char* version; // Defined in global.h

const std::string stateYear = "state_year";
const std::string stateMonth = "state_month";
const std::string stateDay = "state_day";
const std::string stateNLayer = "state_nlayer";
const std::string stateNNode = "state_nnode";

StateIONetCDF::StateIONetCDF(std::string filename, IOType ioType, const ProgramState* state) : StateIO(filename, ioType, state), netCDF(NULL), filename(filename) {
  populateMetaData();
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

  NcFile::FileMode netcdfFileMode = NcFile::write;
  if (ioType == StateIO::Reader) {
    netcdfFileMode = NcFile::read;
  }
  // Do not specify the type of file that will be read here (NcFile::nc4). The library will throw an exception
  // if it is provided for open types read or write (since the file was already created with a certain format).
  netCDF = new NcFile(filename, netcdfFileMode);
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

  verifyGlobalAttributes(*netCDF);    //TODO: const this parameter.

  // Set up the dimensions and variables.

  fprintf(stderr, "Setting up grid dimensions, lat size: %ld, lon size: %ld\n",
      (size_t) state->global_param.gridNumLatDivisions,
      (size_t) state->global_param.gridNumLonDivisions);

  NcDim latDim = netCDF->addDim("lat", (size_t) state->global_param.gridNumLatDivisions);
  NcDim lonDim = netCDF->addDim("lon", (size_t) state->global_param.gridNumLonDivisions);
  NcDim bounds = netCDF->addDim("bnds", 2);

  NcDim layers = netCDF->addDim("layers", state->options.Nlayer);
  NcDim nodes = netCDF->addDim("nodes", state->options.Nnode);
  NcDim lakeNodes = netCDF->addDim("lake_active_nodes", MAX_LAKE_NODES+1);
  NcDim frostLayerAreas = netCDF->addDim("frost_layer_subareas", state->options.Nlayer * FROST_SUBAREAS);
  NcDim frostAreas = netCDF->addDim("frost_subareas", FROST_SUBAREAS);


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
  for (std::map<int, StateVariableMetaData>::iterator it = metaData.begin(); it != metaData.end(); ++it) {

    const std::string varName = it->second.name;
    int id = it->first;

    std::vector<NcDim> dimensions;

    if (id == SOIL_DZ_NODE || id == SOIL_ZSUM_NODE || id == ENERGY_T || id == LAKE_ENERGY_T) {
      const NcDim tempDims [] = { latDim, lonDim, nodes };
      dimensions.assign(tempDims, tempDims + 3);
    } else {
      const NcDim tempDims [] = { latDim, lonDim };
      dimensions.assign(tempDims, tempDims + 2);
    }

    try {
      NcVar data = netCDF->addVar(varName, ncFloat, dimensions);
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

int StateIONetCDF::write(const int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  return 0;
}

int StateIONetCDF::write(const double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  return 0;
}

int StateIONetCDF::write(const char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  return 0;
}

int StateIONetCDF::read(int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  return 0;
}

int StateIONetCDF::read(double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  return 0;
}

int StateIONetCDF::read(char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) {
  return 0;
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

int StateIONetCDF::seekToCell(int cellid, int* nVeg, int* nBand) {
  return 0;
}

void StateIONetCDF::flush() {
  //TODO: implement if applicable
}

void StateIONetCDF::rewindFile() {

}

void StateIONetCDF::populateMetaData() {
  using namespace StateVariables;
  metaData[NONE] = StateVariableMetaData("NONE");
  metaData[GRID_CELL] = StateVariableMetaData("GRID_CELL");
  metaData[VEG_TYPE_NUM] = StateVariableMetaData("VEG_TYPE_NUM");
  metaData[NUM_BANDS] = StateVariableMetaData("NUM_BANDS");
  metaData[SOIL_DZ_NODE] = StateVariableMetaData("SOIL_DZ_NODE");
  metaData[SOIL_ZSUM_NODE] = StateVariableMetaData("SOIL_ZSUM_NODE");
  metaData[SOIL_DEPTH] = StateVariableMetaData("SOIL_DEPTH");
  metaData[SOIL_EFFECTIVE_POROSITY] = StateVariableMetaData("SOIL_EFFECTIVE_POROSITY");
  metaData[SOIL_DP] = StateVariableMetaData("SOIL_DP");
  metaData[PRCP_MU] = StateVariableMetaData("PRCP_MU");
  metaData[INIT_STILL_STORM] = StateVariableMetaData("INIT_STILL_STORM");
  metaData[INIT_DRY_TIME] =   StateVariableMetaData("INIT_DRY_TIME");
  metaData[HRU_VEG_INDEX] = StateVariableMetaData("HRU_VEG_INDEX");
  metaData[HRU_BAND_INDEX] = StateVariableMetaData("HRU_BAND_INDEX");
  metaData[LAYER_MOIST] = StateVariableMetaData("LAYER_MOIST");
  metaData[LAYER_SOIL_ICE] = StateVariableMetaData("LAYER_SOIL_ICE");
  metaData[LAYER_ICE_CONTENT] = StateVariableMetaData("LAYER_ICE_CONTENT");
  metaData[HRU_VEG_VAR_WDEW] = StateVariableMetaData("HRU_VEG_VAR_WDEW");
  metaData[SNOW_LAST_SNOW] = StateVariableMetaData("SNOW_LAST_SNOW");
  metaData[SNOW_MELTING] = StateVariableMetaData("SNOW_MELTING");
  metaData[SNOW_COVERAGE] = StateVariableMetaData("SNOW_COVERAGE");
  metaData[SNOW_SWQ] = StateVariableMetaData("SNOW_SWQ");
  metaData[SNOW_SURF_TEMP] = StateVariableMetaData("SNOW_SURF_TEMP");
  metaData[SNOW_SURF_WATER] = StateVariableMetaData("SNOW_SURF_WATER");
  metaData[SNOW_PACK_TEMP] = StateVariableMetaData("SNOW_PACK_TEMP");
  metaData[SNOW_PACK_WATER] = StateVariableMetaData("SNOW_PACK_WATER");
  metaData[SNOW_DENSITY] = StateVariableMetaData("SNOW_DENSITY");
  metaData[SNOW_COLD_CONTENT] = StateVariableMetaData("SNOW_COLD_CONTENT");
  metaData[SNOW_CANOPY] = StateVariableMetaData("SNOW_CANOPY");
  metaData[ENERGY_T] = StateVariableMetaData("ENERGY_T");
  metaData[LAKE_LAYER_MOIST] = StateVariableMetaData("LAKE_LAYER_MOIST");
  metaData[LAKE_LAYER_SOIL_ICE] = StateVariableMetaData("LAKE_LAYER_SOIL_ICE");
  metaData[LAKE_LAYER_ICE_CONTENT] = StateVariableMetaData("LAKE_LAYER_ICE_CONTENT");
  metaData[LAKE_SNOW_LAST_SNOW] = StateVariableMetaData("LAKE_SNOW_LAST_SNOW");
  metaData[LAKE_SNOW_MELTING] = StateVariableMetaData("LAKE_SNOW_MELTING");
  metaData[LAKE_SNOW_COVERAGE] = StateVariableMetaData("LAKE_SNOW_COVERAGE");
  metaData[LAKE_SNOW_SWQ] = StateVariableMetaData("LAKE_SNOW_SWQ");
  metaData[LAKE_SNOW_SURF_TEMP] = StateVariableMetaData("LAKE_SNOW_SURF_TEMP");
  metaData[LAKE_SNOW_SURF_WATER] = StateVariableMetaData("LAKE_SNOW_SURF_WATER");
  metaData[LAKE_SNOW_PACK_TEMP] = StateVariableMetaData("LAKE_SNOW_PACK_TEMP");
  metaData[LAKE_SNOW_PACK_WATER] = StateVariableMetaData("LAKE_SNOW_PACK_WATER");
  metaData[LAKE_SNOW_DENSITY] = StateVariableMetaData("LAKE_SNOW_DENSITY");
  metaData[LAKE_SNOW_COLD_CONTENT] = StateVariableMetaData("LAKE_SNOW_COLD_CONTENT");
  metaData[LAKE_SNOW_CANOPY] = StateVariableMetaData("LAKE_SNOW_CANOPY");
  metaData[LAKE_ENERGY_T] = StateVariableMetaData("LAKE_ENERGY_T");
  metaData[LAKE_ACTIVENOD] = StateVariableMetaData("LAKE_ACTIVENOD");
  metaData[LAKE_DZ] = StateVariableMetaData("LAKE_DZ");
  metaData[LAKE_SURFDZ] = StateVariableMetaData("LAKE_SURFDZ");
  metaData[LAKE_LDEPTH] = StateVariableMetaData("LAKE_LDEPTH");
  metaData[LAKE_SURFACE] = StateVariableMetaData("LAKE_SURFACE");
  metaData[LAKE_SAREA] = StateVariableMetaData("LAKE_SAREA");
  metaData[LAKE_VOLUME] = StateVariableMetaData("LAKE_VOLUME");
  metaData[LAKE_TEMP] = StateVariableMetaData("LAKE_TEMP");
  metaData[LAKE_TEMPAVG] = StateVariableMetaData("LAKE_TEMPAVG");
  metaData[LAKE_AREAI] = StateVariableMetaData("LAKE_AREAI");
  metaData[LAKE_NEW_ICE_AREA] = StateVariableMetaData("LAKE_NEW_ICE_AREA");
  metaData[LAKE_ICE_WATER_EQ] = StateVariableMetaData("LAKE_ICE_WATER_EQ");
  metaData[LAKE_HICE] = StateVariableMetaData("LAKE_HICE");
  metaData[LAKE_TEMPI] = StateVariableMetaData("LAKE_TEMPI");
  metaData[LAKE_SWE] = StateVariableMetaData("LAKE_SWE");
  metaData[LAKE_SURF_TEMP] = StateVariableMetaData("LAKE_SURF_TEMP");
  metaData[LAKE_PACK_TEMP] = StateVariableMetaData("LAKE_PACK_TEMP");
  metaData[LAKE_SALBEDO] = StateVariableMetaData("LAKE_SALBEDO");
  metaData[LAKE_SDEPTH] = StateVariableMetaData("LAKE_SDEPTH");
}
