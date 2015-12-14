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
  LAKE_SDEPTH,


  STATE_num_veg_types, //new
  STATE_NR, //new
  STATE_NF, //new


  GLOBAL_PARAM_MEASURE_H, //new
  GLOBAL_PARAM_WIND_H, //new
  GLOBAL_PARAM_RESOLUTION, //new
  GLOBAL_PARAM_DT, //new
  GLOBAL_PARAM_OUT_DT, //new
  GLOBAL_PARAM_ENDDAY, //new
  GLOBAL_PARAM_ENDMONTH, //new
  GLOBAL_PARAM_ENDYEAR, //new
  GLOBAL_PARAM_FORCEDAY, //new
  GLOBAL_PARAM_FORCEHOUR, //new
  GLOBAL_PARAM_FORCEMONTH, //new
  GLOBAL_PARAM_FORCESKIP, //new
  GLOBAL_PARAM_FORCEYEAR, //new
  GLOBAL_PARAM_NRECS, //new
  GLOBAL_PARAM_SKIPYEAR, //new
  GLOBAL_PARAM_STARTDAY, //new
  GLOBAL_PARAM_STARTHOUR, //new
  GLOBAL_PARAM_STARTMONTH, //new
  GLOBAL_PARAM_STARTYEAR, //new
  GLOBAL_PARAM_STATEDAY, //new
  GLOBAL_PARAM_STATEMONTH, //new
  GLOBAL_PARAM_STATEYEAR, //new
  GLOBAL_PARAM_GLACIERACCUMSTARTYEAR, //new
  GLOBAL_PARAM_GLACIERACCUMSTARTMONTH, //new
  GLOBAL_PARAM_GLACIERACCUMSTARTDAY, //new
  GLOBAL_PARAM_GLACIERACCUMINTERVAL, //new
  GLOBAL_PARAM_GRIDSTARTLAT, //new
  GLOBAL_PARAM_GRIDSTARTLON, //new
  GLOBAL_PARAM_GRIDSTEPLAT, //new
  GLOBAL_PARAM_GRIDSTEPLON, //new
  GLOBAL_PARAM_GRIDNUMLATDIVISIONS, //new
  GLOBAL_PARAM_GRIDNUMLONDIVISIONS, //new

  VEG_LIB_OVERSTORY, //new
  VEG_LIB_LAI, //new
  VEG_LIB_WDMAX, //new
  VEG_LIB_ALBEDO, //new
  VEG_LIB_DISPLACEMENT, //new
  VEG_LIB_EMISSIVITY, //new
  VEG_LIB_NVEGLIBTYPES, //new
  VEG_LIB_RAD_ATTEN, //new
  VEG_LIB_RARC, //new
  VEG_LIB_RMIN, //new
  VEG_LIB_ROUGHNESS, //new
  VEG_LIB_TRUNK_RATIO, //new
  VEG_LIB_WIND_ATTEN, //new
  VEG_LIB_WIND_H, //new
  VEG_LIB_RGL, //new
  VEG_LIB_VEG_CLASS, //new


  OPTIONS_ABOVETREELINEVEG, //new
  OPTIONS_AERO_RESIST_CANSNOW, //new
  OPTIONS_BLOWING, //new
  OPTIONS_COMPUTE_TREELINE, //new
  OPTIONS_CONTINUEONERROR, //new
  OPTIONS_CORRPREC, //new
  OPTIONS_DIST_PRCP, //new
  OPTIONS_EQUAL_AREA, //new
  OPTIONS_EXP_TRANS, //new
  OPTIONS_FROZEN_SOIL, //new
  OPTIONS_FULL_ENERGY, //new
  OPTIONS_GRND_FLUX_TYPE, //new
  OPTIONS_IMPLICIT, //new
  OPTIONS_JULY_TAVG_SUPPLIED, //new
  OPTIONS_LAKES, //new
  OPTIONS_LW_CLOUD, //new
  OPTIONS_LW_TYPE, //new
  OPTIONS_MIN_WIND_SPEED, //new
  OPTIONS_MTCLIM_SWE_CORR, //new
  OPTIONS_NLAYER, //new
  OPTIONS_NNODE, //new
  OPTIONS_NOFLUX, //new
  OPTIONS_PLAPSE, //new
  OPTIONS_PREC_EXPT, //new
  OPTIONS_ROOT_ZONES, //new
  OPTIONS_QUICK_FLUX, //new
  OPTIONS_QUICK_SOLVE, //new
  OPTIONS_SNOW_ALBEDO, //new
  OPTIONS_SNOW_DENSITY, //new
  OPTIONS_SNOW_BAND, //new
  OPTIONS_SNOW_STEP, //new
  OPTIONS_SW_PREC_THRESH, //new
  OPTIONS_TFALLBACK, //new
  OPTIONS_VP_INTERP, //new
  OPTIONS_VP_ITER, //new
  OPTIONS_TEMP_TH_TYPE, //new MS, 19-10-15
  OPTIONS_NLAKENODE, //new
  OPTIONS_ALMA_INPUT, //new
  OPTIONS_ARC_SOIL, //new
  OPTIONS_BASEFLOW, //new
  OPTIONS_GRID_DECIMAL, //new
  OPTIONS_VEGPARAM_LAI, //new
  OPTIONS_LAI_SRC, //new
  OPTIONS_LAKE_PROFILE, //new
  OPTIONS_ORGANIC_FRACT, //new
  OPTIONS_INIT_STATE, //new
  OPTIONS_SAVE_STATE, //new
  OPTIONS_MAX_MEMORY, //new
  OPTIONS_ALMA_OUTPUT, //new
  OPTIONS_COMPRESS, //new
  OPTIONS_MOISTFRACT, //new
  OPTIONS_NOUTFILES, //new
  OPTIONS_PRT_HEADER, //new
  OPTIONS_PRT_SNOW_BAND, //new
  OPTIONS_NETCDF_FULL_FILE_PATH, //new
  OPTIONS_GLACIER_ID, //new


  PARAM_SET_FORCE_DT, //new
  PARAM_SET_FORCE_ENDIAN, //new
  PARAM_SET_FORCE_FORMAT, //new
  PARAM_SET_FORCE_INDEX, //new
  PARAM_SET_N_TYPES, //new

  FORCE_TYPE_SIGNED, //new
  FORCE_TYPE_SUPPLIED, //new
  FORCE_TYPE_multiplier, //new

  CELL_ERRSTR, //new
  CELL_CV_SUM, //new

  FALLBACK_Tfoliage_fbcount_total, //new
  FALLBACK_Tcanopy_fbcount_total, //new
  FALLBACK_Tsnowsurf_fbcount_total, //new
  FALLBACK_Tsurf_fbcount_total, //new
  FALLBACK_Tsoil_fbcount_total, //new
  FALLBACK_Tglacsurf_fbcount_total, //new

  CELLBALANCE_water_last_storage, //new
  CELLBALANCE_water_cum_error, //new
  CELLBALANCE_water_max_error, //new
  CELLBALANCE_energy_cum_error, //new
  CELLBALANCE_energy_max_error, //new

  SAVE_DATA_total_soil_moist, //new
  SAVE_DATA_surfstor, //new
  SAVE_DATA_swe, //new
  SAVE_DATA_wdew, //new

  ATMOS_AIR_TEMP, //new
  ATMOS_CHANNEL_IN, //new
  ATMOS_DENSITY, //new
  ATMOS_LONGWAVE, //new
  ATMOS_OUT_PREC, //new
  ATMOS_OUT_RAIN, //new
  ATMOS_OUT_SNOW, //new
  ATMOS_PREC, //new
  ATMOS_PRESSURE, //new
  ATMOS_SHORTWAVE, //new
  ATMOS_SNOWFLAG, //new
  ATMOS_TSKC, //new
  ATMOS_VP, //new
  ATMOS_VPD, //new
  ATMOS_WIND, //new

  SOIL_FS_ACTIVE, //new
  SOIL_DS, //new
  SOIL_DSMAX, //new
  SOIL_KSAT, //new
  SOIL_WCR, //new
  SOIL_WPWP, //new
  SOIL_WS, //new
  SOIL_DS_ORIG, //new
  SOIL_DSMAX_ORIG, //new
  SOIL_WS_ORIG, //new
  SOIL_ALPHA, //new
  SOIL_ANNUAL_PREC, //new
  SOIL_AVG_TEMP, //new
  SOIL_AVGJULYAIRTEMP, //new
  SOIL_B_INFILT, //new
  SOIL_BETA, //new
  SOIL_BUBBLE, //new
  SOIL_BUBBLE_NODE, //new
  SOIL_BULK_DENSITY, //new
  SOIL_BULK_DENS_MIN, //new
  SOIL_BULK_DENS_ORG, //new
  SOIL_C, //new
  SOIL_DEPTH_FULL_SNOW_COVER, //new
  SOIL_EXPT, //new
  SOIL_EXPT_NODE, //new
  SOIL_FROST_FRACT, //new
  SOIL_FROST_SLOPE, //new
  SOIL_GAMMA, //new
  SOIL_INIT_MOIST, //new
  SOIL_MAX_INFIL, //new
  SOIL_MAX_MOIST, //new
  SOIL_MAX_MOIST_NODE, //new
  SOIL_PHI_S, //new
  SOIL_POROSITY, //new
  SOIL_QUARTZ, //new
  SOIL_ORGANIC, //new
  SOIL_RESID_MOIST, //new
  SOIL_ROUGH, //new
  SOIL_SNOW_ROUGH, //new
  SOIL_SOIL_DENSITY, //new
  SOIL_SOIL_DENS_MIN, //new
  SOIL_SOIL_DENS_ORG, //new
  SOIL_BANDELEV, //new
  SOIL_AREAFRACT, //new
  SOIL_PFACTOR, //new
  SOIL_TFACTOR, //new
  SOIL_ABOVETREELINE, //new
  SOIL_UFWC_TABLE_LAYER, //new
  SOIL_UFWC_TABLE_NODE, //new
  SOIL_ELEVATION, //new
  SOIL_LAT, //new
  SOIL_LNG, //new
  SOIL_CELL_AREA, //new
  SOIL_TIME_ZONE_LNG, //new
  SOIL_LAYER_NODE_FRACT, //new
  SOIL_GRIDCEL, //new
  SOIL_ZWTVMOIST_ZWT, //new
  SOIL_ZWTVMOIST_MOIST, //new
  SOIL_SLOPE, //new
  SOIL_ASPECT, //new
  SOIL_EHORIZ, //new
  SOIL_WHORIZ, //new
  SOIL_MIN_DEPTH, //new
  SOIL_POROSITY_NODE, //new
  SOIL_EFFECTIVE_POROSITY_NODE, //new
  SOIL_WCR_FRACT, //new
  SOIL_WPWP_FRACT, //new
  SOIL_SUBSIDENCE, //new
  SOIL_NEW_SNOW_ALB, //new
  SOIL_SNOW_ALB_ACCUM_A, //new
  SOIL_SNOW_ALB_ACCUM_B, //new
  SOIL_SNOW_ALB_THAW_A, //new
  SOIL_SNOW_ALB_THAW_B, //new
  SOIL_MIN_RAIN_TEMP, //new
  SOIL_MAX_SNOW_TEMP, //new
  SOIL_PADJ, //new
  SOIL_T_LAPSE, //new
  SOIL_PGRAD, //new
  SOIL_GLAC_SURF_THICK, //new
  SOIL_GLAC_SURF_WE, //new
  SOIL_GLAC_KMIN, //new
  SOIL_GLAC_DK, //new
  SOIL_GLAC_A, //new
  SOIL_GLAC_ALBEDO, //new
  SOIL_GLAC_ROUGH, //new

  HRU_VEG_CON_CV, //new
  HRU_VEG_CON_ROOT, //new
  HRU_VEG_CON_ZONE_DEPTH, //new
  HRU_VEG_CON_ZONE_FRACT, //new
  HRU_VEG_CON_VEGINDEX, //new
//  HRU_VEG_CON_VEGCLASS, //new
  HRU_VEG_CON_SIGMA_SLOPE, //new
  HRU_VEG_CON_LAG_ONE, //new
  HRU_VEG_CON_FETCH, //new
  HRU_VEG_CON_LAKE, //new

  HRU_ISGLACIER, //new
  HRU_ISBARESOIL, //new

  HRUCELL_ASAT, //new
  HRUCELL_BASEFLOW, //new
  HRUCELL_INFLOW, //new
  HRUCELL_POT_EVAP, //new
  HRUCELL_RUNOFF, //new
  HRUCELL_ROOTMOIST, //new
  HRUCELL_WETNESS, //new
  HRUCELL_ZWT, //new
  HRUCELL_ZWT2, //new
  HRUCELL_ZWT3, //new

  AERO_SURFACE, //new
  AERO_OVERSTORY, //new

  LAYER_CS, //new
  LAYER_T, //new
  LAYER_EVAP, //new
  LAYER_KAPPA, //new
  LAYER_PHI, //new
  LAYER_ZWT, //new

  HRU_VEG_VAR_CANOPYEVAP, //new
  HRU_VEG_VAR_THROUGHFALL, //new

  SNOW_ALBEDO, //new
  SNOW_CANOPY_ALBEDO, //new
  SNOW_DEPTH, //new
  SNOW_MAX_SWQ, //new
  SNOW_SNOW, //new
  SNOW_STORE_COVERAGE, //new
  SNOW_STORE_SNOW, //new
  SNOW_STORE_SWQ, //new
  SNOW_SURF_TEMP_FBCOUNT, //new
  SNOW_SURF_TEMP_FBFLAG, //new
  SNOW_SWQ_SLOPE, //new
  SNOW_TMP_INT_STORAGE, //new
  SNOW_BLOWING_FLUX, //new
  SNOW_CANOPY_VAPOR_FLUX, //new
  SNOW_MASS_ERROR, //new
  SNOW_MELT, //new
  SNOW_QNET, //new
  SNOW_SURFACE_FLUX, //new
  SNOW_TRANSPORT, //new
  SNOW_VAPOR_FLUX, //new

  GLAC_COLD_CONTENT, //new
  GLAC_SURF_TEMP, //new
  GLAC_SURF_TEMP_FBCOUNT, //new
  GLAC_SURF_TEMP_FBFLAG, //new
  GLAC_QNET, //new
  GLAC_MASS_BALANCE, //new
  GLAC_ICE_MASS_BALANCE, //new
  GLAC_CUM_MASS_BALANCE, //new
  GLAC_ACCUMULATION, //new
  GLAC_MELT, //new
  GLAC_VAPOR_FLUX, //new
  GLAC_WATER_STORAGE, //new
  GLAC_OUTFLOW, //new
  GLAC_OUTFLOW_COEF, //new
  GLAC_INFLOW, //new

  ENERGY_ALBEDOLAKE, //new
  ENERGY_ALBEDOOVER, //new
  ENERGY_ALBEDOUNDER, //new
  ENERGY_CS, //new
  ENERGY_CS_NODE, //new
  ENERGY_FDEPTH, //new
  ENERGY_FROZEN, //new
  ENERGY_ICE_CONTENT, //new
  ENERGY_KAPPA, //new
  ENERGY_KAPPA_NODE, //new
  ENERGY_MOIST, //new
  ENERGY_NFROST, //new
  ENERGY_NTHAW, //new
  ENERGY_T_FBFLAG, //new
  ENERGY_T_FBCOUNT, //new
  ENERGY_T1_INDEX, //new
  ENERGY_TCANOPY, //new
  ENERGY_TCANOPY_FBFLAG, //new
  ENERGY_TCANOPY_FBCOUNT, //new
  ENERGY_TDEPTH, //new
  ENERGY_TFOLIAGE, //new
  ENERGY_TFOLIAGE_FBFLAG, //new
  ENERGY_TFOLIAGE_FBCOUNT, //new
  ENERGY_TSURF, //new
  ENERGY_TSURF_FBFLAG, //new
  ENERGY_TSURF_FBCOUNT, //new
  ENERGY_UNFROZEN, //new
  ENERGY_ADVECTED_SENSIBLE, //new
  ENERGY_ADVECTION, //new
  ENERGY_ATMOSERROR, //new
  ENERGY_ATMOSLATENT, //new
  ENERGY_ATMOSLATENTSUB, //new
  ENERGY_ATMOSSENSIBLE, //new
  ENERGY_CANOPY_ADVECTION, //new
  ENERGY_CANOPY_LATENT, //new
  ENERGY_CANOPY_LATENT_SUB, //new
  ENERGY_CANOPY_REFREEZE, //new
  ENERGY_CANOPY_SENSIBLE, //new
  ENERGY_DELTACC, //new
  ENERGY_DELTAH, //new
  ENERGY_ERROR, //new
  ENERGY_FUSION, //new
  ENERGY_GRND_FLUX, //new
  ENERGY_LATENT, //new
  ENERGY_LATENT_SUB, //new
  ENERGY_LONGWAVE, //new
  ENERGY_LONGOVERIN, //new
  ENERGY_LONGUNDERIN, //new
  ENERGY_LONGUNDEROUT, //new
  ENERGY_MELT_ENERGY, //new
  ENERGY_NETLONGATMOS, //new
  ENERGY_NETLONGOVER, //new
  ENERGY_NETLONGUNDER, //new
  ENERGY_NETSHORTATMOS, //new
  ENERGY_NETSHORTGRND, //new
  ENERGY_NETSHORTOVER, //new
  ENERGY_NETSHORTUNDER, //new
  ENERGY_OUT_LONG_CANOPY, //new
  ENERGY_OUT_LONG_SURFACE, //new
  ENERGY_REFREEZE_ENERGY, //new
  ENERGY_SENSIBLE, //new
  ENERGY_SHORTWAVE, //new
  ENERGY_SHORTOVERIN, //new
  ENERGY_SHORTUNDERIN, //new
  ENERGY_SNOW_FLUX, //new
  ENERGY_GLACIER_FLUX, //new
  ENERGY_DELTACC_GLAC, //new

	GLAC_MASS_BALANCE_INFO, // to store the Glacier Mass Balance information (grid cell num and equation) for use in VIC-RGM integration
	NUM_GLAC_MASS_BALANCE_INFO_TERMS,
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
  DIST_DIM,
	GLAC_MASS_BALANCE_INFO_DIM,
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
  int process(float* data, int numValues, const StateVariables::StateMetaDataVariableIndices id); //new
  int process(bool* data, int numValues, const StateVariables::StateMetaDataVariableIndices id); //new
  int process(char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(const int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(const double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(const float* data, int numValues, const StateVariables::StateMetaDataVariableIndices id); //new
  int process(const bool* data, int numValues, const StateVariables::StateMetaDataVariableIndices id); //new
  int process(const char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);

  virtual void initializeOutput() = 0;
  virtual int write(const int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int write(const double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int write(const float* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0; //new
  virtual int write(const bool* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0; //new
  virtual int write(const char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int processNewline();
  virtual StateHeader readHeader() = 0;
  virtual void notifyDimensionUpdate(StateVariables::StateVariableDimensionId dimension, int value = -1) {}
  virtual void initializeDimensionIndices() {}
  virtual int getCurrentDimensionIndex(StateVariables::StateVariableDimensionId dimension) { return -1; }
  virtual int seekToCell(int cellid, int* nVeg, int* nBand) = 0;
  virtual int read(int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int read(double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int read(float* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0; //new
  virtual int read(bool* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0; //new
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
