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

  STATE_num_veg_types,
  STATE_NR,
  STATE_NF,

  GLOBAL_PARAM_MEASURE_H,
  GLOBAL_PARAM_WIND_H,
  GLOBAL_PARAM_RESOLUTION,
  GLOBAL_PARAM_DT,
  GLOBAL_PARAM_OUT_DT,
  GLOBAL_PARAM_ENDDAY,
  GLOBAL_PARAM_ENDMONTH,
  GLOBAL_PARAM_ENDYEAR,
  GLOBAL_PARAM_FORCEDAY,
  GLOBAL_PARAM_FORCEHOUR,
  GLOBAL_PARAM_FORCEMONTH,
  GLOBAL_PARAM_FORCESKIP,
  GLOBAL_PARAM_FORCEYEAR,
  GLOBAL_PARAM_NRECS,
  GLOBAL_PARAM_SKIPYEAR,
  GLOBAL_PARAM_STARTDAY,
  GLOBAL_PARAM_STARTHOUR,
  GLOBAL_PARAM_STARTMONTH,
  GLOBAL_PARAM_STARTYEAR,
  GLOBAL_PARAM_STATEDAY,
  GLOBAL_PARAM_STATEMONTH,
  GLOBAL_PARAM_STATEYEAR,
  GLOBAL_PARAM_GLACIERACCUMSTARTYEAR,
  GLOBAL_PARAM_GLACIERACCUMSTARTMONTH,
  GLOBAL_PARAM_GLACIERACCUMSTARTDAY,
  GLOBAL_PARAM_GLACIERACCUMINTERVAL,
  GLOBAL_PARAM_GRIDSTARTLAT,
  GLOBAL_PARAM_GRIDSTARTLON,
  GLOBAL_PARAM_GRIDSTEPLAT,
  GLOBAL_PARAM_GRIDSTEPLON,
  GLOBAL_PARAM_GRIDNUMLATDIVISIONS,
  GLOBAL_PARAM_GRIDNUMLONDIVISIONS,

  VEG_LIB_OVERSTORY,
  VEG_LIB_LAI,
  VEG_LIB_WDMAX,
  VEG_LIB_ALBEDO,
  VEG_LIB_DISPLACEMENT,
  VEG_LIB_EMISSIVITY,
  VEG_LIB_NVEGLIBTYPES,
  VEG_LIB_RAD_ATTEN,
  VEG_LIB_RARC,
  VEG_LIB_RMIN,
  VEG_LIB_ROUGHNESS,
  VEG_LIB_TRUNK_RATIO,
  VEG_LIB_WIND_ATTEN,
  VEG_LIB_WIND_H,
  VEG_LIB_RGL,
  VEG_LIB_VEG_CLASS,

  OPTIONS_ABOVETREELINEVEG,
  OPTIONS_AERO_RESIST_CANSNOW,
  OPTIONS_BLOWING,
  OPTIONS_COMPUTE_TREELINE,
  OPTIONS_CONTINUEONERROR,
  OPTIONS_CORRPREC,
  OPTIONS_DIST_PRCP,
  OPTIONS_EQUAL_AREA,
  OPTIONS_EXP_TRANS,
  OPTIONS_FROZEN_SOIL,
  OPTIONS_FULL_ENERGY,
  OPTIONS_GRND_FLUX_TYPE,
  OPTIONS_IMPLICIT,
  OPTIONS_JULY_TAVG_SUPPLIED,
  OPTIONS_LAKES,
  OPTIONS_LW_CLOUD,
  OPTIONS_LW_TYPE,
  OPTIONS_MIN_WIND_SPEED,
  OPTIONS_MTCLIM_SWE_CORR,
  OPTIONS_NLAYER,
  OPTIONS_NNODE,
  OPTIONS_NOFLUX,
  OPTIONS_PLAPSE,
  OPTIONS_PREC_EXPT,
  OPTIONS_ROOT_ZONES,
  OPTIONS_QUICK_FLUX,
  OPTIONS_QUICK_SOLVE,
  OPTIONS_SNOW_ALBEDO,
  OPTIONS_SNOW_DENSITY,
  OPTIONS_SNOW_BAND,
  OPTIONS_SNOW_STEP,
  OPTIONS_SW_PREC_THRESH,
  OPTIONS_TFALLBACK,
  OPTIONS_VP_INTERP,
  OPTIONS_VP_ITER,
  OPTIONS_TEMP_TH_TYPE,
  OPTIONS_NLAKENODE,
  OPTIONS_ALMA_INPUT,
  OPTIONS_ARC_SOIL,
  OPTIONS_BASEFLOW,
  OPTIONS_GRID_DECIMAL,
  OPTIONS_VEGPARAM_LAI,
  OPTIONS_LAI_SRC,
  OPTIONS_LAKE_PROFILE,
  OPTIONS_ORGANIC_FRACT,
  OPTIONS_INIT_STATE,
  OPTIONS_SAVE_STATE,
  OPTIONS_MAX_MEMORY,
  OPTIONS_ALMA_OUTPUT,
  OPTIONS_COMPRESS,
  OPTIONS_MOISTFRACT,
  OPTIONS_NOUTFILES,
  OPTIONS_PRT_HEADER,
  OPTIONS_PRT_SNOW_BAND,
  OPTIONS_NETCDF_FULL_FILE_PATH,
  OPTIONS_GLACIER_ID,

  PARAM_SET_FORCE_DT,
  PARAM_SET_FORCE_ENDIAN,
  PARAM_SET_FORCE_FORMAT,
  PARAM_SET_FORCE_INDEX,
  PARAM_SET_N_TYPES,

  FORCE_TYPE_SIGNED,
  FORCE_TYPE_SUPPLIED,
  FORCE_TYPE_multiplier,

  CELL_ERRSTR,
  CELL_CV_SUM,

  FALLBACK_Tfoliage_fbcount_total,
  FALLBACK_Tcanopy_fbcount_total,
  FALLBACK_Tsnowsurf_fbcount_total,
  FALLBACK_Tsurf_fbcount_total,
  FALLBACK_Tsoil_fbcount_total,
  FALLBACK_Tglacsurf_fbcount_total,

  CELLBALANCE_water_last_storage,
  CELLBALANCE_water_cum_error,
  CELLBALANCE_water_max_error,
  CELLBALANCE_energy_cum_error,
  CELLBALANCE_energy_max_error,

  SAVE_DATA_total_soil_moist,
  SAVE_DATA_surfstor,
  SAVE_DATA_swe,
  SAVE_DATA_wdew,

  ATMOS_AIR_TEMP,
  ATMOS_CHANNEL_IN,
  ATMOS_DENSITY,
  ATMOS_LONGWAVE,
  ATMOS_OUT_PREC,
  ATMOS_OUT_RAIN,
  ATMOS_OUT_SNOW,
  ATMOS_PREC,
  ATMOS_PRESSURE,
  ATMOS_SHORTWAVE,
  ATMOS_SNOWFLAG,
  ATMOS_TSKC,
  ATMOS_VP,
  ATMOS_VPD,
  ATMOS_WIND,

  SOIL_FS_ACTIVE,
  SOIL_DS,
  SOIL_DSMAX,
  SOIL_KSAT,
  SOIL_WCR,
  SOIL_WPWP,
  SOIL_WS,
  SOIL_DS_ORIG,
  SOIL_DSMAX_ORIG,
  SOIL_WS_ORIG,
  SOIL_ALPHA,
  SOIL_ANNUAL_PREC,
  SOIL_AVG_TEMP,
  SOIL_AVGJULYAIRTEMP,
  SOIL_B_INFILT,
  SOIL_BETA,
  SOIL_BUBBLE,
  SOIL_BUBBLE_NODE,
  SOIL_BULK_DENSITY,
  SOIL_BULK_DENS_MIN,
  SOIL_BULK_DENS_ORG,
  SOIL_C,
  SOIL_DEPTH_FULL_SNOW_COVER,
  SOIL_EXPT,
  SOIL_EXPT_NODE,
  SOIL_FROST_FRACT,
  SOIL_FROST_SLOPE,
  SOIL_GAMMA,
  SOIL_INIT_MOIST,
  SOIL_MAX_INFIL,
  SOIL_MAX_MOIST,
  SOIL_MAX_MOIST_NODE,
  SOIL_PHI_S,
  SOIL_POROSITY,
  SOIL_QUARTZ,
  SOIL_ORGANIC,
  SOIL_RESID_MOIST,
  SOIL_ROUGH,
  SOIL_SNOW_ROUGH,
  SOIL_SOIL_DENSITY,
  SOIL_SOIL_DENS_MIN,
  SOIL_SOIL_DENS_ORG,
  SOIL_BANDELEV,
  SOIL_AREAFRACT,
  SOIL_PFACTOR,
  SOIL_TFACTOR,
  SOIL_ABOVETREELINE,
  SOIL_UFWC_TABLE_LAYER,
  SOIL_UFWC_TABLE_NODE,
  SOIL_ELEVATION,
  SOIL_LAT,
  SOIL_LNG,
  SOIL_CELL_AREA,
  SOIL_TIME_ZONE_LNG,
  SOIL_LAYER_NODE_FRACT,
  SOIL_GRIDCEL,
  SOIL_ZWTVMOIST_ZWT,
  SOIL_ZWTVMOIST_MOIST,
  SOIL_SLOPE,
  SOIL_ASPECT,
  SOIL_EHORIZ,
  SOIL_WHORIZ,
  SOIL_MIN_DEPTH,
  SOIL_POROSITY_NODE,
  SOIL_EFFECTIVE_POROSITY_NODE,
  SOIL_WCR_FRACT,
  SOIL_WPWP_FRACT,
  SOIL_SUBSIDENCE,
  SOIL_NEW_SNOW_ALB,
  SOIL_SNOW_ALB_ACCUM_A,
  SOIL_SNOW_ALB_ACCUM_B,
  SOIL_SNOW_ALB_THAW_A,
  SOIL_SNOW_ALB_THAW_B,
  SOIL_MIN_RAIN_TEMP,
  SOIL_MAX_SNOW_TEMP,
  SOIL_PADJ,
	SOIL_T_LAPSE,
  SOIL_PGRAD,
  SOIL_GLAC_SURF_THICK,
  SOIL_GLAC_SURF_WE,
  SOIL_GLAC_KMIN,
  SOIL_GLAC_DK,
  SOIL_GLAC_A,
  SOIL_GLAC_ALBEDO,
  SOIL_GLAC_ROUGH,

  HRU_VEG_CON_CV,
  HRU_VEG_CON_ROOT,
  HRU_VEG_CON_ZONE_DEPTH,
  HRU_VEG_CON_ZONE_FRACT,
  HRU_VEG_CON_VEGINDEX,
  HRU_VEG_CON_SIGMA_SLOPE,
  HRU_VEG_CON_LAG_ONE,
  HRU_VEG_CON_FETCH,
  HRU_VEG_CON_LAKE,

  HRU_ISGLACIER,
  HRU_ISBARESOIL,

  HRUCELL_ASAT,
  HRUCELL_BASEFLOW,
  HRUCELL_INFLOW,
  HRUCELL_POT_EVAP,
  HRUCELL_RUNOFF,
  HRUCELL_ROOTMOIST,
  HRUCELL_WETNESS,
  HRUCELL_ZWT,
  HRUCELL_ZWT2,
  HRUCELL_ZWT3,

  AERO_SURFACE,
  AERO_OVERSTORY,

  LAYER_CS,
  LAYER_T,
  LAYER_EVAP,
  LAYER_KAPPA,
  LAYER_PHI,
  LAYER_ZWT,

  HRU_VEG_VAR_CANOPYEVAP,
  HRU_VEG_VAR_THROUGHFALL,

  SNOW_ALBEDO,
  SNOW_CANOPY_ALBEDO,
  SNOW_DEPTH,
  SNOW_MAX_SWQ,
  SNOW_SNOW,
  SNOW_STORE_COVERAGE,
  SNOW_STORE_SNOW,
  SNOW_STORE_SWQ,
  SNOW_SURF_TEMP_FBCOUNT,
  SNOW_SURF_TEMP_FBFLAG,
  SNOW_SWQ_SLOPE,
  SNOW_TMP_INT_STORAGE,
  SNOW_BLOWING_FLUX,
  SNOW_CANOPY_VAPOR_FLUX,
  SNOW_MASS_ERROR,
  SNOW_MELT,
  SNOW_QNET,
  SNOW_SURFACE_FLUX,
  SNOW_TRANSPORT,
  SNOW_VAPOR_FLUX,

  GLAC_COLD_CONTENT,
  GLAC_SURF_TEMP,
  GLAC_SURF_TEMP_FBCOUNT,
  GLAC_SURF_TEMP_FBFLAG,
  GLAC_QNET,
  GLAC_MASS_BALANCE,
  GLAC_ICE_MASS_BALANCE,
  GLAC_CUM_MASS_BALANCE,
  GLAC_ACCUMULATION,
  GLAC_MELT,
  GLAC_VAPOR_FLUX,
  GLAC_WATER_STORAGE,
  GLAC_OUTFLOW,
  GLAC_OUTFLOW_COEF,
  GLAC_INFLOW,

  ENERGY_ALBEDOLAKE,
  ENERGY_ALBEDOOVER,
  ENERGY_ALBEDOUNDER,
  ENERGY_CS,
  ENERGY_CS_NODE,
  ENERGY_FDEPTH,
  ENERGY_FROZEN,
  ENERGY_ICE_CONTENT,
  ENERGY_KAPPA,
  ENERGY_KAPPA_NODE,
  ENERGY_MOIST,
  ENERGY_NFROST,
  ENERGY_NTHAW,
  ENERGY_T_FBFLAG,
  ENERGY_T_FBCOUNT,
  ENERGY_T1_INDEX,
  ENERGY_TCANOPY,
  ENERGY_TCANOPY_FBFLAG,
  ENERGY_TCANOPY_FBCOUNT,
  ENERGY_TDEPTH,
  ENERGY_TFOLIAGE,
  ENERGY_TFOLIAGE_FBFLAG,
  ENERGY_TFOLIAGE_FBCOUNT,
  ENERGY_TSURF,
  ENERGY_TSURF_FBFLAG,
  ENERGY_TSURF_FBCOUNT,
  ENERGY_UNFROZEN,
  ENERGY_ADVECTED_SENSIBLE,
  ENERGY_ADVECTION,
  ENERGY_ATMOSERROR,
  ENERGY_ATMOSLATENT,
  ENERGY_ATMOSLATENTSUB,
  ENERGY_ATMOSSENSIBLE,
  ENERGY_CANOPY_ADVECTION,
  ENERGY_CANOPY_LATENT,
  ENERGY_CANOPY_LATENT_SUB,
  ENERGY_CANOPY_REFREEZE,
  ENERGY_CANOPY_SENSIBLE,
  ENERGY_DELTACC,
  ENERGY_DELTAH,
  ENERGY_ERROR,
  ENERGY_FUSION,
  ENERGY_GRND_FLUX,
  ENERGY_LATENT,
  ENERGY_LATENT_SUB,
  ENERGY_LONGWAVE,
  ENERGY_LONGOVERIN,
  ENERGY_LONGUNDERIN,
  ENERGY_LONGUNDEROUT,
  ENERGY_MELT_ENERGY,
  ENERGY_NETLONGATMOS,
  ENERGY_NETLONGOVER,
  ENERGY_NETLONGUNDER,
  ENERGY_NETSHORTATMOS,
  ENERGY_NETSHORTGRND,
  ENERGY_NETSHORTOVER,
  ENERGY_NETSHORTUNDER,
  ENERGY_OUT_LONG_CANOPY,
  ENERGY_OUT_LONG_SURFACE,
  ENERGY_REFREEZE_ENERGY,
  ENERGY_SENSIBLE,
  ENERGY_SHORTWAVE,
  ENERGY_SHORTOVERIN,
  ENERGY_SHORTUNDERIN,
  ENERGY_SNOW_FLUX,
  ENERGY_GLACIER_FLUX,
  ENERGY_DELTACC_GLAC,

	GLAC_MASS_BALANCE_EQN_TERMS,
  NUM_GLAC_MASS_BALANCE_INFO_TERMS,

  SOIL_PADJ_R,
  SOIL_PADJ_S,
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
	GLAC_MASS_BALANCE_EQN_DIM,
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
  int process(float* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(bool* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(const int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(const double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(const float* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(const bool* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);
  int process(const char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id);

  virtual void initializeOutput() = 0;
  virtual int write(const int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int write(const double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int write(const float* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int write(const char* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int write(const bool* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int processNewline();
  virtual StateHeader readHeader() = 0;
  virtual void notifyDimensionUpdate(StateVariables::StateVariableDimensionId dimension, int value = -1) {}
  virtual void initializeDimensionIndices() {}
  virtual int getCurrentDimensionIndex(StateVariables::StateVariableDimensionId dimension) { return -1; }
  virtual int seekToCell(int cellid, int* nVeg, int* nBand) = 0;
  virtual int read(int* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int read(double* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int read(float* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
  virtual int read(bool* data, int numValues, const StateVariables::StateMetaDataVariableIndices id) = 0;
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
