#include "WriteOutputNetCDF.h"
#include "user_def.h"

#if NETCDF_OUTPUT_AVAILABLE

#ifdef __unix__
//the uname function is unix specific
#include <sys/utsname.h>
#endif

#include <netcdf>
#include <ctime>
#include <map>
#include <sstream>

extern const char* version; // Defined in global.h
using namespace netCDF;

WriteOutputNetCDF::WriteOutputNetCDF(const ProgramState* state) : WriteOutputFormat(state), netCDF(NULL) {
  netCDFOutputFileName = state->options.NETCDF_FULL_FILE_PATH;
  mapping = getMapping(state->global_param.out_dt < 24);
  // The divisor will convert the difference to either hours or days respectively.
  timeIndexDivisor = state->global_param.out_dt < 24 ? (60 * 60) : (60 * 60 * 24);
}

WriteOutputNetCDF::~WriteOutputNetCDF() {
  if (netCDF != NULL) {
    delete netCDF;
  }
}

const char* WriteOutputNetCDF::getDescriptionOfOutputType() {
  return "NETCDF";
}

// Called once for each grid cell. This could potentially support parallel writes if a library supports it.
void WriteOutputNetCDF::openFile() {
  // Do not specify the type of file that will be read here (NcFile::nc4). The library will throw an exception
  // if it is provided for open types read or write (since the file was already created with a certain format).
  netCDF = new NcFile(netCDFOutputFileName.c_str(), NcFile::write);
}

//TODO: fill in all variables
std::map<std::string, VariableMetaData> WriteOutputNetCDF::getMapping(bool isHourly) {
  std::map<std::string, VariableMetaData> mapping;
  mapping["OUT_PREC"] =       VariableMetaData("kg m-2 s-1", "pr", "precipitation_flux", "Precipitation", "time: mean");
  mapping["OUT_EVAP"] =       VariableMetaData("kg m-2 s-1", "evspsbl", "water_evaporation_flux", "Evaporation", "time: mean");
  mapping["OUT_RUNOFF"] =     VariableMetaData("kg m-2 s-1", "mrro", "runoff_flux", "Total Runoff", "time: mean area: mean where land ");
  mapping["OUT_BASEFLOW"] =     VariableMetaData("mm", "FIXME_300581926477890014472171338486967809453", "", "baseflow out of the bottom layer", "");
  mapping["OUT_WDEW"] =     VariableMetaData("mm", "FIXME_307214154627712733439929426987234573336", "", "total moisture interception storage in canopy", "");
  mapping["OUT_SOIL_LIQ"] =   VariableMetaData("kg m-2", "mrlsl", "moisture_content_of_soil_layer", "Water Content of Soil Layer", "time: mean area: mean where land");
  mapping["OUT_RAD_TEMP"] =   VariableMetaData("K s-1", "tntr", "tendency_of_air_temperature_due_to_radiative_heating", "Tendency of Air Temperature due to Radiative Heating", "time: point");
  mapping["OUT_NET_SHORT"] =  VariableMetaData("W m-2", "rsntds", "net_downward_shortwave_flux_at_sea_water_surface", "Net Downward Shortwave Radiation at Sea Water Surface", "time: mean area: mean where sea");
  mapping["OUT_R_NET"] =      VariableMetaData("W m-2", "rtmt", "net_downward_radiative_flux_at_top_of_atmosphere_model", "Net Downward Flux at Top of Model", "time: mean");
  mapping["OUT_LATENT"] =     VariableMetaData("W m-2", "hfls", "surface_upward_latent_heat_flux", "Surface Upward Latent Heat Flux", "time: mean");
  mapping["OUT_EVAP_CANOP"] = VariableMetaData("kg m-2 s-1", "evspsblveg", "water_evaporation_flux_from_canopy", "Evaporation from Canopy", "time: mean area: mean where land");
  mapping["OUT_TRANSP_VEG"] = VariableMetaData("kg m-2 s-1", "tran", "transpiration_flux", "Transpiration", "time: mean area: mean where land");
  mapping["OUT_EVAP_BARE"] =  VariableMetaData("kg m-2 s-1", "evspsblsoi", "water_evaporation_flux_from_soil", "Water Evaporation from Soil", "time: mean area: mean where land");
  mapping["OUT_SUB_CANOP"] =    VariableMetaData("mm", "FIXME_317309329056271095262965656219868887316", "", "net sublimation from snow stored in canopy", "");
  mapping["OUT_SUB_SNOW"] =   VariableMetaData("kg m-2 s-1", "sbl", "surface_snow_and_ice_sublimation_flux", "Surface Snow and Ice Sublimation Flux", "time: mean");
  mapping["OUT_SENSIBLE"] =   VariableMetaData("W m-2", "hfss", "surface_upward_sensible_heat_flux", "Surface Upward Sensible Heat Flux", "time: mean");
  mapping["OUT_GRND_FLUX"] =  VariableMetaData("W m-2", "hfls", "surface_downward_latent_heat_flux", "Surface Downward Latent Heat Flux", "time: mean area: mean where ice_free_sea over sea");
  mapping["OUT_DELTAH"] =     VariableMetaData("W/m2", "FIXME_190446702042803383375556966506401492950", "", "rate of change in heat storage", "");
  mapping["OUT_FUSION"] =     VariableMetaData("W/m2", "FIXME_128664751816636936364937099409763952976", "", "net energy used to melt/freeze soil moisture", "");
  mapping["OUT_AERO_RESIST"] =    VariableMetaData("s/m", "FIXME_181616472461972626395272849195556773757", "", "\"scene\"canopy aerodynamic resistance", "");
  mapping["OUT_SURF_TEMP"] =  VariableMetaData("K", "ts", "surface_temperature", "Surface Temperature", "time: mean");
  mapping["OUT_ALBEDO"] =     VariableMetaData("%", "FIXME_332788171489585915540465750057800968512", "", "average surface albedo", "");
  mapping["OUT_REL_HUMID"] =  VariableMetaData("%", "hur", "relative_humidity", "Relative Humidity", "time: mean");
  mapping["OUT_IN_LONG"] =    VariableMetaData("W m-2", "rlds", "surface_downwelling_longwave_flux_in_air", "Surface Downwelling Longwave Radiation", "time: mean");
  mapping["OUT_AIR_TEMP"] =   VariableMetaData("K", "ta", "air_temperature", "Air Temperature", "time: mean");
  mapping["OUT_WIND"] =       VariableMetaData("m s-1", "sfcWind", "wind_speed", "Daily-Mean Near-Surface Wind Speed", "time: mean");
  mapping["OUT_SWE"] =        VariableMetaData("kg m-2 s-1", "prsn", "snowfall_flux", "Snowfall Flux", "time: mean");
  mapping["OUT_SNOW_DEPTH"] = VariableMetaData("m", "snd", "surface_snow_thickness", "Snow Depth", "time: mean area: mean where land");
  mapping["OUT_SNOW_CANOPY"] = VariableMetaData("kg m-2", "snveg", "veg_snow_amount", "Canopy Snow Amount", "time: mean area: mean where land"); // NOTE: PCIC made this up.
  mapping["OUT_SNOW_COVER"] = VariableMetaData("%", "snc", "surface_snow_area_fraction", "Snow Area Fraction", "time: mean");
  mapping["OUT_ADVECTION"] =    VariableMetaData("W/m2", "FIXME_220927444733091406228511039853387618452", "", "advected energy", "");
  mapping["OUT_DELTACC"] =    VariableMetaData("W/m2", "FIXME_182407129209214043804727417808316869671", "", "rate of change in cold content in snow pack", "");
  mapping["OUT_SNOW_FLUX"] =  VariableMetaData("W m-2", "hfdsn", "surface_downward_heat_flux_in_snow", "Downward Heat Flux into Snow Where Land over Land", "time: mean area: mean where land");
  mapping["OUT_RFRZ_ENERGY"] =    VariableMetaData("W/m2", "FIXME_77193432786456792120135396962917087055", "", "net energy used to refreeze liquid water in snowpack", "");
  mapping["OUT_MELT_ENERGY"] =    VariableMetaData("W/m2", "FIXME_110196731027167292658903755456932365807", "", "energy of fusion (melting) in snowpack", "");
  mapping["OUT_ADV_SENS"] =     VariableMetaData("W/m2", "FIXME_239036528020032962676483737753630104419", "", "net sensible flux advected to snow pack", "");
  mapping["OUT_LATENT_SUB"] = VariableMetaData("W m-2", "hfls", "surface_upward_latent_heat_flux", "Surface Upward Latent Heat Flux", "time: mean");
  mapping["OUT_SNOW_SURF_TEMP"] = VariableMetaData("K", "tsn", "temperature_in_surface_snow", "Snow Internal Temperature", "time: mean (with samples weighted by snow mass) area: mean where land");
  mapping["OUT_SNOW_PACK_TEMP"] =     VariableMetaData("C", "FIXME_301461005851069129713119502752676532492", "", "snow pack temperature", "");
  mapping["OUT_SNOW_MELT"] =  VariableMetaData("kg m-2 s-1", "snm", "surface_snow_melt_flux", "Surface Snow Melt", "time: mean area: mean where land");
  mapping["OUT_SUB_BLOWING"] = VariableMetaData("kg m-2 s-1", "sblow", "surface_snow_and_ice_sublimation_flux", "Blowing Snow Sublimation Flux", "time: mean area: mean where land"); // NOTE: PCIC made this up
  mapping["OUT_SUB_SURFACE"] = VariableMetaData("kg m-2 s-1", "sbl", "surface_snow_and_ice_sublimation_flux", "Surface Snow and Ice Sublimation Flux", "time: mean area: mean where land");
  mapping["OUT_SUB_SNOW"] =     VariableMetaData("mm", "FIXME_334274396941012849543418516660137802620", "", "total net sublimation from snow pack (surface and blowing)", "");
  mapping["OUT_FDEPTH"] =     VariableMetaData("cm", "FIXME_15212692387945942409305218575976782174", "", "depth of freezing fronts", "");
  mapping["OUT_TDEPTH"] =     VariableMetaData("cm", "FIXME_272423697841748397841695330645258350907", "", "depth of thawing fronts", "");
  mapping["OUT_SOIL_MOIST"] = VariableMetaData("kg m-2", "mrso", "soil_moisture_content", "Total Soil Moisture Content", "time: mean area: mean where land");
  mapping["OUT_SURF_FROST_FRAC"] =  VariableMetaData("%", "snc", "surface_snow_area_fraction", "Snow Area Fraction", "time: mean");
  mapping["OUT_SWE_BAND"] =         VariableMetaData("kg m-2", "lwsnl", "liquid_water_content_of_snow_layer", "Liquid Water Content of Snow Layer", "time: mean area: mean where land");
  mapping["OUT_SNOW_DEPTH_BAND"] =  VariableMetaData("m", "snd", "surface_snow_thickness", "Snow Depth", "time: mean area: mean where land", 0.01, 0, true);
  mapping["OUT_SNOW_CANOPY_BAND"] =     VariableMetaData("mm", "FIXME_8155395554501831256840682298328810829", "", "snow interception storage in canopy", "");
  mapping["OUT_ADVECTION_BAND"] =     VariableMetaData("W/m2", "FIXME_51785827765996758918063495782475679406", "", "advected energy", "");
  mapping["OUT_DELTACC_BAND"] =     VariableMetaData("W/m2", "FIXME_304856858194992817246188136208663173086", "", "change in cold content in snow pack", "");
  mapping["OUT_SNOW_FLUX_BAND"] =     VariableMetaData("W/m2", "FIXME_283161918025814344471631808155055213687", "", "energy flux through snow pack", "");
  mapping["OUT_RFRZ_ENERGY_BAND"] =     VariableMetaData("W/m2", "FIXME_126566890217291303977428213549908298041", "", "net energy used to refreeze liquid water in snowpack", "");
  mapping["OUT_NET_SHORT_BAND"] =     VariableMetaData("W/m2", "FIXME_134416319948353040325524030031863986330", "", "net downward shortwave flux", "");
  mapping["OUT_NET_LONG_BAND"] =    VariableMetaData("W/m2", "FIXME_258669125626847102505443325968262198066", "", "net downward longwave flux", "");
  mapping["OUT_ALBEDO_BAND"] =    VariableMetaData("%", "FIXME_126857714527222738073588574864979904389", "", "average surface albedo", "");
  mapping["OUT_LATENT_BAND"] =    VariableMetaData("W/m2", "FIXME_213329483764123762233599335367808460433", "", "net upward latent heat flux", "");
  mapping["OUT_SENSIBLE_BAND"] =    VariableMetaData("W/m2", "FIXME_236203365875103882750289090815248082673", "", "net upward sensible heat flux", "");
  mapping["OUT_GRND_FLUX_BAND"] =     VariableMetaData("W/m2", "FIXME_109892049149533581597022191892815756151", "", "net heat flux into ground", "");
  mapping["OUT_LAKE_ICE_TEMP"] =    VariableMetaData("C", "FIXME_34452480666930950419118115353827200627", "", "temperature of lake ice", "");
  mapping["OUT_LAKE_ICE_HEIGHT"] =    VariableMetaData("cm", "FIXME_93178255230432802544224048584307332491", "", "thickness of lake ice", "");
  mapping["OUT_LAKE_ICE_FRACT"] =     VariableMetaData("%", "FIXME_111826502398426868446569769827812489266", "", "fractional coverage of lake ice", "");
  mapping["OUT_LAKE_DEPTH"] =     VariableMetaData("m", "FIXME_84003765381102902266132885055520771964", "", "lake depth (distance between surface and deepest point)", "");
  mapping["OUT_LAKE_SURF_AREA"] =     VariableMetaData("m2", "FIXME_142255252263485917276779849136519162200", "", "lake surface area", "");
  mapping["OUT_LAKE_VOLUME"] =    VariableMetaData("m3", "FIXME_129734662285997138437784820469231183473", "", "lake volume", "");
  mapping["OUT_LAKE_SURF_TEMP"] =     VariableMetaData("C", "FIXME_76614810514209596805125427349782407719", "", "lake surface temperature", "");
  mapping["OUT_LAKE_EVAP"] =    VariableMetaData("mm", "FIXME_210668757305485901088872458979732011810", "", "net evaporation from lake surface", "");
  return mapping;
}

// The source version is set by the makefile dynamically based on the hg version control values.
std::string getSourceVersion() {
#ifdef SOURCE_VERSION
  return std::string(SOURCE_VERSION);
#else
  return std::string("No source version specified");
#endif
}

// The COMPILE_TIME macro is set by the makefile dynamically.
std::string getDateOfCompilation() {
#ifdef COMPILE_TIME
  return std::string(COMPILE_TIME);
#else
  return std::string("No compile time set");
#endif
}

// The MACHINE_INFO macro is set by the makefile dynamically.
std::string getCompilationMachineInfo() {
#ifdef MACHINE_INFO
  return std::string(MACHINE_INFO);
#else
  return std::string("Unspecified");
#endif
}

// Returns a descriptive string about the machine the model is currently being run on.
// This implementation is currently unix dependant because of the "uname" command.
// The code will still compile on another OS but will just return a non-descriptive string.
std::string getRuntimeMachineInfo() {
#ifdef __unix__
  struct utsname info;
  if (uname(&info) == 0) {
    return std::string(info.sysname) + " node: " + info.nodename + " release: " + info.release + " version: " + info.version + " machine: " + info.machine;
  } else {
    return std::string("Unspecified unix system");
  }
#else
  return std::string("Unspecified non-unix system");
#endif
}

void outputGlobalAttributeError(std::string error, const char** requiredAttributes) {
  std::string message = "Error: netCDF global attribute not specified: \"" + error + "\"\n";
  message += "When the output is in NETCDF format, certain global attributes must be specified in the global input file\n";
  message += "Make sure that at least ";
  int index = 0;
  while (requiredAttributes[index] != NULL) {
    message += std::string("\"") + requiredAttributes[index] + "\" ";
    index++;
  }
  message += "are explicitly specified.\n";
  message += "For example, put the following line in the global input file:\n";
  message += "NETCDF_ATTRIBUTE\t" + error + "\tSome string that may contain spaces";
  throw VICException(message);
}

void verifyGlobalAttributes(const NcFile& file) {
  const char* requiredAttributes [] = {"institution", "contact", "references", NULL}; // Null terminated list of attributes that must be specified.
  int index = 0;
  while (requiredAttributes[index] != NULL) {
    if (file.getAtt(std::string(requiredAttributes[index])).isNull()) {
      outputGlobalAttributeError(std::string(requiredAttributes[index]), requiredAttributes);
    }
    index++;
  }
}

void addGlobalAttributes(NcFile* netCDF, const ProgramState* state) {
  // Add global attributes here. (These could potentially be overwritten by inputs from the global file)
  netCDF->putAtt("title", "VIC output.");
  std::string source = std::string("VIC ") + version + ". ";
  source += "(Built from source: " + getSourceVersion() + ").";
  std::string history = "Created by " + source + ".\n";
  history += " Compiled on: " + getDateOfCompilation() + ".\n";
  history += " Compiled by machine: " + getCompilationMachineInfo() + ".\n";
  history += " Model run by machine: " + getRuntimeMachineInfo();

  netCDF->putAtt("model_start_year", netCDF::ncInt, state->global_param.startyear);
  netCDF->putAtt("model_start_month", netCDF::ncInt, state->global_param.startmonth);
  netCDF->putAtt("model_start_day", netCDF::ncInt, state->global_param.startday);
  netCDF->putAtt("model_start_hour", netCDF::ncInt, state->global_param.starthour);

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
}

int WriteOutputNetCDF::getLengthOfTimeDimension(const ProgramState* state) {
  dmy_struct endTime;
  endTime.year = state->global_param.endyear;
  endTime.month = state->global_param.endmonth;
  endTime.day = state->global_param.endday;
  endTime.hour = 0;
  endTime.day_in_year = 0;
  return getTimeIndex(&endTime, timeIndexDivisor, state) + 1;
}

// Called automatically only once at the beginning of the program.
void WriteOutputNetCDF::initializeFile(const ProgramState* state) {
  NcFile ncFile(netCDFOutputFileName.c_str(), NcFile::replace, NcFile::nc4);

  addGlobalAttributes(&ncFile, state);

  ncFile.putAtt("model_end_year", netCDF::ncInt, state->global_param.endyear);
  ncFile.putAtt("model_end_month", netCDF::ncInt, state->global_param.endmonth);
  ncFile.putAtt("model_end_day", netCDF::ncInt, state->global_param.endday);

  verifyGlobalAttributes(ncFile);

  // Set up the dimensions and variables.
  int timeSize = getLengthOfTimeDimension(state);
  int valuesSize = MAX_BANDS;
  fprintf(stderr, "Setting up grid dimensions, lat size: %ld, lon size: %ld, time: %d\n", (size_t)state->global_param.gridNumLatDivisions, (size_t)state->global_param.gridNumLonDivisions, timeSize);

  NcDim latDim = ncFile.addDim("lat", (size_t)state->global_param.gridNumLatDivisions);
  NcDim lonDim = ncFile.addDim("lon", (size_t)state->global_param.gridNumLonDivisions);
  NcDim bounds = ncFile.addDim("bnds", 2);
  NcDim timeDim = ncFile.addDim("time", timeSize);
  NcDim valuesDim = ncFile.addDim("depth", valuesSize);  // This dimension allows for variables which are actually arrays of values.


  // Define the coordinate variables.
  NcVar latVar = ncFile.addVar("lat", ncFloat, latDim);
  NcVar lonVar = ncFile.addVar("lon", ncFloat, lonDim);
  NcVar timeVar = ncFile.addVar("time", ncFloat, timeDim);
  NcVar valuesVar = ncFile.addVar("depth", ncFloat, valuesDim);

  //FIXME: add the bounds variables for each of time, lat, lon
  latVar.putAtt("axis", "Y");
  latVar.putAtt("units", "degrees_north");
  latVar.putAtt("standard name", "latitude");
  latVar.putAtt("long name", "latitude");
  latVar.putAtt("bounds", "lat_bnds");
  for (int i = 0; i < state->global_param.gridNumLatDivisions; i++) {
    std::vector<size_t> start, count;
    start.push_back(i);
    count.push_back(1);
    float value = state->global_param.gridStartLat + (i * state->global_param.gridStepLat);
    latVar.putVar(start, count, &value);
  }

  lonVar.putAtt("axis", "X");
  lonVar.putAtt("units", "degrees_east");
  lonVar.putAtt("standard name", "longitude");
  lonVar.putAtt("long name", "longitude");
  lonVar.putAtt("bounds", "lon_bnds");
  for (int i = 0; i < state->global_param.gridNumLonDivisions; i++) {
    std::vector<size_t> start, count;
    start.push_back(i);
    count.push_back(1);
    float value = state->global_param.gridStartLon + (i * state->global_param.gridStepLon);
    lonVar.putVar(start, count, &value);
  }

  std::stringstream ss;
  if (state->global_param.out_dt < 24) {
    ss << "hours since " << state->global_param.starthour << ":00 " << "on ";
  } else {
    ss << "days since ";
  }
  ss << state->global_param.startday << "/" << state->global_param.startmonth << "/" << state->global_param.startyear;

  timeVar.putAtt("axis", "T");
  timeVar.putAtt("standard name", "time");
  timeVar.putAtt("long name", "time");
  timeVar.putAtt("units", ss.str());
  timeVar.putAtt("bounds", "time_bnds");
  timeVar.putAtt("calendar", "gregorian");
  std::vector<size_t> start, count;
  start.push_back(0);
  count.push_back(1);
  for (int i = 0; i < timeSize; i++) {
    start[0] = i;
    float index = i;
    timeVar.putVar(start, count, &index);
  }

  valuesVar.putAtt("standard name", "z_dim");
  valuesVar.putAtt("long name", "array values");
  valuesVar.putAtt("units", "z_dim");
  for (int i = 0; i < valuesSize; i++) {
    start[0] = i;
    float index = i;
    valuesVar.putVar(start, count, &index);
  }


  out_data_struct* out_data_defaults = create_output_list(state);
  out_data_file_struct* out_file = set_output_defaults(out_data_defaults, state);

  // Define dimension orders.
  // If you change the ordering, make sure you also change the order that variables are written in the WriteOutputNetCDF::write_data() method.
  const NcDim dim3Vals [] = { latDim, lonDim, timeDim };
  const NcDim dim4Vals [] = { valuesDim, latDim, lonDim, timeDim };
  std::vector<NcDim> dimensions3(dim3Vals, dim3Vals + 3);
  std::vector<NcDim> dimensions4(dim4Vals, dim4Vals + 4);

  // Define a netCDF variable. For example, fluxes, snow.
  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {
    fprintf(stderr, "parent variable: %s\n", dataFiles[file_idx]->prefix);
    for (int var_idx = 0; var_idx < dataFiles[file_idx]->nvars; var_idx++) {
      const std::string varName = out_data_defaults[dataFiles[file_idx]->varid[var_idx]].varname;
      bool use4Dimensions = out_data_defaults[dataFiles[file_idx]->varid[var_idx]].nelem > 1;
      fprintf(stderr, "WriteOutputNetCDF initialize: adding variable: %s\n", varName.c_str());
      if (mapping.find(varName) == mapping.end()) {
        throw VICException("Error: could not find variable in mapping: " + varName);
      }
      VariableMetaData metaData = mapping[varName];
      if (metaData.isBands) {
        throw VICException("Error: bands not supported yet on variable mapping " + varName);
      } else {
        try {
          NcVar data = ncFile.addVar(metaData.name.c_str(), ncFloat, (use4Dimensions ? dimensions4 : dimensions3) );
          data.putAtt("long_name", metaData.longName.c_str());
          data.putAtt("units", metaData.units.c_str());
          data.putAtt("standard_name", metaData.standardName.c_str());
          data.putAtt("cell_methods", metaData.cellMethods.c_str());
          data.putAtt("_FillValue", ncFloat, 1e20f);
          data.putAtt("internal_vic_name", varName); // This is basically just for reference.
          data.putAtt("category", dataFiles[file_idx]->prefix);
          if (state->options.COMPRESS) {
            data.setCompression(false, true, 1); // Some reasonable compression level - not too intensive.
          }
        } catch (const netCDF::exceptions::NcException& except) {
          fprintf(stderr, "Error adding variable: %s with name: %s Internal netCDF exception\n", varName.c_str(), metaData.name.c_str());
          throw;
        }
      }
    }
  }
  free_out_data(&out_data_defaults);

  // The file will be automatically close when the NcFile object goes
  // out of scope. This frees up any internal netCDF resources
  // associated with the file, and flushes any buffers.
}

// Returns either the number of hours since the start time (if dt < 24)
// Or returns the number of days since the start time (if dt >= 24)
// Returns INVALID_INT on error.
int WriteOutputNetCDF::getTimeIndex(const dmy_struct* curTime, const int timeIndexDivisor, const ProgramState* state) {
  struct std::tm curTm = { 0, 0, curTime->hour, curTime->day, curTime->month - 1, curTime->year - 1900, 0, 0, 0, 0, "UTC" };  // Current date.
  struct std::tm startTm = { 0, 0, state->global_param.starthour, state->global_param.startday, state->global_param.startmonth - 1, state->global_param.startyear - 1900, 0, 0, 0, 0, "UTC" }; // Simulation start date.
  std::time_t curTimeT = std::mktime(&curTm);
  std::time_t startTimeT = std::mktime(&startTm);
  if (curTimeT != (std::time_t) (-1) && startTimeT != (std::time_t) (-1)) {
    // The function std::difftime returns the difference (curTimeT-startTimeT) in seconds.
    return std::difftime(curTimeT, startTimeT) / (timeIndexDivisor);
  }
  return INVALID_INT;
}

// This is called per time record.
void WriteOutputNetCDF::write_data(out_data_struct* out_data, const dmy_struct* dmy, int dt, const ProgramState* state) {

  if (netCDF == NULL) {
    fprintf(stderr, "Warning: could not write to netCDF file. Record %04i\t%02i\t%02i\t%02i\t Lat: %f, Lon: %f. File: \"%s\".\n",
        dmy->year, dmy->month, dmy->day, dmy->hour, lat, lon, netCDFOutputFileName.c_str());
    return;
  }

  const size_t timeIndex = getTimeIndex(dmy, timeIndexDivisor, state);
  const size_t lonIndex = longitudeToIndex(this->lon, state);
  const size_t latIndex = latitudeToIndex(this->lat, state);

  if (IS_INVALID((int)timeIndex) || IS_INVALID((int)lonIndex) || IS_INVALID((int)latIndex)) {
    std::stringstream s;
    s << "Error: Invalid index. timeIndex=" << timeIndex << ", lonIndex=" << lonIndex << ", latIndex=" << latIndex;
    throw VICException(s.str());
  }

  // Defines the dimension order of how variables are written. Only variables which have more than one value have the extra values dimension.
  // If you change the dimension ordering here, make sure that it is also changed in the WriteOutputNetCDF::initializeFile() method.
  // If the z dimension position changes also change the count4 vector update inside the nested loop below.
  const size_t start3Vals [] = { latIndex, lonIndex, timeIndex };     // (y, x, t)
  const size_t count3Vals [] = { 1,1,1 };
  const size_t start4Vals [] = { 0, latIndex, lonIndex, timeIndex };  // (z, y, x, t)
  const size_t count4Vals [] = { 1,1,1,1 };
  std::vector<size_t> start3(start3Vals, start3Vals + 3), count3(count3Vals, count3Vals + 3);
  std::vector<size_t> start4(start4Vals, start4Vals + 4), count4(count4Vals, count4Vals + 4);

  std::multimap<std::string, netCDF::NcVar> allVars = netCDF->getVars();

  // Loop over output files
  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {
    // Loop over this output file's data variables
    for (int var_idx = 0; var_idx < dataFiles[file_idx]->nvars; var_idx++) {
      // Loop over this variable's elements
      bool use4Dimensions = out_data[dataFiles[file_idx]->varid[var_idx]].nelem > 1;
      count4.at(0) = out_data[dataFiles[file_idx]->varid[var_idx]].nelem; // Change the number of values to write to the z dimension (array of values).

      try {
        NcVar variable = allVars.find(mapping[out_data[dataFiles[file_idx]->varid[var_idx]].varname].name)->second;
        if (use4Dimensions) {
          variable.putVar(start4, count4, out_data[dataFiles[file_idx]->varid[var_idx]].aggdata);
        } else {
          variable.putVar(start3, count3, out_data[dataFiles[file_idx]->varid[var_idx]].aggdata);
        }
      } catch (std::exception& e) {
        fprintf(stderr, "Error writing variable: %s, at latIndex: %d, lonIndex: %d, timeIndex: %d\n", mapping[out_data[dataFiles[file_idx]->varid[var_idx]].varname].name.c_str(), (int)latIndex, (int)lonIndex, (int)timeIndex);
        throw;
      }
    }
  }
}

void WriteOutputNetCDF::compressFiles() {
  // NetCDF output is not compressed at this point, if the compression option is on, then
  // each NetCDF variable will be compressed internally (built-in compression is a feature of NetCDF).
  // This method is intentionally empty.
}

void WriteOutputNetCDF::write_header(out_data_struct* out_data, const dmy_struct* dmy,
    const ProgramState* state) {
  // This is not applicable for netCDF output format. Intentionally empty.
}

#endif // NETCDF_OUTPUT_AVAILABLE
