#include "WriteOutputNetCDF.h"

#include <netcdf>
#include <ctime>
#include <map>
#include <sstream>

extern const char* version; // Defined in global.h
using namespace netCDF;

WriteOutputNetCDF::WriteOutputNetCDF(const ProgramState* state) : WriteOutputFormat(state), netCDF(NULL) {
  netCDFOutputFileName = state->options.NETCDF_FULL_FILE_PATH;
  mapping = getMapping(state->global_param.out_dt < 24);
}

WriteOutputNetCDF::~WriteOutputNetCDF() {
  if (netCDF != NULL) {
    fprintf(stderr, "NetCDF output file closed from destructor\n");
    delete netCDF;
  }
}

const char* WriteOutputNetCDF::getDescriptionOfOutputType() {
  return "NETCDF";
}

// Called once for each grid cell. This could potentially support parallel writes if a library supports it.
void WriteOutputNetCDF::openFile() {
  fprintf(stderr, "WriteOutputNetCDF::openFile() opening file: %s \n", netCDFOutputFileName.c_str());
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

void WriteOutputNetCDF::outputGlobalAttributeError(std::string error) {
  std::string message = "Error: netCDF global attribute not specified: \"" + error + "\"\n";
  message += "When the output is in NETCDF format, certain global attributes must be specified in the global input file\n";
  message += "For example, put the following line in the global:\n";
  message += "NETCDF_ATTRIBUTE\t" + error + "\tSome string that may contain spaces";
  throw VICException(message);
}

void WriteOutputNetCDF::verifyGlobalAttributes(NcFile& file) {
  if (file.getAtt("institution").isNull()) {
    outputGlobalAttributeError("institution");
  } else if (file.getAtt("contact").isNull()) {
    outputGlobalAttributeError("contact");
  } else if (file.getAtt("references").isNull()) {
    outputGlobalAttributeError("references");
  }
}

// Called automatically only once at the beginning of the program.
void WriteOutputNetCDF::initializeFile(const ProgramState* state) {
  NcFile ncFile(netCDFOutputFileName.c_str(), NcFile::replace, NcFile::nc4);

  // Add global attributes here. (These could potentially be overwritten by inputs from the global file)
  ncFile.putAtt("title", "VIC output.");

  //Add global attributes specified in the global input file.
  for (std::vector<std::pair<std::string, std::string> >::const_iterator it = state->global_param.netCDFGlobalAttributes.begin();
      it != state->global_param.netCDFGlobalAttributes.end(); ++it) {
    ncFile.putAtt(it->first, it->second);
  }

  std::string source = "VIC " + version + ". ";
  source += "(Built from source: " + getSourceVersion() + ").";
  std::string history = "Created by " + source;
  history += " on " + getDateOfCompilation() + ".";

  ncFile.putAtt("source", source.c_str());
  ncFile.putAtt("history", history.c_str());//TODO: arguments
  ncFile.putAtt("frequency", state->global_param.out_dt < 24 ? "hour" : "day");
  ncFile.putAtt("Conventions", "CF-1.6");

  verifyGlobalAttributes(ncFile);

  // Set up the dimensions and variables.

  fprintf(stderr, "Setting up grid dimensions, lat size: %ld, lon size: %ld\n", (size_t)state->global_param.gridNumLatDivisions, (size_t)state->global_param.gridNumLonDivisions);

  NcDim latDim = ncFile.addDim("lat", (size_t)state->global_param.gridNumLatDivisions);
  NcDim lonDim = ncFile.addDim("lon", (size_t)state->global_param.gridNumLonDivisions);
  NcDim bounds = ncFile.addDim("bnds", 2);
  NcDim timeDim = ncFile.addDim("time");    // Time is unlimited.


  // Define the coordinate variables.
  NcVar latVar = ncFile.addVar("lat", ncFloat, latDim);
  NcVar lonVar = ncFile.addVar("lon", ncFloat, lonDim);
  NcVar timeVar = ncFile.addVar("time", ncFloat, timeDim);
  //FIXME: add the bounds variables for each of time, lat, lon
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


  out_data_struct* out_data_defaults = create_output_list(state);
  out_data_file_struct* out_file = set_output_defaults(out_data_defaults, state);
  std::vector<NcDim> dimensions;
  dimensions.push_back(latDim);
  dimensions.push_back(lonDim);
  dimensions.push_back(timeDim);
  // Define a netCDF variable. For example, fluxes, snow.
  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {
    fprintf(stderr, "parent variable: %s\n", dataFiles[file_idx]->prefix);
    for (int var_idx = 0; var_idx < dataFiles[file_idx]->nvars; var_idx++) {
      const std::string varName = out_data_defaults[dataFiles[file_idx]->varid[var_idx]].varname;
      fprintf(stderr, "WriteOutputNetCDF initialize: adding variable: %s\n", varName.c_str());
      if (mapping.find(varName) == mapping.end()) {
        throw VICException("Error: could not find variable in mapping: " + varName);
      }
      VariableMetaData metaData = mapping[varName];
      if (metaData.isBands) {
        throw VICException("Error: bands not supported yet on variable mapping " + varName);
      } else {
        try {
          NcVar data = ncFile.addVar(metaData.name.c_str(), ncFloat, dimensions);
          data.putAtt("long_name", metaData.longName.c_str());
          data.putAtt("units", metaData.units.c_str());
          data.putAtt("standard_name", metaData.standardName.c_str());
          data.putAtt("cell_methods", metaData.cellMethods.c_str());
          data.putAtt("_FillValue", ncFloat, 1e20f);
          data.putAtt("internal_vic_name", varName); // This is basically just for reference.
          data.putAtt("category", dataFiles[file_idx]->prefix);
          if (state->options.COMPRESS) {
            data.setCompression(true, true, 3); // Some reasonable compression level - not too intensive.
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

//TODO: thoroughly test this function.
// Returns either the number of hours since the start time (if dt < 24)
// Or returns the number of days since the start time (if dt >= 24)
// Returns INVALID_INT on error.
int getTimeIndex(const dmy_struct* curTime, int dt, const ProgramState* state) {
  struct std::tm curTm = { 0, 0, curTime->hour, curTime->day, curTime->month, curTime->year - 1900, 0, 0, 0, 0, "UTC" };  // Current date.
  struct std::tm startTm = { 0, 0, state->global_param.starthour, state->global_param.startday, state->global_param.startmonth, state->global_param.startyear - 1900, 0, 0, 0, 0, "UTC" }; // Simulation start date.
  std::time_t curTimeT = std::mktime(&curTm);
  std::time_t startTimeT = std::mktime(&startTm);
  if (curTimeT != (std::time_t) (-1) && startTimeT != (std::time_t) (-1)) {
    // The divisor will convert the difference to either hours or days respectively.
    double divisor = dt < 24 ? (60 * 60) : (60 * 60 * 24);
    // The function std::difftime returns the difference (curTimeT-startTimeT) in seconds.
    return std::difftime(curTimeT, startTimeT) / (divisor);
  }
  return INVALID_INT;
}

int latitudeToIndex(double lat, const ProgramState* state) {
  return std::abs(lat - state->global_param.gridStartLat) / state->global_param.gridStepLat;
}

int longitudeToIndex(double lon, const ProgramState* state) {
  return std::abs(lon - state->global_param.gridStartLon) / state->global_param.gridStepLon;
}

//TODO: exceptions on error (terminate)
// This is called per time record.
void WriteOutputNetCDF::write_data(out_data_struct* out_data, const dmy_struct* dmy, int dt, const ProgramState* state) {

  if (netCDF == NULL) {
    fprintf(stderr, "Warning: could not write to netCDF file. Record %04i\t%02i\t%02i\t%02i\t Lat: %f, Lon: %f. File: \"%s\".\n",
        dmy->year, dmy->month, dmy->day, dmy->hour, lat, lon, netCDFOutputFileName.c_str());
    return;
  }

  int timeIndex = getTimeIndex(dmy, dt, state);
  int lonIndex = longitudeToIndex(this->lon, state);
  int latIndex = latitudeToIndex(this->lat, state);

  if (IS_INVALID(timeIndex) || IS_INVALID(lonIndex) || IS_INVALID(latIndex)) {
    std::stringstream s;
    s << "Error: Invalid index. timeIndex=" << timeIndex << ", lonIndex=" << lonIndex << ", latIndex=" << latIndex;
    throw VICException(s.str());
  }

  std::vector<size_t> start, count;
  start.push_back(latIndex);
  start.push_back(lonIndex);
  start.push_back(timeIndex);
  start.push_back(0);
  count.push_back(1);
  count.push_back(1);
  count.push_back(1);
  count.push_back(0);

  std::multimap<std::string, netCDF::NcVar> allVars = netCDF->getVars();

  // Loop over output files
  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {
    // Loop over this output file's data variables
    for (int var_idx = 0; var_idx < dataFiles[file_idx]->nvars; var_idx++) {
      // Loop over this variable's elements
      //for (int elem_idx = 0; elem_idx < out_data[dataFiles[file_idx]->varid[var_idx]].nelem; elem_idx++) {
      count[3] = out_data[dataFiles[file_idx]->varid[var_idx]].nelem;
      try {
        NcVar variable = allVars.find(mapping[out_data[dataFiles[file_idx]->varid[var_idx]].varname].name)->second;//netCDF->getVar(mapping[out_data[dataFiles[file_idx]->varid[var_idx]].varname].name);
        variable.putVar(start, count, out_data[dataFiles[file_idx]->varid[var_idx]].aggdata);
      } catch (std::exception& e) {
        fprintf(stderr, "Error writing variable: %s, at latIndex: %d, lonIndex: %d, timeIndex: %d\n", mapping[out_data[dataFiles[file_idx]->varid[var_idx]].varname].name.c_str(), latIndex, lonIndex, timeIndex);
        throw;
      }
      //variable.putVar(out_data[dataFiles[file_idx]->varid[var_idx]].aggdata, latIndex, lonIndex, timeIndex);
      //}
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


