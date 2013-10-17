#include "WriteOutputNetCDF.h"
#include <netcdf>  //TODO: switch to netcdf-4 library

#include <ctime>
#include <map>
#include <sstream>

using namespace netCDF;

WriteOutputNetCDF::WriteOutputNetCDF(const ProgramState* state) : WriteOutputFormat(state), netCDF(NULL) {
  netCDFOutputFileName = state->options.NETCDF_FULL_FILE_PATH;
}

WriteOutputNetCDF::~WriteOutputNetCDF() {
  fprintf(stderr, "NetCDF output file closed from destructor\n");
  if (netCDF != NULL) {
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

using std::string;
class VariableMetaData {
public:
  VariableMetaData() {}
  VariableMetaData(string units, string name, string standardName, string longName, string cellMethods, bool isBands = false): units(units), name(name), standardName(standardName), longName(longName), cellMethods(cellMethods), isBands(isBands) {}
  string units, name, standardName, longName, cellMethods;
  bool isBands;
};

//TODO: fill in all variables
std::map<std::string, VariableMetaData> getMapping() {
  std::map<std::string, VariableMetaData> mapping;
  mapping["OUT_PREC"] = VariableMetaData("OUT_PREC", "pr", "precipitation_flux", "Precipitation", "time: mean");
  return mapping;
}

//TODO: these attributes should be configurable from the global file.
// TODO: throw exception on error
// Called automatically only once at the beginning of the program.
void WriteOutputNetCDF::initializeFile(const ProgramState* state) {
  NcFile ncFile(netCDFOutputFileName.c_str(), NcFile::replace, NcFile::nc4);

  // Add global attributes here.
  ncFile.putAtt("title", "VIC output.");
  std::string source = "VIC "; //TODO: grab the source version here
  //source += version;            //TODO: get actual version
  ncFile.putAtt("source", source.c_str());
  ncFile.putAtt("history", ("Created by " + source).c_str());//TODO: arguments, creation date and time
  ncFile.putAtt("references", "http://www.pacificclimate.org");  //TODO: point to an actual page for this source code
  ncFile.putAtt("institution", "Pacific Climate Impacts Consortium (PCIC), Victoria, BC, www.pacificclimate.org");
  ncFile.putAtt("contact", "jstone@uvic.ca");
  ncFile.putAtt("frequency", state->global_param.out_dt < 24 ? "hour" : "day");
  ncFile.putAtt("Conventions", "CF-1.6");
  // When we create netCDF dimensions, we get back a pointer to an
  // NcDim for each one.
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
  timeVar.putAtt("units", ss.str().c_str());
  timeVar.putAtt("bounds", "time_bnds");
  timeVar.putAtt("calendar", "gregorian");


  std::map<std::string, VariableMetaData> mapping = getMapping();
  out_data_struct* out_data_defaults = create_output_list(state);
  out_data_file_struct* out_file = set_output_defaults(out_data_defaults, state);
  std::vector<NcDim> dimensions;
  dimensions.push_back(latDim);
  dimensions.push_back(lonDim);
  dimensions.push_back(timeDim);
  // Define a netCDF variable. For example, fluxes, snow.
  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {
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
        NcVar data = ncFile.addVar(metaData.name.c_str(), ncFloat, dimensions);
        data.putAtt("long_name", metaData.longName.c_str());
        data.putAtt("units", metaData.units.c_str());
        data.putAtt("standard_name", metaData.standardName.c_str());
        data.putAtt("cell_methods", metaData.cellMethods.c_str());
        data.putAtt("_FillValue", ncFloat, 1e20f);
        if (state->options.COMPRESS) {
          data.setCompression(true, true, 3); // Some reasonable compression level.
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
  start.push_back(0);
  start.push_back(latIndex);
  start.push_back(lonIndex);
  start.push_back(timeIndex);
  count.push_back(0); // This changes in the loop below.
  count.push_back(1);
  count.push_back(1);
  count.push_back(1);

  // Loop over output files
  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {
    // Loop over this output file's data variables
    for (int var_idx = 0; var_idx < dataFiles[file_idx]->nvars; var_idx++) {
      // Loop over this variable's elements
      //for (int elem_idx = 0; elem_idx < out_data[dataFiles[file_idx]->varid[var_idx]].nelem; elem_idx++) {
        count[0] = out_data[dataFiles[file_idx]->varid[var_idx]].nelem;
        NcVar variable = netCDF->getVar(out_data[dataFiles[file_idx]->varid[var_idx]].varname);
        variable.putVar(start, count, out_data[dataFiles[file_idx]->varid[var_idx]].aggdata);
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


