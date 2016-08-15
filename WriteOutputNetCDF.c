#include "WriteOutputNetCDF.h"
#include "user_def.h"

#if NETCDF_OUTPUT_AVAILABLE

#ifdef __unix__
//the uname function is unix specific
#include <sys/utsname.h>
#endif

#include <netcdf>
#include <ctime>
//#include <map>
#include <sstream>

using namespace netCDF;

WriteOutputNetCDF::WriteOutputNetCDF(const ProgramState* state) : WriteOutputFormat(state), netCDF(NULL) {
  netCDFOutputFileName = state->options.NETCDF_FULL_FILE_PATH;
  // The divisor will convert the difference to sub-daily (e.g. hourly, 3/4/6/8/12-hourly) or daily, respectively.
  timeIndexDivisor = state->global_param.out_dt < 24 ? (60 * 60 * state->global_param.out_dt) : (60 * 60 * 24); //new (*state->global_param.dt)
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

// The source version is set by the makefile at compile time based on the hg version control values.
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
  if(state->options.OUTPUT_FORCE){
	  netCDF->putAtt("title", "VIC meteorological forcing disaggregator mode output.");
  }
  else {
	  netCDF->putAtt("title", "VIC model run output.");
  }
  std::string source = std::string("VIC ");
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
  std::string frequency = std::to_string(state->global_param.out_dt);
  netCDF->putAtt("frequency", state->global_param.out_dt < 24 ? frequency+=" hour" : "day");
  netCDF->putAtt("Conventions", "CF-1.6");
}

int WriteOutputNetCDF::getLengthOfTimeDimension(const ProgramState* state) {
  dmy_struct endTime;
  endTime.year = state->global_param.endyear;
  endTime.month = state->global_param.endmonth;
  endTime.day = state->global_param.endday;
//  endTime.hour = 0;
  endTime.hour = 23; //new
  endTime.day_in_year = 0;
  return getTimeIndex(&endTime, timeIndexDivisor, state) + 1;
}

// Called automatically only once at the beginning of the program.
void WriteOutputNetCDF::initializeFile(const ProgramState* state, const out_data_struct* out_data_defaults) {

	NcFile ncFile(netCDFOutputFileName.c_str(), NcFile::replace, NcFile::nc4);

  addGlobalAttributes(&ncFile, state);

  ncFile.putAtt("model_end_year", netCDF::ncInt, state->global_param.endyear);
  ncFile.putAtt("model_end_month", netCDF::ncInt, state->global_param.endmonth);
  ncFile.putAtt("model_end_day", netCDF::ncInt, state->global_param.endday);

  verifyGlobalAttributes(ncFile);

  // Set up the dimensions and variables.
  int timeSize = getLengthOfTimeDimension(state);
  //int timeSize = state->global_param.nrecs; //new
  int valuesSize = MAX_BANDS;
  fprintf(stderr, "Setting up grid dimensions, lat size: %ld, lon size: %ld, time: %d\n", (size_t)state->global_param.gridNumLatDivisions, (size_t)state->global_param.gridNumLonDivisions, timeSize);

  NcDim latDim = ncFile.addDim("lat", (size_t)state->global_param.gridNumLatDivisions);
  NcDim lonDim = ncFile.addDim("lon", (size_t)state->global_param.gridNumLonDivisions);
  NcDim bounds = ncFile.addDim("bnds", 2);
  NcDim timeDim = ncFile.addDim("time", timeSize);
  NcDim valuesDim = ncFile.addDim("depth", valuesSize);  // This dimension allows for variables which are actually arrays of values.


  // Define the coordinate variables.
  NcVar latVar = ncFile.addVar("lat", ncDouble, latDim);
  NcVar lonVar = ncFile.addVar("lon", ncDouble, lonDim);
  NcVar timeVar = ncFile.addVar("time", ncFloat, timeDim);
  NcVar valuesVar = ncFile.addVar("depth", ncFloat, valuesDim);

  latVar.putAtt("axis", "Y");
  latVar.putAtt("units", "degrees_north");
  latVar.putAtt("standard name", "latitude");
  latVar.putAtt("long name", "latitude");
  latVar.putAtt("bounds", "lat_bnds");

  for (int i = 0; i < state->global_param.gridNumLatDivisions; i++) {
    std::vector<size_t> start, count;
    start.push_back(i);
    count.push_back(1);
    double value = state->global_param.gridStartLat + (i * state->global_param.gridStepLat);
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
    double value = state->global_param.gridStartLon + (i * state->global_param.gridStepLon);
    lonVar.putVar(start, count, &value);
  }

  std::stringstream ss;
  if (state->global_param.out_dt < 24) {
    ss << "hours since ";
  } else {
    ss << "days since ";
  }
  ss << state->global_param.startyear << "-" << state->global_param.startmonth << "-" << state->global_param.startday;
  if (state->global_param.out_dt < 24)
    ss << " " << state->global_param.starthour << ":00";

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
    float index = state->global_param.out_dt < 24 ? (i * state->global_param.out_dt) : i;
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

  // Define dimension orders.
  // If you change the ordering, make sure you also change the order that variables are written in the WriteOutputNetCDF::write_data() method.
  const NcDim dim3Vals [] = { timeDim, latDim, lonDim };
  const NcDim dim4Vals [] = { timeDim, valuesDim, latDim, lonDim };
  std::vector<NcDim> dimensions3(dim3Vals, dim3Vals + 3);
  std::vector<NcDim> dimensions4(dim4Vals, dim4Vals + 4);

  // Define a netCDF variable. For example, fluxes, snow.
  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {
    for (int var_idx = 0; var_idx < dataFiles[file_idx]->nvars; var_idx++) {
      const std::string varName = out_data_defaults[dataFiles[file_idx]->varid[var_idx]].varname;
      bool use4Dimensions = out_data_defaults[dataFiles[file_idx]->varid[var_idx]].nelem > 1;

      if (state->output_mapping.find(varName) == state->output_mapping.end()) {
        throw VICException("Error: could not find variable in output_mapping: " + varName);
      }
      VariableMetaData metaData = state->output_mapping.at(varName);
#if VERBOSE
      fprintf(stderr, "WriteOutputNetCDF initializeFile: adding variable: %s (NetCDF output variable name: %s)\n", varName.c_str(), metaData.name.c_str());
#endif
      if (metaData.isBands) {
        throw VICException("Error: bands not supported yet on variable output mapping " + varName);
      } else {
        try {
          NcVar data = ncFile.addVar(metaData.name.c_str(), ncFloat, (use4Dimensions ? dimensions4 : dimensions3) );
          data.putAtt("long_name", metaData.longName.c_str());
          data.putAtt("units", metaData.units.c_str());
          data.putAtt("standard_name", metaData.standardName.c_str());
          data.putAtt("cell_methods", metaData.cellMethods.c_str());
          data.putAtt("_FillValue", ncFloat, NETCDF_FILL_VALUE);
          data.putAtt("internal_vic_name", varName); // This should be the same as that given next to OUTVAR in the global file (and the key used to find variable metadata in state->mapping)
          data.putAtt("category", dataFiles[file_idx]->prefix);
          if (state->options.COMPRESS) {
            data.setCompression(false, true, 5); // Some reasonable compression level - not too intensive.
          }
        } catch (const netCDF::exceptions::NcException& except) {
          fprintf(stderr, "Error adding variable: %s with name: %s Internal netCDF exception\n", varName.c_str(), metaData.name.c_str());
          throw;
        }
      }
    }
  }
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
void WriteOutputNetCDF::write_data_one_cell(out_data_struct* out_data, const dmy_struct* dmy, int numTimeRecords, const ProgramState* state) {

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
  const size_t start3Vals [] = { timeIndex, latIndex, lonIndex };     // (t, y, x)
  const size_t count3Vals [] = { numTimeRecords,1,1 };
  const size_t start4Vals [] = { timeIndex, 0, latIndex, lonIndex };  // (t, z, y, x)
  const size_t count4Vals [] = { numTimeRecords,1,1,1 };
  std::vector<size_t> start3(start3Vals, start3Vals + 3), count3(count3Vals, count3Vals + 3);
  std::vector<size_t> start4(start4Vals, start4Vals + 4), count4(count4Vals, count4Vals + 4);

  std::multimap<std::string, netCDF::NcVar> allVars = netCDF->getVars();

  // Loop over output files
  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {
    // Loop over this output file's data variables
    for (int var_idx = 0; var_idx < dataFiles[file_idx]->nvars; var_idx++) {
      // Loop over this variable's elements
      bool use4Dimensions = out_data[dataFiles[file_idx]->varid[var_idx]].nelem > 1;
      count4.at(1) = out_data[dataFiles[file_idx]->varid[var_idx]].nelem; // Change the number of values to write to the z dimension (array of values).

      try {
      	NcVar variable = allVars.find(state->output_mapping.at(out_data[dataFiles[file_idx]->varid[var_idx]].varname).name)->second;

        if (use4Dimensions) {
          variable.putVar(start4, count4, out_data[dataFiles[file_idx]->varid[var_idx]].aggdata);
        } else {
          variable.putVar(start3, count3, out_data[dataFiles[file_idx]->varid[var_idx]].aggdata);
        }
      } catch (std::exception& e) {
        fprintf(stderr, "Error writing variable: %s, at latIndex: %d, lonIndex: %d, timeIndex: %d\n", state->output_mapping.at(out_data[dataFiles[file_idx]->varid[var_idx]].varname).name.c_str(), (int)latIndex, (int)lonIndex, (int)timeIndex);
        throw;
      }
    }
  }
}

void WriteOutputNetCDF::write_data_all_cells(std::vector<out_data_struct*>& all_out_data, out_data_file_struct *out_data_files_template, const dmy_struct *dmy, int dt, const ProgramState* state) {

	if (netCDF == NULL) {
		fprintf(stderr, "Warning: could not write to netCDF file. Record %04i\t%02i\t%02i\t%02i\t. File: \"%s\".\n",
				dmy->year, dmy->month, dmy->day, dmy->hour, state->options.NETCDF_FULL_FILE_PATH);
		return;
	}

	const size_t timeIndex = getTimeIndex(dmy, timeIndexDivisor, state);

	const size_t start3Vals [] = { timeIndex, 0, 0 };     // (t, y, x)
	const size_t count3Vals [] = { 1,(size_t) state->global_param.gridNumLatDivisions,(size_t) state->global_param.gridNumLonDivisions };
	const size_t start4Vals [] = { timeIndex, 0, 0, 0 };  // (t, z, y, x)
	const size_t count4Vals [] = { 1,1,(size_t) state->global_param.gridNumLatDivisions,(size_t) state->global_param.gridNumLonDivisions };
	std::vector<size_t> start3(start3Vals, start3Vals + 3), count3(count3Vals, count3Vals + 3);
	std::vector<size_t> start4(start4Vals, start4Vals + 4), count4(count4Vals, count4Vals + 4);

	std::multimap<std::string, netCDF::NcVar> allVars = netCDF->getVars();
	int num_cells = state->global_param.gridNumLatDivisions * state->global_param.gridNumLonDivisions;
	bool *modeled_cell_mask_ptr;

	// Loop through (legacy) out_data_files_template for listing of output variables
	for (int file_idx = 0; file_idx < state->options.Noutfiles; file_idx++) {
		// Loop over this output file's data variables
		for (int var_idx = 0; var_idx < out_data_files_template[file_idx].nvars; var_idx++) {

//			  		fprintf(stderr, "WriteOutputAllCells::write_data: varname = %s,  timeIndex = %d\n",all_out_data[0][out_data_files_template[file_idx].varid[var_idx]].varname, timeIndex);

			// Create temporary array of data for this variable across all cells
			float *vardata, *vardata_ptr;
			int varnumelem = all_out_data[0][out_data_files_template[file_idx].varid[var_idx]].nelem;
			bool use4Dimensions = varnumelem > 1;
			int vardatasize = varnumelem * num_cells;
			int modeled_cell_idx = 0; // index among modeled cells

			vardata = new float [vardatasize];
			vardata_ptr = vardata;

			// Interleave data from each cell for this variable in temporary array vardata
			for (int elem=0; elem<varnumelem; elem++) {
				modeled_cell_mask_ptr = state->modeled_cell_mask;
				modeled_cell_idx = 0;
				for (int cell_idx = 0; cell_idx < num_cells; cell_idx++) {
					if (*modeled_cell_mask_ptr) { // Check if this cell is marked as modeled in modeled_cell_mask
						*vardata_ptr = all_out_data[modeled_cell_idx][out_data_files_template[file_idx].varid[var_idx]].aggdata[elem];
						modeled_cell_idx++;
					}
					else {
						*vardata_ptr = NETCDF_FILL_VALUE;
					}
					modeled_cell_mask_ptr++;
					vardata_ptr++;
				}
			}
			// Write data to file for this variable
			try {
				NcVar variable = allVars.find(state->output_mapping.at(all_out_data[0][out_data_files_template[file_idx].varid[var_idx]].varname).name)->second;

				if (use4Dimensions) {
					count4.at(1) = varnumelem;  // Set the number of values to write to the z dimension
					variable.putVar(start4, count4, vardata);
				} else {
					variable.putVar(start3, count3, vardata);
				}
			} catch (std::exception& e) {
				fprintf(stderr, "Error writing variable: %s, at timeIndex: %d\n", all_out_data[0][out_data_files_template[file_idx].varid[var_idx]].varname, (int)timeIndex);

				throw;
			}
			delete [] vardata;
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
