#include "WriteOutputAllCells.h"

#include <netcdf>
#include <ctime>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace netCDF;

//namespace netCDF {
//  class NcFile;
//}


//WriteOutputAllCells::WriteOutputAllCells(const ProgramState* state) : WriteOutputFormat(state), netCDF(NULL) {
WriteOutputAllCells::WriteOutputAllCells(const ProgramState* state) : netCDF(NULL) {

	netCDFOutputFileName = state->options.NETCDF_FULL_FILE_PATH;
  // The divisor will convert the difference to sub-daily (e.g. hourly, 3/4/6/8/12-hourly) or daily, respectively.
  timeIndexDivisor = state->global_param.out_dt < 24 ? (60 * 60 * state->global_param.out_dt) : (60 * 60 * 24); //new (*state->global_param.dt)
}

// Returns either the number of hours since the start time (if dt < 24)
// Or returns the number of days since the start time (if dt >= 24)
// Returns INVALID_INT on error.
int WriteOutputAllCells::getTimeIndex(const dmy_struct* curTime, const int timeIndexDivisor, const ProgramState* state){
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

void WriteOutputAllCells::openFile() {
  // Do not specify the type of file that will be read here (NcFile::nc4). The library will throw an exception
  // if it is provided for open types read or write (since the file was already created with a certain format).
  netCDF = new NcFile(netCDFOutputFileName.c_str(), NcFile::write);
}

//void WriteOutputAllCells::write_all_cells_output(out_data_struct *all_out_data, , const dmy_struct* dmy, int dt, const ProgramState* state) {
void WriteOutputAllCells::write_all_cells_output(out_data_struct *all_out_data, out_data_file_struct* out_data_file, const dmy_struct* dmy, int dt, const ProgramState* state) {

  if (netCDF == NULL) {
    fprintf(stderr, "Warning: could not write to netCDF file. Record %04i\t%02i\t%02i\t%02i\t. File: \"%s\".\n",
        dmy->year, dmy->month, dmy->day, dmy->hour, state->options.NETCDF_FULL_FILE_PATH);
    return;
  }

  const size_t timeIndex = getTimeIndex(dmy, timeIndexDivisor, state);
//  const size_t lonIndex = longitudeToIndex(this->lon, state);
//  const size_t latIndex = latitudeToIndex(this->lat, state);
//
//  if (IS_INVALID((int)timeIndex) || IS_INVALID((int)lonIndex) || IS_INVALID((int)latIndex)) {
//    std::stringstream s;
//    s << "Error: Invalid index. timeIndex=" << timeIndex << ", lonIndex=" << lonIndex << ", latIndex=" << latIndex;
//    throw VICException(s.str());
//  }

  // Defines the dimension order of how variables are written. Only variables which have more than one value have the extra values dimension.
  // If you change the dimension ordering here, make sure that it is also changed in the WriteOutputNetCDF::initializeFile() method.
  // If the z dimension position changes also change the count4 vector update inside the nested loop below.
//  const size_t start3Vals [] = { timeIndex, latIndex, lonIndex };     // (t, y, x)
//  const size_t count3Vals [] = { 1,1,1 };
//  const size_t start4Vals [] = { timeIndex, 0, latIndex, lonIndex };  // (t, z, y, x)
//  const size_t count4Vals [] = { 1,1,1,1 };
  const size_t start3Vals [] = { timeIndex, 0, 0 };     // (t, y, x)
  const size_t count3Vals [] = { 1,(size_t) state->global_param.gridNumLatDivisions,(size_t) state->global_param.gridNumLonDivisions };
  const size_t start4Vals [] = { timeIndex, 0, 0, 0 };  // (t, z, y, x)
  // NOTE: the second element of count4Vals will change according to the variable (see line 412);
  // how do we account for this when trying to write out all variables for all cells for one time step?
  const size_t count4Vals [] = { 1,1,(size_t) state->global_param.gridNumLonDivisions,(size_t) state->global_param.gridNumLonDivisions };
  std::vector<size_t> start3(start3Vals, start3Vals + 3), count3(count3Vals, count3Vals + 3);
  std::vector<size_t> start4(start4Vals, start4Vals + 4), count4(count4Vals, count4Vals + 4);

  std::multimap<std::string, netCDF::NcVar> allVars = netCDF->getVars();

  // Loop over output files NOTE: DOES THIS MULTI-CELL WRITE BREAK ASCII MODE?
//  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {
    // Loop over this output file's data variables
    for (int var_idx = 0; var_idx < out_data_file->nvars; var_idx++) {
      // Loop over this variable's elements
     // bool use4Dimensions = all_out_data[dataFiles[file_idx]->varid[var_idx]].nelem > 1;
     // count4.at(1) = all_out_data[dataFiles[file_idx]->varid[var_idx]].nelem; // Change the number of values to write to the z dimension (array of values).
    	int varnumelem = all_out_data[out_data_file->varid[var_idx]].nelem;
      bool use4Dimensions = varnumelem > 1;
    	// Assign the length of the data to be written for this variable
//      unsigned int num_cells = sizeof(all_out_data) / sizeof(*all_out_data);
      unsigned int num_cells = 2;

    	size_t vardatasize = varnumelem * num_cells;
    	count4.at(1) = vardatasize;
      // Create temporary array of data for this variable across all cells
      double *vardata = new double[vardatasize];
    	// Fill vardata
    	for (unsigned int cellidx = 0; cellidx < num_cells; cellidx++) {
    		vardata[cellidx*varnumelem] = *all_out_data[out_data_file->varid[var_idx]].aggdata;
    	}

    	// Write data to file for this variable
      try {
      	NcVar variable = allVars.find(state->output_mapping.at(all_out_data[out_data_file->varid[var_idx]].varname).name)->second;

        if (use4Dimensions) {
//          variable.putVar(start4, count4, all_out_data[dataFiles[file_idx]->varid[var_idx]].aggdata);
          variable.putVar(start4, count4, vardata);

        } else {
//          variable.putVar(start3, count3, all_out_data[dataFiles[file_idx]->varid[var_idx]].aggdata);
          variable.putVar(start3, count3, vardata);
        }
      } catch (std::exception& e) {
        fprintf(stderr, "Error writing variable: %s, at timeIndex: %d\n", state->output_mapping.at(all_out_data[out_data_file->varid[var_idx]].varname).name.c_str(), (int)timeIndex);
        throw;
      }
    }
  }
//}

// OLD VERSION:
//void write_all_cells_output(WriteOutputFormat *outputFormat, out_data_struct *all_out_data, const dmy_struct* dmy, int dt, const ProgramState* state) {
//
//  if (netCDF == NULL) {
//    fprintf(stderr, "Warning: could not write to netCDF file. Record %04i\t%02i\t%02i\t%02i\t. File: \"%s\".\n",
//        dmy->year, dmy->month, dmy->day, dmy->hour, state->options.NETCDF_FULL_FILE_PATH);
//    return;
//  }
//
//  const size_t timeIndex = getTimeIndex(dmy, state);
////  const size_t lonIndex = longitudeToIndex(this->lon, state);
////  const size_t latIndex = latitudeToIndex(this->lat, state);
////
////  if (IS_INVALID((int)timeIndex) || IS_INVALID((int)lonIndex) || IS_INVALID((int)latIndex)) {
////    std::stringstream s;
////    s << "Error: Invalid index. timeIndex=" << timeIndex << ", lonIndex=" << lonIndex << ", latIndex=" << latIndex;
////    throw VICException(s.str());
////  }
//
//  // Defines the dimension order of how variables are written. Only variables which have more than one value have the extra values dimension.
//  // If you change the dimension ordering here, make sure that it is also changed in the WriteOutputNetCDF::initializeFile() method.
//  // If the z dimension position changes also change the count4 vector update inside the nested loop below.
////  const size_t start3Vals [] = { timeIndex, latIndex, lonIndex };     // (t, y, x)
////  const size_t count3Vals [] = { 1,1,1 };
////  const size_t start4Vals [] = { timeIndex, 0, latIndex, lonIndex };  // (t, z, y, x)
////  const size_t count4Vals [] = { 1,1,1,1 };
//  const size_t start3Vals [] = { timeIndex, 0, 0 };     // (t, y, x)
//  const size_t count3Vals [] = { 1,state->global_param.gridNumLatDivisions,state->global_param.gridNumLonDivisions };
//  const size_t start4Vals [] = { timeIndex, 0, 0, 0 };  // (t, z, y, x)
//  // NOTE: the second element of count4Vals will change according to the variable (see line 412);
//  // how do we account for this when trying to write out all variables for all cells for one time step?
//  const size_t count4Vals [] = { 1,1,state->global_param.gridNumLonDivisions,state->global_param.gridNumLonDivisions };
//  std::vector<size_t> start3(start3Vals, start3Vals + 3), count3(count3Vals, count3Vals + 3);
//  std::vector<size_t> start4(start4Vals, start4Vals + 4), count4(count4Vals, count4Vals + 4);
//
//  std::multimap<std::string, netCDF::NcVar> allVars = netCDF->getVars();
//
//  // Loop over output files NOTE: DOES THIS MULTI-CELL WRITE BREAK ASCII MODE?
//  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {
//    // Loop over this output file's data variables
//    for (int var_idx = 0; var_idx < dataFiles[file_idx]->nvars; var_idx++) {
//      // Loop over this variable's elements
//     // bool use4Dimensions = all_out_data[dataFiles[file_idx]->varid[var_idx]].nelem > 1;
//     // count4.at(1) = all_out_data[dataFiles[file_idx]->varid[var_idx]].nelem; // Change the number of values to write to the z dimension (array of values).
//    	int varnumelem = all_out_data[dataFiles[file_idx]->varid[var_idx]].nelem;
//      bool use4Dimensions = varnumelem > 1;
//    	// Assign the length of the data to be written for this variable
//    	int vardatasize = varnumelem * all_out_data.size();
//    	count4.at(1) = vardatasize;
//      // Create temporary array of data for this variable across all cells
//      double *vardata = new double[vardatasize];
//    	// Fill vardata
//    	for (unsigned int cellidx = 0; cellidx < all_out_data.size(); cellidx++) {
//    		vardata[cellidx*varnumelem] = all_out_data[dataFiles[file_idx]->varid[var_idx]].aggdata;
//    	}
//
//    	// Write data to file for this variable
//      try {
//      	NcVar variable = allVars.find(state->output_mapping.at(all_out_data[dataFiles[file_idx]->varid[var_idx]].varname).name)->second;
//
//        if (use4Dimensions) {
////          variable.putVar(start4, count4, all_out_data[dataFiles[file_idx]->varid[var_idx]].aggdata);
//          variable.putVar(start4, count4, vardata);
//
//        } else {
////          variable.putVar(start3, count3, all_out_data[dataFiles[file_idx]->varid[var_idx]].aggdata);
//          variable.putVar(start3, count3, vardata);
//        }
//      } catch (std::exception& e) {
//        fprintf(stderr, "Error writing variable: %s, at latIndex: %d, lonIndex: %d, timeIndex: %d\n", state->output_mapping.at(all_out_data[dataFiles[file_idx]->varid[var_idx]].varname).name.c_str(), (int)latIndex, (int)lonIndex, (int)timeIndex);
//        throw;
//      }
//    }
//  }
//}
