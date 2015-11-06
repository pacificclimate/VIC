#include "WriteOutputAllCells.h"

#include <netcdf>
#include <ctime>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace netCDF;

WriteOutputAllCells::WriteOutputAllCells(const ProgramState* state) : WriteOutputNetCDF(state) {

}

WriteOutputAllCells::~WriteOutputAllCells() {

}

void WriteOutputAllCells::write_data(std::vector<out_data_struct*>& all_out_data, out_data_file_struct *out_data_files_template, const dmy_struct* dmy, int dt, const ProgramState* state) {

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

  fprintf(stderr, "WriteOutputAllCells::write_data: all_out_data.size() = %d\n", all_out_data.size());

  for (unsigned int cell_idx = 0; cell_idx < all_out_data.size(); cell_idx++) {
  	for(unsigned int varnum = 0; varnum < 45; varnum++){
//  		all_out_data[cell_idx]->
  		fprintf(stderr, "WriteOutputAllCells::write_data: all_out_data[%d][%d].varname = %s\n",cell_idx, varnum, all_out_data[cell_idx][varnum].varname);
  		fprintf(stderr, "WriteOutputAllCells::write_data: all_out_data[%d][%d].nelem = %d\n",cell_idx, varnum, all_out_data[cell_idx][varnum].nelem);
  		fprintf(stderr, "WriteOutputAllCells::write_data: all_out_data[%d][%d].data[0] = %f\n",cell_idx, varnum, all_out_data[cell_idx][varnum].data[0]);
  	}
  }

//  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {
//  for (unsigned int file_idx = 0; file_idx < out_data_files_template.size(); file_idx++) {
  for (unsigned int file_idx = 0; file_idx < 2; file_idx++) {

    // Loop over this output file's data variables
//    for (int var_idx = 0; var_idx < out_data_files_template[file_idx]->nvars; var_idx++) {
    for (int var_idx = 0; var_idx < out_data_files_template[file_idx].nvars; var_idx++) {

      // Loop over this variable's elements
     // bool use4Dimensions = all_out_data[dataFiles[file_idx]->varid[var_idx]].nelem > 1;
     // count4.at(1) = all_out_data[dataFiles[file_idx]->varid[var_idx]].nelem; // Change the number of values to write to the z dimension (array of values).
//    	unsigned int varnumelem = all_out_data[out_data_files_template[file_idx]->varid[var_idx]]->nelem;
    	unsigned int varnumelem = all_out_data[out_data_files_template[file_idx].varid[var_idx]]->nelem;

      bool use4Dimensions = varnumelem > 1;
    	// Assign the length of the data to be written for this variable
      unsigned int num_cells = all_out_data.size();

    	size_t vardatasize = varnumelem * num_cells;
    	count4.at(1) = vardatasize;
      // Create temporary array of data for this variable across all cells
    	double *vardata;
      vardata = new double [vardatasize];
    	// Fill vardata
    	for (unsigned int cellidx = 0; cellidx < num_cells; cellidx++) {
//    		vardata[cellidx*varnumelem] = *(all_out_data[out_data_files_template[file_idx]->varid[var_idx]]->aggdata);
//    		vardata[cellidx*varnumelem] = *all_out_data[out_data_files_template[file_idx].varid[var_idx]].aggdata;
//    		vardata[cellidx*varnumelem] = all_out_data[out_data_files_template[file_idx]->varid[var_idx]]->aggdata;

//    		vardata[cellidx*varnumelem] = all_out_data[out_data_files_template[file_idx].varid[var_idx]]->aggdata;

    	}

    	// Write data to file for this variable
      try {
      	//      		NcVar variable = allVars.find(state->output_mapping.at(all_out_data[out_data_files_template[file_idx].varid[var_idx]].varname).name)->second;

//      		NcVar variable = allVars.find(state->output_mapping.at(all_out_data[out_data_files_template[file_idx]->varid[var_idx]]->varname).name)->second;
    		NcVar variable = allVars.find(state->output_mapping.at(all_out_data[out_data_files_template[file_idx].varid[var_idx]]->varname).name)->second;

        if (use4Dimensions) {
        	// NOTE: commented out temporarily for test:
  //        variable.putVar(start4, count4, vardata);

        } else {
        	// NOTE: commented out temporarily for test:
  //        variable.putVar(start3, count3, vardata);
        }
      } catch (std::exception& e) {
      	//      		fprintf(stderr, "Error writing variable: %s, at timeIndex: %d\n", state->output_mapping.at(all_out_data[out_data_files_template[file_idx].varid[var_idx]].varname).name.c_str(), (int)timeIndex);
//        fprintf(stderr, "Error writing variable: %s, at timeIndex: %d\n", state->output_mapping.at(all_out_data[out_data_files_template[file_idx]->varid[var_idx]]->varname).name.c_str(), (int)timeIndex);

        fprintf(stderr, "Error writing variable: %s, at timeIndex: %d\n", state->output_mapping.at(all_out_data[out_data_files_template[file_idx].varid[var_idx]]->varname).name.c_str(), (int)timeIndex);
        throw;
      }
    }
  }
}
