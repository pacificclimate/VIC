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
//  const size_t startLatIndex = latitudeToIndex(state->global_param.gridStartLat, state);
//  const size_t startLonIndex = longitudeToIndex(state->global_param.gridStartLon, state);
//
//
//  if (IS_INVALID((int)timeIndex) || IS_INVALID((int)lonIndex) || IS_INVALID((int)latIndex)) {
//    std::stringstream s;
//    s << "Error: Invalid index. timeIndex=" << timeIndex << ", lonIndex=" << lonIndex << ", latIndex=" << latIndex;
//    throw VICException(s.str());
//  }

  // Defines the dimension order of how variables are written. Only variables which have more than one value have the extra values dimension.
  // If you change the dimension ordering here, make sure that it is also changed in the WriteOutputNetCDF::initializeFile() method.
  // If the z dimension position changes also change the count4 vector update inside the nested loop below.
//  const size_t start3Vals [] = { timeIndex, startLatIndex, startLonIndex };     // (t, y, x)
//  const size_t count3Vals [] = { 1,1,1 };
//  const size_t start4Vals [] = { timeIndex, 0, startLatIndex, startLonIndex };  // (t, z, y, x)
//  const size_t count4Vals [] = { 1,1,1,1 };

  const size_t start3Vals [] = { timeIndex, 0, 0 };     // (t, y, x)
  const size_t count3Vals [] = { 1,(size_t) state->global_param.gridNumLatDivisions,(size_t) state->global_param.gridNumLonDivisions };
  const size_t start4Vals [] = { timeIndex, 0, 0, 0 };  // (t, z, y, x)
  const size_t count4Vals [] = { 1,1,(size_t) state->global_param.gridNumLatDivisions,(size_t) state->global_param.gridNumLonDivisions };
  std::vector<size_t> start3(start3Vals, start3Vals + 3), count3(count3Vals, count3Vals + 3);
  std::vector<size_t> start4(start4Vals, start4Vals + 4), count4(count4Vals, count4Vals + 4);

  std::multimap<std::string, netCDF::NcVar> allVars = netCDF->getVars();
  int num_cells = all_out_data.size();

  fprintf(stderr, "WriteOutputAllCells::write_data: all_out_data.size() = %d\n", all_out_data.size());

//  for (unsigned int cell_idx = 0; cell_idx < all_out_data.size(); cell_idx++) {
//  	for(unsigned int varnum = 0; varnum < 180; varnum++){
//  		fprintf(stderr, "WriteOutputAllCells::write_data: all_out_data[%d][%d].varname = %s\n",cell_idx, varnum, all_out_data[cell_idx][varnum].varname);
//  		fprintf(stderr, "WriteOutputAllCells::write_data: all_out_data[%d][%d].nelem = %d\n",cell_idx, varnum, all_out_data[cell_idx][varnum].nelem);
//  		fprintf(stderr, "WriteOutputAllCells::write_data: all_out_data[%d][%d].aggdata[0] = %f\n",cell_idx, varnum, all_out_data[cell_idx][varnum].aggdata[0]);
//  	}
//  }


	// Loop through (legacy) out_data_files_template for listing of output variables
	for (int file_idx = 0; file_idx < state->options.Noutfiles; file_idx++) {

		// Loop over this output file's data variables
		for (int var_idx = 0; var_idx < out_data_files_template[file_idx].nvars; var_idx++) {

      // Create temporary array of data for this variable across all cells
    	double *vardata, *vardata_ptr;
    	unsigned int varnumelem = all_out_data[0][out_data_files_template[file_idx].varid[var_idx]].nelem;
			bool use4Dimensions = varnumelem > 1;
    	size_t vardatasize = varnumelem * num_cells;
      vardata = new double [vardatasize];
      vardata_ptr = vardata;

			// Get data from each cell for this variable and put into temporary array vardata
			for (int cell_idx = 0; cell_idx < num_cells; cell_idx++) {

				// Loop over this variable's elements and append to temporary vardata array
				//vardata = all_out_data[cell_idx][out_data_files_template[file_idx].varid[var_idx]].aggdata;
				std::copy(all_out_data[cell_idx][out_data_files_template[file_idx].varid[var_idx]].aggdata + 0, \
						all_out_data[cell_idx][out_data_files_template[file_idx].varid[var_idx]].aggdata + varnumelem, \
						vardata_ptr);

				vardata_ptr += varnumelem;
			}

    	// Write data to file for this variable
      try {
    		NcVar variable = allVars.find(state->output_mapping.at(all_out_data[0][out_data_files_template[file_idx].varid[var_idx]].varname).name)->second;
    		fprintf(stderr, "WriteOutputAllCells::write_data: varname = %s\n",all_out_data[0][out_data_files_template[file_idx].varid[var_idx]].varname);

        if (use4Dimensions) {
        	count4.at(1) = vardatasize;  // Set the number of values to write to the z dimension
//        	variable.putVar(start4, count4, vardata);
        } else {
        	variable.putVar(start3, count3, vardata);
        }
      } catch (std::exception& e) {
      	//      		fprintf(stderr, "Error writing variable: %s, at timeIndex: %d\n", state->output_mapping.at(all_out_data[out_data_files_template[file_idx].varid[var_idx]].varname).name.c_str(), (int)timeIndex);
//        fprintf(stderr, "Error writing variable: %s, at timeIndex: %d\n", state->output_mapping.at(all_out_data[out_data_files_template[file_idx]->varid[var_idx]]->varname).name.c_str(), (int)timeIndex);

        fprintf(stderr, "Error writing variable: %s, at timeIndex: %d\n", state->output_mapping.at(all_out_data[out_data_files_template[file_idx].varid[var_idx]]->varname).name.c_str(), (int)timeIndex);
        throw;
      }
      delete [] vardata;
    }
  }
}
