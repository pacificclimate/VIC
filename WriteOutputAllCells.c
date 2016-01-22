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

void WriteOutputAllCells::write_data(std::vector<out_data_struct*>& all_out_data, out_data_file_struct *out_data_files_template, const dmy_struct *dmy, int dt, const ProgramState* state) {

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

			  		fprintf(stderr, "WriteOutputAllCells::write_data: varname = %s,  timeIndex = %d\n",all_out_data[0][out_data_files_template[file_idx].varid[var_idx]].varname, timeIndex);

			// Create temporary array of data for this variable across all cells
			float *vardata, *vardata_ptr;
			int varnumelem = all_out_data[0][out_data_files_template[file_idx].varid[var_idx]].nelem;
			bool use4Dimensions = varnumelem > 1;
			int vardatasize = varnumelem * num_cells;
			int modeled_cell_idx = 0; // index among modeled cells

			vardata = new float [vardatasize];
			vardata_ptr = vardata;
			modeled_cell_mask_ptr = state->modeled_cell_mask;

			// Interleave data from each cell for this variable in temporary array vardata
			for (int elem=0; elem<varnumelem; elem++) {
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

