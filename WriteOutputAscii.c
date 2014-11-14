#include "WriteOutputAscii.h"

void WriteOutputAscii::openFile() {
  for (unsigned int i = 0; i < dataFiles.size(); i++) {
    dataFiles[i]->fh = open_file(dataFiles[i]->filename, "w");
  }
}


void WriteOutputAscii::write_data(out_data_struct* out_data, const dmy_struct* dmy, int dt,
    const ProgramState* state) {
  // Loop over output files
  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {

  if (!state->options.OUTPUT_FORCE) {
    // Write the date
    if (dt < 24) {
      // Write year, month, day, and hour
      fprintf(dataFiles[file_idx]->fh, "%04i\t%02i\t%02i\t%02i\t",
              dmy->year, dmy->month, dmy->day, dmy->hour);
    }
    else {
      // Only write year, month, and day
      fprintf(dataFiles[file_idx]->fh, "%04i\t%02i\t%02i\t",
              dmy->year, dmy->month, dmy->day);
    }
  }

    // Loop over this output file's data variables
    for (int var_idx = 0; var_idx < dataFiles[file_idx]->nvars; var_idx++) {
      // Loop over this variable's elements
      for (int elem_idx = 0; elem_idx < out_data[dataFiles[file_idx]->varid[var_idx]].nelem; elem_idx++) {
        if (!(var_idx == 0 && elem_idx == 0)) {
          fprintf(dataFiles[file_idx]->fh, "\t ");
        }
        fprintf(dataFiles[file_idx]->fh, out_data[dataFiles[file_idx]->varid[var_idx]].format, out_data[dataFiles[file_idx]->varid[var_idx]].aggdata[elem_idx]);
      }
    }
    fprintf(dataFiles[file_idx]->fh, "\n");
  }
}

const char* WriteOutputAscii::getDescriptionOfOutputType() {
  return "ASCII";
}

void WriteOutputAscii::compressFiles() {
  for (unsigned int i = 0; i < dataFiles.size(); i++) {
    compress_files(dataFiles[i]->filename);
  }
}

/* ASCII header format:

 # NRECS: (nrecs)
 # DT: (dt)
 # STARTDATE: yyyy-mm-dd hh:00:00
 # ALMA_OUTPUT: (0 or 1)
 # NVARS: (Nvars)
 # VARNAME    VARNAME   VARNAME   ...

 where
    nrecs       = Number of records in the file
    dt          = Output time step length in hours
    start date  = Date and time of first record of file
    ALMA_OUTPUT = Indicates units of the variables; 0 = standard VIC units; 1 = ALMA units
    Nvars       = Number of variables in the file, including date fields*/
void WriteOutputAscii::write_header(out_data_struct* out_data, const dmy_struct* dmy, const ProgramState* state) {
  char                tmp_ALMA_OUTPUT;
  char                Nvars;

  if (state->options.ALMA_OUTPUT)
    tmp_ALMA_OUTPUT = 1;
  else
    tmp_ALMA_OUTPUT = 0;

  // Loop over output files
  for (unsigned int file_idx = 0; file_idx < dataFiles.size(); file_idx++) {

    // Header part 1: Global attributes
    Nvars = dataFiles[file_idx]->nvars;
    if (!state->options.OUTPUT_FORCE) {
    	if (state->global_param.out_dt < 24)
    		Nvars += 4;
    	else
    		Nvars += 3;
    }
    fprintf(dataFiles[file_idx]->fh, "# NRECS: %d\n", state->global_param.nrecs);
    fprintf(dataFiles[file_idx]->fh, "# DT: %d\n", state->global_param.out_dt);
    fprintf(dataFiles[file_idx]->fh, "# STARTDATE: %04d-%02d-%02d %02d:00:00\n",
      dmy->year, dmy->month, dmy->day, dmy->hour);
    fprintf(dataFiles[file_idx]->fh, "# ALMA_OUTPUT: %d\n", tmp_ALMA_OUTPUT);
    fprintf(dataFiles[file_idx]->fh, "# NVARS: %d\n", Nvars);

    // Header part 2: Variables
    fprintf(dataFiles[file_idx]->fh, "# ");

    if (!state->options.OUTPUT_FORCE) {
    	// Write the date
    	if (state->global_param.out_dt < 24) {
    		// Write year, month, day, and hour
    		fprintf(dataFiles[file_idx]->fh, "YEAR\tMONTH\tDAY\tHOUR\t");
    	}
    	else {
    		// Only write year, month, and day
    		fprintf(dataFiles[file_idx]->fh, "YEAR\tMONTH\tDAY\t");
    	}
    }

    // Loop over this output file's data variables
    for (int var_idx = 0; var_idx < dataFiles[file_idx]->nvars; var_idx++) {
      // Loop over this variable's elements
      for (int elem_idx = 0; elem_idx < out_data[dataFiles[file_idx]->varid[var_idx]].nelem; elem_idx++) {
        if (!(var_idx == 0 && elem_idx == 0)) {
          fprintf(dataFiles[file_idx]->fh, "\t ");
        }
        fprintf(dataFiles[file_idx]->fh, "%s", out_data[dataFiles[file_idx]->varid[var_idx]].varname);
        if (out_data[dataFiles[file_idx]->varid[var_idx]].nelem > 1) {
          fprintf(dataFiles[file_idx]->fh, "_%d", elem_idx);
        }
      }
    }
    fprintf(dataFiles[file_idx]->fh, "\n");

  }

}


