#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include "WriteOutputNetCDF.h"

static char vcid[] = "$Id$";

void initializeNetCDFOutput(const filenames_struct *fnames, const out_data_file_struct* outFiles, ProgramState *state) {
  // Initialise the netcdf full path name.
  strcpy(state->options.NETCDF_FULL_FILE_PATH, fnames->result_dir);
  strcat(state->options.NETCDF_FULL_FILE_PATH, "/");
  strcat(state->options.NETCDF_FULL_FILE_PATH, fnames->netCDFOutputFileName);

  if (state->options.OUTPUT_FORMAT == OutputFormat::NETCDF_FORMAT) {
#if NETCDF_OUTPUT_AVAILABLE
    WriteOutputNetCDF output(state);
    copy_data_file_format(outFiles, output.dataFiles, state);
    output.initializeFile(state);  // This is only done once per invocation of VIC. It creates a fresh netcdf output file.
#endif /* NETCDF_OUTPUT_AVAILABLE */
  }
}

filep_struct get_files(const filenames_struct *fnames, ProgramState* state)
/**********************************************************************
	check_files		Dag Lohmann		January 1996

  This routine opens files for soil, vegetation, and global parameters.

  Modifcations:
  02-27-01 Added controls for lake model parameter file    KAC
  2005-Apr-13 Added logic for OUTPUT_FORCE option.			TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct.	TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB

**********************************************************************/
{
  filep_struct file_pointers;

  file_pointers.soilparam   = open_file(fnames->soil, "r");

  if (!state->options.OUTPUT_FORCE) {
    file_pointers.veglib      = open_file(fnames->veglib, "r");
    file_pointers.vegparam    = open_file(fnames->veg, "r");
    if(state->options.SNOW_BAND>1)
      file_pointers.snowband    = open_file(fnames->snowband, "r");
    if ( state->options.LAKES )
      file_pointers.lakeparam = open_file(fnames->lakeparam,"r");
  }

  return file_pointers;
}


