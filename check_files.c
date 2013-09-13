#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

filep_struct get_files(const filenames_struct *fnames)
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
  extern option_struct  options;

  filep_struct file_pointers;

  file_pointers.soilparam   = open_file(fnames->soil, "r");
#if !OUTPUT_FORCE
  file_pointers.veglib      = open_file(fnames->veglib, "r");
  file_pointers.vegparam    = open_file(fnames->veg, "r");
  if(options.SNOW_BAND>1)
    file_pointers.snowband    = open_file(fnames->snowband, "r");
  if ( options.LAKES )
    file_pointers.lakeparam = open_file(fnames->lakeparam,"r");
#endif /* !OUTPUT_FORCE */

  return file_pointers;
}


