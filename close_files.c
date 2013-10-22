#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vicNl.h"
#include <netcdf.h>
 
static char vcid[] = "$Id$";

void close_files(const filep_struct         *filep,
                 filenames_struct     *fnames,
                 bool                 compress,
                 const ProgramState   *state)
/**********************************************************************
	close_files	Dag Lohmann		January 1996

  This routine closes all forcing data files, and output files.

  Modifications:
  7-19-96  Files are now gzipped when they are closed.  This
	   was added to save space when using large volumes
	   of data.						KAC
  02-27-01 Now closes files opened for lake model applications  KAC
  11-18-02 Now closes lake debugging file.                      LCB
  29-Oct-03 Distinguishing between input lakeparam file and output
	    lake file.						TJB
  2005-Mar-24 Added support for ALMA output files.		TJB
  2005-Apr-10 Added logic for OUTPUT_FORCE option.		TJB
  2006-Sep-23 Implemented flexible output configuration; uses new
	      out_data_files structure. TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct. TJB

**********************************************************************/
{
  int filenum;

  /**********************
    Close All Input Files
    **********************/

  if(state->param_set.FORCE_FORMAT[0] != NETCDF)
    fclose(filep->forcing[0]);
  else /* FIXME this should only happen ONCE per invokation of VIC! (below too) */
    nc_close(filep->forcing_ncid[0]);
  if(compress) compress_files(fnames->forcing[0]);
  if(filep->forcing[1]!=NULL) {
    if(state->param_set.FORCE_FORMAT[1] != NETCDF)
      fclose(filep->forcing[1]);
    else
      nc_close(filep->forcing_ncid[1]);      
    if(compress) compress_files(fnames->forcing[1]);
  }

#if !OUTPUT_FORCE

  /*******************************
    Close All Used Debugging Files
    *******************************/ 

#if LINK_DEBUG
  if(state->debug.DEBUG || state->debug.PRT_TEMP) {
    fclose(state->debug.fg_temp);
  }
  if(state->debug.DEBUG || state->debug.PRT_MOIST) {
    fclose(state->debug.fg_moist);
  }
  if(state->debug.DEBUG || state->debug.PRT_KAPPA) {
    fclose(state->debug.fg_kappa);
  }
  if(state->debug.DEBUG || state->debug.PRT_LAKE) {
    fclose(state->debug.fg_lake);
  }
  if(state->debug.DEBUG || state->debug.PRT_BALANCE) {
    fclose(state->debug.fg_balance);
  }
  if(state->debug.DEBUG || state->debug.PRT_FLUX) {
    fclose(state->debug.fg_energy);
  }
  if(state->debug.DEBUG || state->debug.PRT_SNOW) {
    fclose(state->debug.fg_snow);
  }
  if(state->debug.DEBUG || state->debug.PRT_GRID) {
    fclose(state->debug.fg_grid);
  }
#endif /* LINK_DEBUG */
#endif /* !OUTPUT_FORCE */

}
