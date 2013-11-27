#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include "vicNl.h"
#include "StateIOContext.h"

static char vcid[] = "$Id$";

void read_initial_model_state(const char* initStateFilename, cell_info_struct *cell, int Nveg, int Ndist, const ProgramState *state)
/*********************************************************************
  read_initial_model_state   Keith Cherkauer         April 14, 2000

  This subroutine initializes the model state at hour 0 of the date 
  defined in the given state file.  

  Soil moisture, soil thermal, and snowpack variables  are stored 
  for each vegetation type and snow band.  However moisture variables
  from the distributed precipitation model are averaged so that the
  model is restarted with mu = 1.

  Modifications:
  04-10-03 Rewritten to handle updates to vicNl_def.h and to write
           the file as binary to minimize write time and differences
           with simulations started with the state file.         KAC
  04-10-03 Modified to read storm parameters from the state file.  KAC
  06-03-03 Modified to read ASCII as well as BINARY state file.  KAC
  06-06-03 It is not necessary to initialize ice content from the
           model state file as the model recomutes it in subsequent
           steps.                                               KAC
  06-06-03 Added check to make sure that soil moisture obtained from
           the state file does not exceed the maximum moisture 
           content.                                             KAC
  06-07-03 Added check to verify that the sum of the defined nodes
           equals the damping depth.                            KAC
  03-Oct-03 Modified to loop over tmp_Nveg and tmp_Nband when searching
            for desired cellnum in ASCII file, rather than over Nveg
            and Nbands.  As we skip over other records in the state
            file while searching for the desired record, the loop
            must parse each undesired record differently, according
            to how many veg classes and snow bands exist in the
            record (tmp_Nveg and tmp_Nband, respectively), rather
            than the number of veg classes and snow bands in the
            desired record (Nveg and Nbands, respectively).			TJB
  01-Nov-04 Modified to read state files containing SPATIAL_FROST
	    and LAKE_MODEL state variables.					TJB
  02-Nov-04 Added a few more lake state variables.				TJB
  03-Nov-04 Now reads extra_veg from state file.				TJB
  2005-Apr-10 Fixed incorrect check on soil node depths.			TJB
  2005-Jan-10 modified to read lake nodal variables for each of the
	      active nodes.							JCA
  2006-Jun-16 Skip reading if areafract < 0.					GCT
  2006-Aug-23 Changed order of fread/fwrite statements from ...1, sizeof...
	      to ...sizeof, 1,...						GCT
  2006-Sep-07 Changed "Skip reading if areafract < 0" to "<=0".			GCT
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      This included moving infiles.statefile to filep.init_state.	TJB
  2006-Nov-07 Removed LAKE_MODEL option.					TJB
  2007-Apr-28 modified to read Zsum_node.					JCA
  2007-May-07 Fixed fread checks to make sure correct number of items were
	      read in rather than the size of the item read in.			JCA
  2007-May-07 Nsum and sum removed from declaration.				JCA
  2007-Aug-24 Added features for EXCESS_ICE option.				JCA
  2007-Sep-14 Fixed bug for read-in during EXCESS_ICE option.			JCA
  2007-Sep-18 Check for soil moist exceeding max moist moved from
	      here to initialize_model_state.					JCA
  2007-Nov-06 New list of lake state variables.					LCB via TJB
  2009-Jul-31 Removed extra lake/wetland veg tile; updated set of lake state
	      variables.							TJB
  2009-Aug-27 Now once again expects to read data for all bands, regardless of
	      whether they have area > 0.  This makes it much easier to ensure
	      that the value of Nbands stored in the state file matches the number
	      of bands actually stored in the state file.			TJB
  2009-Sep-28 Now stores soil, snow, and energy states from lake separately
	      from wetland.							TJB
  2010-Jan-10 Corrected typo in condition for checking Wdew.			TJB
  2012-Jan-01 Removed lake area condition from logic determining whether to read
	      lake state data.  Now, if options.LAKES is TRUE, every grid cell
	      will save lake state data.  If no lake is present, default NULL
	      values will be stored.						TJB
*********************************************************************/
{
  char   tmpstr[MAXSTRING];
  char   ErrStr[MAXSTRING];
  char   tmpchar;
  double tmpval;
  double depth_node[MAX_NODES];
  int    tmp_cellnum;
  int    tmp_Nveg;
  int    tmp_Nband;
  int    tmp_char;
  int    byte, Nbytes;
  int    tmp_int, node;
  int    frost_area;
  
  StateIOContext context(initStateFilename, StateIO::Reader, state);
  StateIO* reader = context.stream;

#if !NO_REWIND 
  reader->rewindFile();
  /* skip header */
  StateHeader header = reader->readHeader();

#else
#error // NO_REWIND is an untested code path. Continue at your own risk!
#endif

  int cellNVeg, cellNBand;
  if (reader->seekToCell(cell->soil_con.gridcel, &cellNVeg, &cellNBand) < 0) {
    std::stringstream ss;
    ss << "Requested grid cell (" << cell->soil_con.gridcel << ") is not in the model state file.";
    throw VICException(ss.str());
  }

  if (cellNVeg != Nveg) {
    std::stringstream ss;
    ss << "The number of vegetation types in cell " << cell->soil_con.gridcel << " (" << cellNVeg << ") does not equal that defined in vegetation parameter file (" << Nveg << ").  Check your input files.";
    throw VICException(ss.str());
  }
  if (cellNBand != state->options.SNOW_BAND) {
    std::stringstream ss;
    ss << "The number of snow bands in cell " << cell->soil_con.gridcel << " (" << cellNBand << ") does not equal that defined in the snow band file (" << state->options.SNOW_BAND << ").  Check your input files.";
    throw VICException(ss.str());
  }

  processCellForStateFile(cell, reader, state);

  //Read specific stuff follows:
  if ( state->options.Nnode == 1 ) cell->soil_con.dz_node[0] = 0;

  if ( state->options.Nnode == 1 ) cell->soil_con.Zsum_node[0] = 0;
  if ( cell->soil_con.Zsum_node[state->options.Nnode-1] - cell->soil_con.dp > SMALL) {
    fprintf( stderr, "WARNING: Sum of soil nodes (%f) exceeds defined damping depth (%f).  Resetting damping depth.\n", cell->soil_con.Zsum_node[state->options.Nnode-1], cell->soil_con.dp );
    cell->soil_con.dp = cell->soil_con.Zsum_node[state->options.Nnode-1];
  }

  for (std::vector<HRU>::iterator element = cell->prcp.hruList.begin(); element != cell->prcp.hruList.end(); ++element) {
    if (element->snow.density > 0.)
      element->snow.depth = 1000. * element->snow.swq / element->snow.density;
  }

  if (state->options.LAKES == TRUE) {
    if (cell->prcp.lake_var.snow.density > 0.)
      cell->prcp.lake_var.snow.depth = 1000. * cell->prcp.lake_var.snow.swq / cell->prcp.lake_var.snow.density;
  }

}
