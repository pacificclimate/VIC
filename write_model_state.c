#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vicNl.h"
#include "StateIOContext.h"

static char vcid[] = "$Id$";

void write_model_state(cell_info_struct* cell, const char* filename, const ProgramState  *state)
/*********************************************************************
  write_model_state      Keith Cherkauer           April 14, 2000

  This subroutine saves the model state at hour 0 of the date 
  defined in the global control file using STATEDAY, STATEMONTH,
  and STATEYEAR.  The saved files can then be used to initialize 
  the model to the same state as when the files were created.

  Soil moisture, soil thermal, and snowpack variables  are stored 
  for each vegetation type and snow band.  However moisture variables
  from the distributed precipitation model are averaged so that the
  model is restarted with mu = 1.

  Modifications:
  04-10-03 Rewritten to handle updates to vicNl_def.h and to write
           the file as binary to minimize write time and differences
           with simulations started with the state file.         KAC
  04-10-03 Model is now restarted with the correct values for mu
           and LAST_STORM
  06-03-03 Modified to create ASCII as well as BINARY state file.  KAC
  06-06-03 It is not necessary to store the current ice content as
           it is recomputed in initialize_model_state.         KAC
  09-Oct-03 Removed initial space on veg/band info line in ASCII file.		TJB
  26-Oct-04 Changed calculation of Nbytes in binary state file to
	    account for bare soil values (extra veg class per grid
	    cell).  Without this fix, attempts to skip grid cells
	    fail.								TJB
  01-Nov-04 Added storage of state variables for SPATIAL_FROST and
	    LAKE_MODEL.								TJB
  02-Nov-04 Added a few more lake state variables.				TJB
  03-Nov-04 Now outputs extra_veg to aid other programs in parsing
	    state files.							TJB
  2005-Dec-07 STATE_FILE option is set in global file.				GCT
  2005-Jan-10 writes temp[0] instead of tp_in for lake skin surface
	      temperature.							JCA
  2005-Jan-10 modified to write lake nodal variables for each of the
	      active nodes.							JCA
  2006-Jun-16 Skip writing snow band if areafract < 0.				GCT
  2006-Aug-23 Changed order of fread/fwrite statements from ...1, sizeof...
	      to ...sizeof, 1,...						GCT
  2006-Sep-07 Changed "Skip writing snow band if areafract < 0" to "<=0".	GCT
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.					TJB
  2007-Apr-24 Modified to write Zsum_node.					JCA
  2007-Apr-25 Removed variable Nsum.						JCA
  2007-Aug-24 Added features for EXCESS_ICE option.				JCA
  2007-Nov-06 New list of lake state variables.					LCB via TJB
  2009-Jul-31 Removed extra lake/wetland veg tile; updated set of lake state
	      variables.							TJB
  2009-Aug-27 Now once again writes data for all bands, regardless of
	      whether they have area > 0.  This makes it much easier to ensure
	      that the value of Nbands stored in the state file matches the number
	      of bands actually stored in the state file.			TJB
  2009-Sep-28 Now stores soil, snow, and energy states from lake separately
	      from wetland.							TJB
  2010-Mar-05 Fixed typo in writing of number of lake active nodes.		TJB
  2012-Jan-01 Removed lake area condition from logic determining whether to write
	      lake state data.  Now, if options.LAKES is TRUE, every grid cell
	      will save lake state data.  If no lake is present, default NULL
	      values will be stored.						TJB
*********************************************************************/
{
  int    Ndist;
  int    Nbands;

  if(state->options.DIST_PRCP)
    Ndist = 2;
  else 
    Ndist = 1;
  Nbands = state->options.SNOW_BAND;

  StateIOContext context(filename, state);
  StateIO* writer = context.stream;

  /* write cell information */
  writer->write(&cell->soil_con.gridcel, 1, NULL);
  writer->write(&cell->veg_con[0].vegetat_type_num, 1, NULL);
  writer->write(&Nbands, 1, NULL);
  
  /* Write soil thermal node deltas */
  writer->write(cell->soil_con.dz_node, state->options.Nnode, NULL);

  /* Write soil thermal node depths */
  writer->write(cell->soil_con.Zsum_node, state->options.Nnode, NULL);

  writer->writeNewline();
  
  /* Write dynamic soil properties */
#if EXCESS_ICE
  /* Write soil depth */
  writer->write(cell->soil_con.depth, state->options.Nlayer, NULL);
  
  /* Write effective porosity */
  writer->write(cell->soil_con.effective_porosity, state->options.Nlayer, NULL);
  
  /* Write damping depth */
  writer->write(&cell->soil_con.dp, 1, NULL);
  writer->writeNewline();
#endif
  
  /* Output for all vegetation types */
  for (std::vector<HRU>::iterator it = cell->prcp.hruList.begin(); it != cell->prcp.hruList.end(); ++it) {
    
    // Do the following only once per vegetation type
    if (it->bandIndex == 0) {
      // Store distributed precipitation fraction
      writer->write(&cell->prcp.mu[it->vegIndex], 1, NULL);

      // Store distributed precipitation variables
      writer->write(&cell->init_STILL_STORM[it->vegIndex], 1, NULL);
      writer->write(&cell->init_DRY_TIME[it->vegIndex], 1, NULL);
      writer->writeNewline();
    }

    /* Output for all snow bands */
    /* Write cell identification information */
      writer->write(&it->vegIndex, 1, NULL);
      writer->write(&it->bandIndex, 1, NULL);

    for (int dist = 0; dist < Ndist; dist++) {
      hru_data_struct& cellRef = it->cell[dist];
      /* Write total soil moisture */
      std::vector<double> soilMoisture;
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        double tmpval = cellRef.layer[lidx].moist; /* MPN */
        soilMoisture.push_back(tmpval);
      }
      writer->write(&(soilMoisture[0]), soilMoisture.size(), NULL);

      /* Write average ice content */
      std::vector<double> avgIceContent;
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
#if SPATIAL_FROST
#error
        for (int frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
          tmpval = cellRef.layer[lidx].soil_ice[frost_area];
          avgIceContent.push_back(tmpval);
        }
#else
        double tmpval = cellRef.layer[lidx].soil_ice;
        avgIceContent.push_back(tmpval);
#endif // SPATIAL_FROST
      }
      writer->write(&(avgIceContent[0]), avgIceContent.size(), NULL);

      /* Write dew storage */
      if (it->vegIndex < cell->veg_con[0].vegetat_type_num) {
        double tmpval = it->veg_var[dist].Wdew;
        writer->write(&tmpval, 1, NULL);
      }
    }

    /* Write snow data */
    writer->write(&it->snow.last_snow, 1, NULL);
    writer->write(&it->snow.MELTING, 1, NULL);
    writer->write(&it->snow.coverage, 1, NULL);
    writer->write(&it->snow.swq, 1, NULL);
    writer->write(&it->snow.surf_temp, 1, NULL);
    writer->write(&it->snow.surf_water, 1, NULL);
    writer->write(&it->snow.pack_temp, 1, NULL);
    writer->write(&it->snow.pack_water, 1, NULL);
    writer->write(&it->snow.density, 1, NULL);
    writer->write(&it->snow.coldcontent, 1, NULL);
    writer->write(&it->snow.snow_canopy, 1, NULL);

    /* Write soil thermal node temperatures */
    writer->write(it->energy.T, state->options.Nnode, NULL);

    writer->writeNewline();
  }

  if (state->options.LAKES) {
    for (int dist = 0; dist < Ndist; dist++) {
      // Store both wet and dry fractions if using distributed precipitation

      /* Write total soil moisture */
      std::vector<double> moistValues;  // Convert to array for single write operation.
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        moistValues.push_back(cell->prcp.lake_var.soil.layer[lidx].moist);
      }
      writer->write(&(moistValues[0]), moistValues.size(), NULL);


      /* Write average ice content */
      std::vector<double> avgIceContent;
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
#if SPATIAL_FROST
        for (int frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
          avgIceContent.push_back(cell->prcp.lake_var.soil.layer[lidx].soil_ice[frost_area]);
        }
#else
        avgIceContent.push_back(cell->prcp.lake_var.soil.layer[lidx].soil_ice);
#endif // SPATIAL_FROST
      }
      writer->write(&(avgIceContent[0]), avgIceContent.size(), NULL);
    }

    /* Write snow data */
    writer->write(&cell->prcp.lake_var.snow.last_snow,   1, NULL);
    writer->write(&cell->prcp.lake_var.snow.MELTING,     1, NULL);
    writer->write(&cell->prcp.lake_var.snow.coverage,    1, NULL);
    writer->write(&cell->prcp.lake_var.snow.swq,         1, NULL);
    writer->write(&cell->prcp.lake_var.snow.surf_temp,   1, NULL);
    writer->write(&cell->prcp.lake_var.snow.surf_water,  1, NULL);
    writer->write(&cell->prcp.lake_var.snow.pack_temp,   1, NULL);
    writer->write(&cell->prcp.lake_var.snow.pack_water,  1, NULL);
    writer->write(&cell->prcp.lake_var.snow.density,     1, NULL);
    writer->write(&cell->prcp.lake_var.snow.coldcontent, 1, NULL);
    writer->write(&cell->prcp.lake_var.snow.snow_canopy, 1, NULL);

    /* Write soil thermal node temperatures */
    writer->write(cell->prcp.lake_var.energy.T, state->options.Nnode, NULL);

    /* Write lake-specific variables */
    writer->write(&cell->prcp.lake_var.activenod, 1, NULL);
    writer->write(&cell->prcp.lake_var.dz, 1, NULL);
    writer->write(&cell->prcp.lake_var.surfdz, 1, NULL);
    writer->write(&cell->prcp.lake_var.ldepth, 1, NULL);
    writer->write(cell->prcp.lake_var.surface, cell->prcp.lake_var.activenod, NULL);
    writer->write(&cell->prcp.lake_var.sarea, 1, NULL);
    writer->write(&cell->prcp.lake_var.volume, 1, NULL);
    writer->write(cell->prcp.lake_var.temp, cell->prcp.lake_var.activenod, NULL);
    writer->write(&cell->prcp.lake_var.tempavg, 1, NULL);
    writer->write(&cell->prcp.lake_var.areai, 1, NULL);
    writer->write(&cell->prcp.lake_var.new_ice_area, 1, NULL);
    writer->write(&cell->prcp.lake_var.ice_water_eq, 1, NULL);
    writer->write(&cell->prcp.lake_var.hice, 1, NULL);
    writer->write(&cell->prcp.lake_var.tempi, 1, NULL);
    writer->write(&cell->prcp.lake_var.swe, 1, NULL);
    writer->write(&cell->prcp.lake_var.surf_temp, 1, NULL);
    writer->write(&cell->prcp.lake_var.pack_temp, 1, NULL);
    writer->write(&cell->prcp.lake_var.coldcontent, 1, NULL);
    writer->write(&cell->prcp.lake_var.surf_water, 1, NULL);
    writer->write(&cell->prcp.lake_var.pack_water, 1, NULL);
    writer->write(&cell->prcp.lake_var.SAlbedo, 1, NULL);
    writer->write(&cell->prcp.lake_var.sdepth, 1, NULL);

    writer->writeNewline();

  }
  /* Force file to be written */
  writer->flush();
}

