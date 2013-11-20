#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>

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
  int    Nbands;

  Nbands = state->options.SNOW_BAND;

  StateIOContext context(filename, StateIO::Writer, state);
  StateIO* writer = context.stream;

  /* write cell information */
  writer->write(&cell->soil_con.gridcel, 1, StateVariables::GRID_CELL);
  writer->write(&cell->veg_con[0].vegetat_type_num, 1, StateVariables::VEG_TYPE_NUM);
  writer->write(&Nbands, 1, StateVariables::NUM_BANDS);

  processCellForStateFile(cell, writer, state);
  
}

void processCellForStateFile(cell_info_struct* cell, StateIO* stream, const ProgramState *state) {

  using namespace StateVariables;

  int Ndist;
  if(state->options.DIST_PRCP)
    Ndist = 2;
  else
    Ndist = 1;

  /* Write soil thermal node deltas */
  stream->process(cell->soil_con.dz_node, state->options.Nnode, SOIL_DZ_NODE);

  /* Write soil thermal node depths */
  stream->process(cell->soil_con.Zsum_node, state->options.Nnode, SOIL_ZSUM_NODE);

  stream->processNewline();
  
  /* Write dynamic soil properties */
#if EXCESS_ICE
  /* Write soil depth */
  stream->process(cell->soil_con.depth, state->options.Nlayer, SOIL_DEPTH);
  
  /* Write effective porosity */
  stream->process(cell->soil_con.effective_porosity, state->options.Nlayer, SOIL_EFFECTIVE_POROSITY);
  
  /* Write damping depth */
  stream->process(&cell->soil_con.dp, 1, SOIL_DP);
  stream->processNewline();
#endif
  
  /* Output for all vegetation types */
  for (std::vector<HRU>::iterator it = cell->prcp.hruList.begin(); it != cell->prcp.hruList.end(); ++it) {
    
    // Do the following only once per vegetation type
    if (it->bandIndex == 0) {
      // Store distributed precipitation fraction
      stream->process(&cell->prcp.mu[it->vegIndex], 1, PRCP_MU);

      // Store distributed precipitation variables
      stream->process(&cell->init_STILL_STORM[it->vegIndex], 1, INIT_STILL_STORM);
      stream->process(&cell->init_DRY_TIME[it->vegIndex], 1, INIT_DRY_TIME);
      stream->processNewline();
    }

    /* Output for all snow bands */
    /* Write cell identification information */
    int originalVeg = it->vegIndex;
    int originalBand = it->bandIndex;

    stream->process(&it->vegIndex, 1, HRU_VEG_INDEX);
    stream->process(&it->bandIndex, 1, HRU_BAND_INDEX);

    // The following is read specific and there is nothing we can do about it.
    if (stream->getType() == StateIO::Reader) {
      if (originalVeg != it->vegIndex || originalBand != it->bandIndex) {
        std::stringstream ss;
        ss << "The vegetation and snow band indices in the model state file (veg = " << originalVeg << ", band = " << originalBand << ")";
        ss << "do not match those currently requested (veg = " << it->vegIndex << ", band = " << it->bandIndex << "). ";
        ss << "Model state file must be stored with variables for all vegetation indexed by variables for all snow bands.\n";
        throw new VICException(ss.str());
      }
    }

    for (int dist = 0; dist < Ndist; dist++) {
      hru_data_struct& cellRef = it->cell[dist];
      /* Write total soil moisture */
      double soilMoisture [state->options.Nlayer];
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        soilMoisture[lidx] = (cellRef.layer[lidx].moist);   // This is write specific (no-op for read).
      }
      stream->process(soilMoisture, state->options.Nlayer, LAYER_MOIST);
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        cellRef.layer[lidx].moist = soilMoisture[lidx];   // This is read specific (no-op for write).
      }

      /* Write average ice content */

#if SPATIAL_FROST
#error
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        stream->process(cellRef.layer[lidx].soil_ice, FROST_SUBAREAS, LAYER_SOIL_ICE);
      }
#else
      double iceContent [state->options.Nlayer];
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        iceContent[lidx] = cellRef.layer[lidx].soil_ice;        // Write specific.
      }
      stream->process(iceContent, state->options.Nlayer, LAYER_ICE_CONTENT);
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        cellRef.layer[lidx].soil_ice = iceContent[lidx];        // Read specific.
      }
#endif // SPATIAL_FROST

      /* Write dew storage */
      if (it->vegIndex < cell->veg_con[0].vegetat_type_num) {
        stream->process(&it->veg_var[dist].Wdew, 1, HRU_VEG_VAR_WDEW);
      }
    }

    /* Write snow data */
    stream->process(&it->snow.last_snow, 1, SNOW_LAST_SNOW);
    stream->process(&it->snow.MELTING, 1, SNOW_MELTING);
    stream->process(&it->snow.coverage, 1, SNOW_COVERAGE);
    stream->process(&it->snow.swq, 1, SNOW_SWQ);
    stream->process(&it->snow.surf_temp, 1, SNOW_SURF_TEMP);
    stream->process(&it->snow.surf_water, 1, SNOW_SURF_WATER);
    stream->process(&it->snow.pack_temp, 1, SNOW_PACK_TEMP);
    stream->process(&it->snow.pack_water, 1, SNOW_PACK_WATER);
    stream->process(&it->snow.density, 1, SNOW_DENSITY);
    stream->process(&it->snow.coldcontent, 1, SNOW_COLD_CONTENT);
    stream->process(&it->snow.snow_canopy, 1, SNOW_CANOPY);

    /* Write soil thermal node temperatures */
    stream->process(it->energy.T, state->options.Nnode, ENERGY_T);

    stream->processNewline();
  }

  if (state->options.LAKES) {
    for (int dist = 0; dist < Ndist; dist++) {
      // Store both wet and dry fractions if using distributed precipitation

      /* Write total soil moisture */
      double moistValues [state->options.Nlayer];
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        moistValues[lidx] = cell->prcp.lake_var.soil.layer[lidx].moist;   // Write specific.
      }
      stream->process(moistValues, state->options.Nlayer, LAKE_LAYER_MOIST);
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        cell->prcp.lake_var.soil.layer[lidx].moist = moistValues[lidx];   // Read specific.
      }

      /* Write average ice content */
#if SPATIAL_FROST
      double avgIceContent [state->options.Nlayer * FROST_SUBAREAS];
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        for (int frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
          avgIceContent[(lidx * FROST_SUBAREAS) + frost_area] = cell->prcp.lake_var.soil.layer[lidx].soil_ice[frost_area];  // Write specific.
        }
      }
      stream->process(avgIceContent, state->options.Nlayer * FROST_SUBAREAS, LAKE_LAYER_SOIL_ICE);
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        for (int frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
          cell->prcp.lake_var.soil.layer[lidx].soil_ice[frost_area] = avgIceContent[(lidx * FROST_SUBAREAS) + frost_area];  // Read specific.
        }
      }
#else
      double avgIceContent [state->options.Nlayer];
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        avgIceContent[lidx] = cell->prcp.lake_var.soil.layer[lidx].soil_ice;  // Write specific.
      }
      stream->process(avgIceContent, state->options.Nlayer, LAKE_LAYER_ICE_CONTENT);
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        cell->prcp.lake_var.soil.layer[lidx].soil_ice = avgIceContent[lidx];  // Read specific.
      }
#endif // SPATIAL_FROST
    }

    /* Write snow data */
    stream->process(&cell->prcp.lake_var.snow.last_snow,   1, LAKE_SNOW_LAST_SNOW);
    stream->process(&cell->prcp.lake_var.snow.MELTING,     1, LAKE_SNOW_MELTING);
    stream->process(&cell->prcp.lake_var.snow.coverage,    1, LAKE_SNOW_COVERAGE);
    stream->process(&cell->prcp.lake_var.snow.swq,         1, LAKE_SNOW_SWQ);
    stream->process(&cell->prcp.lake_var.snow.surf_temp,   1, LAKE_SNOW_SURF_TEMP);
    stream->process(&cell->prcp.lake_var.snow.surf_water,  1, LAKE_SNOW_SURF_WATER);
    stream->process(&cell->prcp.lake_var.snow.pack_temp,   1, LAKE_SNOW_PACK_TEMP);
    stream->process(&cell->prcp.lake_var.snow.pack_water,  1, LAKE_SNOW_PACK_WATER);
    stream->process(&cell->prcp.lake_var.snow.density,     1, LAKE_SNOW_DENSITY);
    stream->process(&cell->prcp.lake_var.snow.coldcontent, 1, LAKE_SNOW_COLD_CONTENT);
    stream->process(&cell->prcp.lake_var.snow.snow_canopy, 1, LAKE_SNOW_CANOPY);

    /* Write soil thermal node temperatures */
    stream->process(cell->prcp.lake_var.energy.T, state->options.Nnode, LAKE_ENERGY_T);

    /* Write lake-specific variables */
    stream->process(&cell->prcp.lake_var.activenod, 1, LAKE_ACTIVENOD);
    stream->process(&cell->prcp.lake_var.dz, 1, LAKE_DZ);
    stream->process(&cell->prcp.lake_var.surfdz, 1, LAKE_SURFDZ);
    stream->process(&cell->prcp.lake_var.ldepth, 1, LAKE_LDEPTH);
    stream->process(cell->prcp.lake_var.surface, cell->prcp.lake_var.activenod, LAKE_SURFACE);
    stream->process(&cell->prcp.lake_var.sarea, 1, LAKE_SAREA);
    stream->process(&cell->prcp.lake_var.volume, 1, LAKE_VOLUME);
    stream->process(cell->prcp.lake_var.temp, cell->prcp.lake_var.activenod, LAKE_TEMP);
    stream->process(&cell->prcp.lake_var.tempavg, 1, LAKE_TEMPAVG);
    stream->process(&cell->prcp.lake_var.areai, 1, LAKE_AREAI);
    stream->process(&cell->prcp.lake_var.new_ice_area, 1, LAKE_NEW_ICE_AREA);
    stream->process(&cell->prcp.lake_var.ice_water_eq, 1, LAKE_ICE_WATER_EQ);
    stream->process(&cell->prcp.lake_var.hice, 1, LAKE_HICE);
    stream->process(&cell->prcp.lake_var.tempi, 1, LAKE_TEMPI);
    stream->process(&cell->prcp.lake_var.swe, 1, LAKE_SWE);
    stream->process(&cell->prcp.lake_var.surf_temp, 1, LAKE_SURF_TEMP);
    stream->process(&cell->prcp.lake_var.pack_temp, 1, LAKE_PACK_TEMP);
    stream->process(&cell->prcp.lake_var.SAlbedo, 1, LAKE_SALBEDO);
    stream->process(&cell->prcp.lake_var.sdepth, 1, LAKE_SDEPTH);

    stream->processNewline();

  }
  /* Force file to be written */
  stream->flush();
}
