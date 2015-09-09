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
  int Nbands = state->options.SNOW_BAND;
  int numHRUs = cell->prcp.hruList.size();
  //int numGMBterms = state->num_gmb_terms;

  StateIOContext context(filename, StateIO::Writer, state);
  StateIO* writer = context.stream;

  /* write cell information */
  writer->initializeDimensionIndices();
  writer->notifyDimensionUpdate(StateVariables::LAT_DIM, latitudeToIndex(cell->soil_con.lat, state));
  writer->notifyDimensionUpdate(StateVariables::LON_DIM, longitudeToIndex(cell->soil_con.lng, state));

  writer->write(&cell->soil_con.gridcel, 1, StateVariables::GRID_CELL);
  writer->write(&numHRUs, 1, StateVariables::VEG_TYPE_NUM);
  writer->write(&Nbands, 1, StateVariables::NUM_BANDS);
 // writer->write(&numGMBterms, 1, StateVariables::NUM_GLAC_MASS_BALANCE_INFO_TERMS);

  processCellForStateFile(cell, writer, state);
  
}

/*
 * The processCellForStateFile function is used for reading and writing state files (depending on the type of StateIO stream).
 * This method is also generic for each different state format type (binary, ascii, netCDF). This means that adding a variable
 * will make it written and read in the same place for each format. To add a new variable, use the following process as a guide:
 *
 * 1) Add it to the correct location below. For example if the variable is different per HRU then add it inside the HRU loop.
 * 2) Create a new entry in the StateMetaDataVariableIndices enum to uniquely describe this variable (see StateIO.h).
 * 3) Update the StateIONetCDF::populateMetaData() method to include a descriptive name for this variable and its dimensions.
 * 4) If a new dimension is needed then update the StateVariableDimension enum and add an entry in StateIONetCDF::populateMetaDimensions()
 *    (this step is probably not needed unless you added another loop below which does not iterate along an existing dimension).
 * 5) If you added another loop below, then make sure that stream->notifyDimensionUpdate is called appropriately for your dimension.
 * 6) Test. Test read and write for at least one format. The results of ASCII and NetCDF should be easily verifiable (human readable).
 */
void processCellForStateFile(cell_info_struct* cell, StateIO* stream, const ProgramState *state) {

  using namespace StateVariables;

  int Ndist;
  if(state->options.DIST_PRCP)
    Ndist = 2;
  else
    Ndist = 1;


////////////////

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

  stream->process(cell->soil_con.min_depth, state->options.Nlayer, SOIL_MIN_DEPTH); //new
  stream->process(cell->soil_con.porosity_node, state->options.Nnode, SOIL_POROSITY_NODE); //new
  stream->process(cell->soil_con.effective_porosity_node, state->options.Nnode, SOIL_EFFECTIVE_POROSITY_NODE); //new
  stream->process(cell->soil_con.subsidence, state->options.Nlayer, SOIL_SUBSIDENCE); //new
#endif


#if SPATIAL_SNOW
  stream->process(&cell->soil_con.depth_full_snow_cover, 1, SOIL_DEPTH_FULL_SNOW_COVER); //new
#endif

  /* Write Glacier Mass Balance equation for the grid cell */
  float gmbInfo[state->num_gmb_terms];
  // the following assumes that num_gmb_terms = 5
  gmbInfo[0] = cell->soil_con.gridcel;
  gmbInfo[1] = cell->gmbEquation.b0;
  gmbInfo[2] = cell->gmbEquation.b1;
  gmbInfo[3] = cell->gmbEquation.b2;
  gmbInfo[4] = cell->gmbEquation.fitError;
  stream->process(gmbInfo, state->num_gmb_terms, GLAC_MASS_BALANCE_INFO);

  /* Output for all vegetation types */
	int dimcount = 0;
  for (std::vector<HRU>::iterator it = cell->prcp.hruList.begin(); it != cell->prcp.hruList.end(); ++it) {
  	stream->notifyDimensionUpdate(HRU_DIM, dimcount);
  	dimcount++; // for the next iteration

    /* Output for all snow bands */
    /* Write cell identification information */
    int originalVeg = it->veg_con.vegClass;
    int originalBand = it->bandIndex;

    stream->process(&(it->veg_con.vegClass), 1, HRU_VEG_INDEX); //not necessary
    stream->process(&(it->bandIndex), 1, HRU_BAND_INDEX); //not necessary

    // The following is read specific and there is nothing we can do about it.
    if (stream->getType() == StateIO::Reader) {
      it->veg_con.vegIndex = getVegIndex(it->veg_con.vegClass, state);
      if (originalVeg != it->veg_con.vegClass || originalBand != it->bandIndex) {
        std::stringstream ss;
        ss << "The vegetation and snow band indices in the model state file (veg = " << originalVeg << ", band = " << originalBand << ")\n";
        ss << "do not match those currently requested (veg = " << it->veg_con.vegClass << ", band = " << it->bandIndex << "). \n";
        ss << "At HRU dimension index: " << stream->getCurrentDimensionIndex(HRU_DIM) << "\n";
        ss << "Model state file must be stored with variables for all vegetation indexed by variables for all snow bands.\n";
        throw VICException(ss.str());
      }
    }

    stream->process(&(it->mu), 1, PRCP_MU);
    stream->process(&(it->init_STILL_STORM), 1, INIT_STILL_STORM);
    stream->process(&(it->init_DRY_TIME), 1, INIT_DRY_TIME);

    stream->processNewline();

    for (int dist = 0; dist < Ndist; dist++) {
      stream->notifyDimensionUpdate(DIST_DIM, dist);


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
      double iceContent [state->options.Nlayer * FROST_SUBAREAS];
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        for (int frost = 0; frost < FROST_SUBAREAS; frost++) {
         iceContent[lidx * state->options.Nlayer + frost] = cellRef.layer[lidx].soil_ice[frost];  // Write specific.
        }
      }
      stream->process(iceContent, state->options.Nlayer * FROST_SUBAREAS, LAYER_SOIL_ICE);
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        for (int frost = 0; frost < FROST_SUBAREAS; frost++) {
         cellRef.layer[lidx].soil_ice[frost] = iceContent[lidx * state->options.Nlayer + frost];  // Read specific.
        }
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
      if (it->isArtificialBareSoil == false) {
        stream->process(&it->veg_var[dist].Wdew, 1, HRU_VEG_VAR_WDEW);
      }

    }

    /* Write snow data */

    stream->process(&it->snow.albedo, 1, SNOW_ALBEDO); //new
	  stream->process(&it->snow.canopy_albedo, 1, SNOW_CANOPY_ALBEDO); //new
	  stream->process(&it->snow.coldcontent, 1, SNOW_COLD_CONTENT);
    stream->process(&it->snow.coverage, 1, SNOW_COVERAGE);
    stream->process(&it->snow.density, 1, SNOW_DENSITY);
    stream->process(&it->snow.depth, 1, SNOW_DEPTH); //new
    stream->process(&it->snow.last_snow, 1, SNOW_LAST_SNOW);
    stream->process(&it->snow.MELTING, 1, SNOW_MELTING); //not necessary
    stream->process(&it->snow.pack_temp, 1, SNOW_PACK_TEMP);
    stream->process(&it->snow.pack_water, 1, SNOW_PACK_WATER);
    stream->process(&it->snow.snow_canopy, 1, SNOW_CANOPY);
    stream->process(&it->snow.surf_temp, 1, SNOW_SURF_TEMP);
    stream->process(&it->snow.surf_temp_fbcount, 1, SNOW_SURF_TEMP_FBCOUNT); //new
	  stream->process(&it->snow.surf_temp_fbflag, 1, SNOW_SURF_TEMP_FBFLAG); //new
    stream->process(&it->snow.surf_water, 1, SNOW_SURF_WATER);
    stream->process(&it->snow.swq, 1, SNOW_SWQ);
    stream->process(&it->snow.tmp_int_storage, 1, SNOW_TMP_INT_STORAGE); //new
    stream->process(&it->snow.surface_flux, 1, SNOW_SURFACE_FLUX); //new
	  stream->process(&it->snow.vapor_flux, 1, SNOW_VAPOR_FLUX); //new


    /* Write glacier data */
    stream->process(&it->glacier.surf_temp, 1, GLAC_SURF_TEMP); //new
    stream->process(&it->glacier.surf_temp_fbcount, 1, GLAC_SURF_TEMP_FBCOUNT); //new
    stream->process(&it->glacier.surf_temp_fbflag, 1, GLAC_SURF_TEMP_FBFLAG); //new
    stream->process(&it->glacier.Qnet, 1, GLAC_QNET); //new
    stream->process(&it->glacier.cum_mass_balance, 1, GLAC_CUM_MASS_BALANCE); //new
    stream->process(&it->glacier.vapor_flux, 1, GLAC_VAPOR_FLUX); //new
    stream->process(&it->glacier.water_storage, 1, GLAC_WATER_STORAGE); //new


    /* Write soil thermal node temperatures */
    stream->process(it->energy.T, state->options.Nnode, ENERGY_T);
    stream->process(it->energy.T_fbcount, state->options.Nnode, ENERGY_T_FBCOUNT); //new  -- is the error here?
    stream->process(&it->energy.Tcanopy_fbcount, 1, ENERGY_TCANOPY_FBCOUNT); //new
    stream->process(&it->energy.Tfoliage, 1, ENERGY_TFOLIAGE); //new
    stream->process(&it->energy.Tfoliage_fbcount, 1, ENERGY_TFOLIAGE_FBCOUNT); //new
    stream->process(&it->energy.Tsurf_fbcount, 1, ENERGY_TSURF_FBCOUNT); //new

    stream->processNewline();
  }

  if (state->options.LAKES) {
    for (int dist = 0; dist < Ndist; dist++) {
      stream->notifyDimensionUpdate(DIST_DIM, dist);
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
