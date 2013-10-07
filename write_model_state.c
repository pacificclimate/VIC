#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";


void write_model_state(dist_prcp_struct    *prcp,
		       int                  Nveg,
		       int                  cellnum,
		       filep_struct        *filep,
		       soil_con_struct     *soil_con,
		       char                *STILL_STORM,
		       int                 *DRY_TIME,
		       lake_con_struct      lake_con,
		       const ProgramState  *state)
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
  int    byte, Nbytes;
  int    frost_area;

  if(state->options.DIST_PRCP)
    Ndist = 2;
  else 
    Ndist = 1;
  Nbands = state->options.SNOW_BAND;

  /* write cell information */
  if ( state->options.BINARY_STATE_FILE ) {
    fwrite( &cellnum, sizeof(int), 1, filep->statefile );
    fwrite( &Nveg, sizeof(int), 1, filep->statefile );
    fwrite( &Nbands, sizeof(int), 1, filep->statefile );
  }
  else {
    fprintf( filep->statefile, "%i %i %i", cellnum, Nveg, Nbands );
  }
  // This stores the number of bytes from after this value to the end 
  // of the line.  DO NOT CHANGE unless you have changed the values
  // written to the state file.
  // IF YOU EDIT THIS FILE: UPDATE THIS VALUE!
  if ( state->options.BINARY_STATE_FILE ) {
    Nbytes = ( state->options.Nnode * sizeof(double) // dz_node
	       + state->options.Nnode * sizeof(double) // Zsum_node
#if EXCESS_ICE
	       + state->options.Nlayer * sizeof(double) //soil depth
	       + state->options.Nlayer * sizeof(double) //effective porosity
	       + sizeof(double) //damping depth
#endif
	       + (Nveg+1) * sizeof(double) // mu
	       + (Nveg+1) * sizeof(char) // STILL_STORM
	       + (Nveg+1) * sizeof(int) // DRY_TIME
	       + (Nveg+1) * Nbands * 2 * sizeof(int) // veg & band
	       + (Nveg+1) * Nbands * Ndist * state->options.Nlayer * sizeof(double) // soil moisture
#if SPATIAL_FROST
	       + (Nveg+1) * Nbands * Ndist * state->options.Nlayer * FROST_SUBAREAS * sizeof(double) // soil ice
#else
	       + (Nveg+1) * Nbands * Ndist * state->options.Nlayer * sizeof(double) // soil ice
#endif // SPATIAL_FROST
	       + Nveg * Nbands * Ndist * sizeof(double) // dew
	       + (Nveg+1) * Nbands * sizeof(int) // last_snow
	       + (Nveg+1) * Nbands * sizeof(char) // MELTING
	       + (Nveg+1) * Nbands * sizeof(double) * 9 // other snow parameters
	       + (Nveg+1) * Nbands * state->options.Nnode * sizeof(double) ); // soil temperatures
    if ( state->options.LAKES ) {
      /* Lake/wetland tiles have lake-specific state vars */
      Nbytes += sizeof(int) // activenod
	+ sizeof(double) // dz
	+ sizeof(double) // surfdz
	+ sizeof(double) // ldepth
	+ (prcp->lake_var.activenod+1) * sizeof(double) // surface
	+ sizeof(double) // sarea
	+ sizeof(double) // volume
	+ prcp->lake_var.activenod * sizeof(double) // temp
	+ sizeof(double) // tempavg
	+ sizeof(double) // areai
	+ sizeof(double) // new_ice_area
	+ sizeof(double) // ice_water_eq
	+ sizeof(double) // hice
	+ sizeof(double) // tempi
	+ sizeof(double) // swe
	+ sizeof(double) // surf_temp
	+ sizeof(double) // pack_temp
	+ sizeof(double) // coldcontent
	+ sizeof(double) // surf_water
	+ sizeof(double) // pack_water
	+ sizeof(double) // SAlbedo
	+ sizeof(double) // sdepth
	+ Ndist * state->options.Nlayer * sizeof(double) // soil moisture
#if SPATIAL_FROST
	+ Ndist * state->options.Nlayer * FROST_SUBAREAS * sizeof(double) // soil ice
#else
	+ Ndist * state->options.Nlayer * sizeof(double) // soil ice
#endif // SPATIAL_FROST
	+ sizeof(int) // last_snow
	+ sizeof(char) // MELTING
	+ sizeof(double) * 9 // other snow parameters
	+ state->options.Nnode * sizeof(double) // soil temperatures
	;
    }
    fwrite( &Nbytes, sizeof(int), 1, filep->statefile );
  }
  
  /* Write soil thermal node deltas */
  for (int nidx = 0; nidx < state->options.Nnode; nidx++ ) {
    if ( state->options.BINARY_STATE_FILE )
      fwrite( &soil_con->dz_node[nidx], sizeof(double), 1,
	      filep->statefile );
    else
      fprintf( filep->statefile, " %f ", soil_con->dz_node[nidx] );
  } 
  /* Write soil thermal node depths */
  for (int nidx = 0; nidx < state->options.Nnode; nidx++ ) {
    if ( state->options.BINARY_STATE_FILE )
      fwrite( &soil_con->Zsum_node[nidx], sizeof(double), 1, 
	      filep->statefile );
    else
      fprintf( filep->statefile, " %f ", soil_con->Zsum_node[nidx] );
  }    
  if ( !state->options.BINARY_STATE_FILE )
    fprintf( filep->statefile, "\n" );
  
  /* Write dynamic soil properties */
#if EXCESS_ICE
  /* Write soil depth */
  for (int lidx = 0; lidx < state->options.Nlayer; lidx++ ) {
    if ( state->options.BINARY_STATE_FILE )
      fwrite( &soil_con->depth[lidx], sizeof(double), 1,
	      filep->statefile );
    else
      fprintf( filep->statefile, "%f ", soil_con->depth[lidx] );
  }
  
  /* Write effective porosity */
  for (int lidx = 0; lidx < state->options.Nlayer; lidx++ ) {
    if ( state->options.BINARY_STATE_FILE )
      fwrite( &soil_con->effective_porosity[lidx], sizeof(double), 1,
	      filep->statefile );
    else
      fprintf( filep->statefile, "%f ", soil_con->effective_porosity[lidx] );
  }
  
  /* Write damping depth */
  if ( state->options.BINARY_STATE_FILE )
    fwrite( &soil_con->dp, sizeof(double), 1,
	    filep->statefile );
  else
    fprintf( filep->statefile, "%f\n", soil_con->dp );
#endif
  
  /* Output for all vegetation types */
  for (std::vector<HRU>::iterator it = prcp->hruList.begin(); it != prcp->hruList.end(); ++it) {
    
    // Do the following only once per vegetation type
    if (it->bandIndex == 0) {
      // Store distributed precipitation fraction
      if (state->options.BINARY_STATE_FILE)
        fwrite(&prcp->mu[it->vegIndex], sizeof(double), 1, filep->statefile);
      else
        fprintf(filep->statefile, "%f", prcp->mu[it->vegIndex]);

      // Store distributed precipitation variables
      if (state->options.BINARY_STATE_FILE) {
        fwrite(&STILL_STORM[it->vegIndex], sizeof(char), 1, filep->statefile);
        fwrite(&DRY_TIME[it->vegIndex], sizeof(int), 1, filep->statefile);
      } else {
        fprintf(filep->statefile, " %i %i\n", (int) STILL_STORM[it->vegIndex],
            DRY_TIME[it->vegIndex]);
      }
    }
    

    /* Output for all snow bands */
    /* Write cell identification information */
    if (state->options.BINARY_STATE_FILE) {
      fwrite(&it->vegIndex, sizeof(int), 1, filep->statefile);
      fwrite(&it->bandIndex, sizeof(int), 1, filep->statefile);
    } else {
      fprintf(filep->statefile, "%i %i", it->vegIndex, it->bandIndex);
    }

    for (int dist = 0; dist < Ndist; dist++) {
      cell_data_struct& cellRef = it->cell[dist];
      /* Write total soil moisture */
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
        double tmpval = cellRef.layer[lidx].moist; /* MPN */
        if (state->options.BINARY_STATE_FILE)
          fwrite(&tmpval, sizeof(double), 1, filep->statefile);
        else
          fprintf(filep->statefile, " %f", tmpval);
      }

      /* Write average ice content */
      for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
#if SPATIAL_FROST
#error
        for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
          tmpval = cellRef.layer[lidx].soil_ice[frost_area];
          if ( state->options.BINARY_STATE_FILE ) {
            fwrite( &tmpval, sizeof(double), 1, filep->statefile );
          }
          else {
            fprintf( filep->statefile, " %f", tmpval );
          }
        }
#else
        double tmpval = cellRef.layer[lidx].soil_ice;
        if (state->options.BINARY_STATE_FILE) {
          fwrite(&tmpval, sizeof(double), 1, filep->statefile);
        } else {
          fprintf(filep->statefile, " %f", tmpval);
        }
#endif // SPATIAL_FROST
      }

      /* Write dew storage */
      if (it->vegIndex < Nveg) {
        double tmpval = it->veg_var[dist].Wdew;
        if (state->options.BINARY_STATE_FILE)
          fwrite(&tmpval, sizeof(double), 1, filep->statefile);
        else
          fprintf(filep->statefile, " %f", tmpval);
      }
    }

    /* Write snow data */
    if (state->options.BINARY_STATE_FILE) {
      fwrite(&it->snow.last_snow, sizeof(int), 1, filep->statefile);
      fwrite(&it->snow.MELTING, sizeof(char), 1, filep->statefile);
      fwrite(&it->snow.coverage, sizeof(double), 1, filep->statefile);
      fwrite(&it->snow.swq, sizeof(double), 1, filep->statefile);
      fwrite(&it->snow.surf_temp, sizeof(double), 1, filep->statefile);
      fwrite(&it->snow.surf_water, sizeof(double), 1, filep->statefile);
      fwrite(&it->snow.pack_temp, sizeof(double), 1, filep->statefile);
      fwrite(&it->snow.pack_water, sizeof(double), 1, filep->statefile);
      fwrite(&it->snow.density, sizeof(double), 1, filep->statefile);
      fwrite(&it->snow.coldcontent, sizeof(double), 1, filep->statefile);
      fwrite(&it->snow.snow_canopy, sizeof(double), 1, filep->statefile);
    } else {
      fprintf(filep->statefile, " %i %i %f %f %f %f %f %f %f %f %f",
          it->snow.last_snow, (int) it->snow.MELTING,
          it->snow.coverage, it->snow.swq,
          it->snow.surf_temp, it->snow.surf_water,
          it->snow.pack_temp, it->snow.pack_water,
          it->snow.density, it->snow.coldcontent,
          it->snow.snow_canopy);
    }

    /* Write soil thermal node temperatures */
    for (int nidx = 0; nidx < state->options.Nnode; nidx++)
      if (state->options.BINARY_STATE_FILE)
        fwrite(&it->energy.T[nidx], sizeof(double), 1, filep->statefile);
      else
        fprintf(filep->statefile, " %f", it->energy.T[nidx]);

    if (!state->options.BINARY_STATE_FILE)
      fprintf(filep->statefile, "\n");
  }

  if (state->options.LAKES) {
    if (state->options.BINARY_STATE_FILE) {
      for (int dist = 0; dist < Ndist; dist++) {
        // Store both wet and dry fractions if using distributed precipitation

        /* Write total soil moisture */
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
          fwrite(&prcp->lake_var.soil.layer[lidx].moist, sizeof(double), 1,
              filep->statefile);
        }

        /* Write average ice content */
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
#if SPATIAL_FROST
          for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
            fwrite( &prcp->lake_var.soil.layer[lidx].soil_ice[frost_area], sizeof(double), 1, filep->statefile );
          }
#else
          fwrite(&prcp->lake_var.soil.layer[lidx].soil_ice, sizeof(double), 1,
              filep->statefile);
#endif // SPATIAL_FROST
        }
      }
      /* Write snow data */
      fwrite( &prcp->lake_var.snow.last_snow, sizeof(int), 1, filep->statefile );
      fwrite( &prcp->lake_var.snow.MELTING, sizeof(char), 1, filep->statefile );
      fwrite( &prcp->lake_var.snow.coverage, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.snow.swq, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.snow.surf_temp, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.snow.surf_water, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.snow.pack_temp, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.snow.pack_water, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.snow.density, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.snow.coldcontent, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.snow.snow_canopy, sizeof(double), 1, filep->statefile );
      
      /* Write soil thermal node temperatures */
      for (int nidx = 0; nidx < state->options.Nnode; nidx++ )
        fwrite( &prcp->lake_var.energy.T[nidx], sizeof(double), 1, filep->statefile );

      /* Write lake-specific variables */
      fwrite(&prcp->lake_var.activenod, sizeof(int), 1, filep->statefile);
      fwrite(&prcp->lake_var.dz, sizeof(double), 1, filep->statefile);
      fwrite(&prcp->lake_var.surfdz, sizeof(double), 1, filep->statefile);
      fwrite(&prcp->lake_var.ldepth, sizeof(double), 1, filep->statefile);
      for (int node = 0; node <= prcp->lake_var.activenod; node++) {
        fwrite(&prcp->lake_var.surface[node], sizeof(double), 1,
            filep->statefile);
      }
      fwrite( &prcp->lake_var.sarea, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.volume, sizeof(double), 1, filep->statefile );
      for (int node = 0; node < prcp->lake_var.activenod; node++ ) {
        fwrite( &prcp->lake_var.temp[node], sizeof(double), 1, filep->statefile );
      }
      fwrite( &prcp->lake_var.tempavg, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.areai, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.new_ice_area, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.ice_water_eq, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.hice, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.tempi, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.swe, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.surf_temp, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.pack_temp, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.coldcontent, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.surf_water, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.pack_water, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.SAlbedo, sizeof(double), 1, filep->statefile );
      fwrite( &prcp->lake_var.sdepth, sizeof(double), 1, filep->statefile );

    } else {

      for (int dist = 0; dist < Ndist; dist++) {
        // Store both wet and dry fractions if using distributed precipitation

        /* Write total soil moisture */
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
          fprintf(filep->statefile, " %f",
              prcp->lake_var.soil.layer[lidx].moist);
        }

        /* Write average ice content */
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
#if SPATIAL_FROST
          for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
            fprintf( filep->statefile, " %f", prcp->lake_var.soil.layer[lidx].soil_ice[frost_area] );
          }
#else
          fprintf(filep->statefile, " %f",
              prcp->lake_var.soil.layer[lidx].soil_ice);
#endif // SPATIAL_FROST
        }
      }

      /* Write snow data */
      fprintf( filep->statefile, " %i %i %f %f %f %f %f %f %f %f %f", 
		 prcp->lake_var.snow.last_snow, (int)prcp->lake_var.snow.MELTING,
		 prcp->lake_var.snow.coverage, prcp->lake_var.snow.swq,
		 prcp->lake_var.snow.surf_temp, prcp->lake_var.snow.surf_water,
		 prcp->lake_var.snow.pack_temp, prcp->lake_var.snow.pack_water,
		 prcp->lake_var.snow.density, prcp->lake_var.snow.coldcontent,
		 prcp->lake_var.snow.snow_canopy );
      
      /* Write soil thermal node temperatures */
      for (int nidx = 0; nidx < state->options.Nnode; nidx++ )
	fprintf( filep->statefile, " %f", prcp->lake_var.energy.T[nidx] );
      
      /* Write lake-specific variables */
      fprintf( filep->statefile, " %d", prcp->lake_var.activenod );
      fprintf( filep->statefile, " %f", prcp->lake_var.dz );
      fprintf( filep->statefile, " %f", prcp->lake_var.surfdz );
      fprintf( filep->statefile, " %f", prcp->lake_var.ldepth );
      for (int node = 0; node <= prcp->lake_var.activenod; node++ ) {
        fprintf( filep->statefile, " %f", prcp->lake_var.surface[node] );
      }
      fprintf( filep->statefile, " %f", prcp->lake_var.sarea );
      fprintf( filep->statefile, " %f", prcp->lake_var.volume );
      for (int node = 0; node < prcp->lake_var.activenod; node++ ) {
        fprintf( filep->statefile, " %f", prcp->lake_var.temp[node] );
      }
      fprintf( filep->statefile, " %f", prcp->lake_var.tempavg );
      fprintf( filep->statefile, " %f", prcp->lake_var.areai );
      fprintf( filep->statefile, " %f", prcp->lake_var.new_ice_area );
      fprintf( filep->statefile, " %f", prcp->lake_var.ice_water_eq );
      fprintf( filep->statefile, " %f", prcp->lake_var.hice );
      fprintf( filep->statefile, " %f", prcp->lake_var.tempi );
      fprintf( filep->statefile, " %f", prcp->lake_var.swe );
      fprintf( filep->statefile, " %f", prcp->lake_var.surf_temp );
      fprintf( filep->statefile, " %f", prcp->lake_var.pack_temp );
      fprintf( filep->statefile, " %f", prcp->lake_var.coldcontent );
      fprintf( filep->statefile, " %f", prcp->lake_var.surf_water );
      fprintf( filep->statefile, " %f", prcp->lake_var.pack_water );
      fprintf( filep->statefile, " %f", prcp->lake_var.SAlbedo );
      fprintf( filep->statefile, " %f", prcp->lake_var.sdepth );

      fprintf( filep->statefile, "\n" );
    }
  }
  /* Force file to be written */
  fflush(filep->statefile);

}

