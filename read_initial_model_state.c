#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void read_initial_model_state(FILE    *init_state,
			            dist_prcp_struct    *prcp,
			            int                  Nveg,
			            int                  Nbands,
			            int                  cellnum,
			            soil_con_struct     *soil_con,
			            int                  Ndist,
			            char                *init_STILL_STORM,
			            int                 *init_DRY_TIME,
			            lake_con_struct      lake_con,
			            const ProgramState  *state)
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
  
#if !NO_REWIND 
  rewind(init_state);
  
  /* skip header */
  if ( state->options.BINARY_STATE_FILE )
    fread(&tmpstr, sizeof(int)*5, 1, init_state);
  else {
    fgets(tmpstr, MAXSTRING, init_state);
    fgets(tmpstr, MAXSTRING, init_state);
  }
#else
#error // NO_REWIND is an untested code path. Continue at your own risk!
#endif
  
  /* read cell information */
  if ( state->options.BINARY_STATE_FILE ) {
    fread( &tmp_cellnum, sizeof(int), 1, init_state );
    fread( &tmp_Nveg, sizeof(int), 1, init_state );
    fread( &tmp_Nband, sizeof(int), 1, init_state );
    fread( &Nbytes, sizeof(int), 1, init_state );
  }
  else 
    fscanf( init_state, "%d %d %d", &tmp_cellnum, &tmp_Nveg, &tmp_Nband );
  // Skip over unused cell information
  while ( tmp_cellnum != cellnum && !feof(init_state) ) {
    if ( state->options.BINARY_STATE_FILE ) {
      // skip rest of current cells info
      for ( byte = 0; byte < Nbytes; byte++ ) 
	fread ( &tmpchar, 1, 1, init_state);
      // read info for next cell
      fread( &tmp_cellnum, sizeof(int), 1, init_state );
      fread( &tmp_Nveg, sizeof(int), 1, init_state );
      fread( &tmp_Nband, sizeof(int), 1, init_state );
      fread( &Nbytes, sizeof(int), 1, init_state );
    }
    else {
      // skip rest of current cells info
      fgets(tmpstr, MAXSTRING, init_state); // skip rest of general cell info
#if EXCESS_ICE      
      fgets(tmpstr, MAXSTRING, init_state); //excess ice info
#endif
      for (int veg = 0; veg <= tmp_Nveg; veg++ ) {
	fgets(tmpstr, MAXSTRING, init_state); // skip dist precip info
	for (int band = 0; band < tmp_Nband; band++ )
	  fgets(tmpstr, MAXSTRING, init_state); // skip snowband info
      }
      if ( state->options.LAKES ) {
        fgets(tmpstr, MAXSTRING, init_state); // skip lake info
      }
      // read info for next cell
      fscanf( init_state, "%d %d %d", &tmp_cellnum, &tmp_Nveg, &tmp_Nband );
    }//end if
  }//end while
  
  if ( feof(init_state) ) {
    sprintf(ErrStr, "Requested grid cell (%d) is not in the model state file.", 
	    cellnum);
    nrerror(ErrStr);
  }
  
  if ( tmp_Nveg != Nveg ) {
    sprintf(ErrStr,"The number of vegetation types in cell %d (%d) does not equal that defined in vegetation parameter file (%d).  Check your input files.", cellnum, tmp_Nveg, Nveg);
    nrerror(ErrStr);
  }
  if ( tmp_Nband != Nbands ) {
    sprintf(ErrStr,"The number of snow bands in cell %d (%d) does not equal that defined in the snow band file (%d).  Check your input files.", cellnum, tmp_Nband, Nbands);
    nrerror(ErrStr);
  }
 
  /* Read soil thermal node deltas */
  for (int nidx = 0; nidx < state->options.Nnode; nidx++ ) {
    if ( state->options.BINARY_STATE_FILE )
      fread( &soil_con->dz_node[nidx], sizeof(double), 1, init_state );
    else 
      fscanf( init_state, "%lf", &soil_con->dz_node[nidx] );
  }
  if ( state->options.Nnode == 1 ) soil_con->dz_node[0] = 0;
  
  /* Read soil thermal node depths */
  for (int nidx = 0; nidx < state->options.Nnode; nidx++ ) {
    if ( state->options.BINARY_STATE_FILE )
      fread( &soil_con->Zsum_node[nidx], sizeof(double), 1, init_state );
    else 
      fscanf( init_state, "%lf", &soil_con->Zsum_node[nidx] );
  }
  if ( state->options.Nnode == 1 ) soil_con->Zsum_node[0] = 0;
  if ( soil_con->Zsum_node[state->options.Nnode-1] - soil_con->dp > SMALL) {
    fprintf( stderr, "WARNING: Sum of soil nodes (%f) exceeds defined damping depth (%f).  Resetting damping depth.\n", soil_con->Zsum_node[state->options.Nnode-1], soil_con->dp );
    soil_con->dp = soil_con->Zsum_node[state->options.Nnode-1];
  }
  
  /* Read dynamic soil properties */
#if EXCESS_ICE
  /* Read soil depth */
  for (int lidx = 0; lidx < state->options.Nlayer; lidx++ ) {
    if ( state->options.BINARY_STATE_FILE )
      fread( &soil_con->depth[lidx], sizeof(double), 1, init_state );
    else 
      fscanf( init_state, "%lf", &soil_con->depth[lidx] );
  }
  
  /* Read effective porosity */
  for (int lidx = 0; lidx < state->options.Nlayer; lidx++ ) {
    if ( state->options.BINARY_STATE_FILE )
      fread( &soil_con->effective_porosity[lidx], sizeof(double), 1, init_state );
    else 
      fscanf( init_state, "%lf", &soil_con->effective_porosity[lidx] );
  }
  
  /* Reading damping depth */
  if ( state->options.BINARY_STATE_FILE )
    fread( &soil_con->dp, sizeof(double), 1, init_state );
  else 
    fscanf( init_state, "%lf", &soil_con->dp );
#endif //EXCESS_ICE
  
  /* Input for all vegetation types */
  for (int veg = 0; veg <= Nveg; veg++ ) {
    
    // read distributed precipitation variables
    if ( state->options.BINARY_STATE_FILE ) {
      fread( &prcp->mu[veg], sizeof(double), 1, init_state );
      fread( &init_STILL_STORM[veg], sizeof(char), 1, init_state );
      fread( &init_DRY_TIME[veg], sizeof(int), 1, init_state );
    }
    else {
      fscanf( init_state, "%lf %d %d", &prcp->mu[veg], &tmp_char, 
	      &init_DRY_TIME[veg] );
      init_STILL_STORM[veg] = (char)tmp_char;
    }
 
    /* Input for all snow bands */
    for (int band = 0; band < Nbands; band++ ) {
      /* Read cell identification information */
      int iveg = 0;
      int iband = 0;
      if ( state->options.BINARY_STATE_FILE ) {
	if ( fread( &iveg, sizeof(int), 1, init_state) != 1 ) 
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &iband, sizeof(int), 1, init_state) != 1 ) 
	  nrerror("End of model state file found unexpectedly");
      }
      else {
	if ( fscanf(init_state,"%d %d", &iveg, &iband) == EOF ) 
	  nrerror("End of model state file found unexpectedly");
      }
      if ( iveg != veg || iband != band ) {
	fprintf(stderr,"The vegetation and snow band indices in the model state file (veg = %d, band = %d) do not match those currently requested (veg = %d , band = %d).  Model state file must be stored with variables for all vegetation indexed by variables for all snow bands.\n", iveg, iband, veg, band);
	nrerror(ErrStr);
      }
      HRU* element = prcp->getHRUElement(veg, band);
      // Read both wet and dry fractions if using distributed precipitation
      for (int dist = 0; dist < Ndist; dist ++ ) {
        cell_data_struct& cellRef = element->cell[dist];
        /* Read total soil moisture */
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
          if (state->options.BINARY_STATE_FILE) {
            if (fread(&cellRef.layer[lidx].moist,
                sizeof(double), 1, init_state) != 1)
              nrerror("End of model state file found unexpectedly");
          } else {
            if (fscanf(init_state, " %lf",
                &cellRef.layer[lidx].moist) == EOF)
              nrerror("End of model state file found unexpectedly");
          }
        }

        /* Read average ice content */
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
#if SPATIAL_FROST
#error // SPATIAL_FROST is an untested code path. Continue at your own risk!
          for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
            if ( state->options.BINARY_STATE_FILE ) {
              if ( fread( &cellRef.layer[lidx].soil_ice[frost_area],
                      sizeof(double), 1, init_state ) != 1 )
              nrerror("End of model state file found unexpectedly");
            }
            else {
              if ( fscanf(init_state," %lf",
                      &cellRef.layer[lidx].soil_ice[frost_area]) == EOF )
              nrerror("End of model state file found unexpectedly");
            }
          }
#else
          if (state->options.BINARY_STATE_FILE) {
            if (fread(&cellRef.layer[lidx].soil_ice,
                sizeof(double), 1, init_state) != 1)
              nrerror("End of model state file found unexpectedly");
          } else {
            if (fscanf(init_state, " %lf",
                &cellRef.layer[lidx].soil_ice) == EOF)
              nrerror("End of model state file found unexpectedly");
          }
#endif // SPATIAL_FROST
        }

        /* Read dew storage */
        if (veg < Nveg) {
          if (state->options.BINARY_STATE_FILE) {
            if (fread(&element->veg_var[dist].Wdew, sizeof(double), 1,
                init_state) != 1)
              nrerror("End of model state file found unexpectedly");
          } else {
            if (fscanf(init_state, " %lf",
                &element->veg_var[dist].Wdew) == EOF)
              nrerror("End of model state file found unexpectedly");
          }
        }
      }
      
      /* Read snow data */
      if ( state->options.BINARY_STATE_FILE ) {
	if ( fread( &element->snow.last_snow, sizeof(int), 1,
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &element->snow.MELTING, sizeof(char), 1,
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &element->snow.coverage, sizeof(double), 1,
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &element->snow.swq, sizeof(double), 1,
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &element->snow.surf_temp, sizeof(double), 1,
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &element->snow.surf_water, sizeof(double), 1,
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &element->snow.pack_temp, sizeof(double), 1,
		    init_state ) != 1 ) 
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &element->snow.pack_water, sizeof(double), 1,
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &element->snow.density, sizeof(double), 1,
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &element->snow.coldcontent, sizeof(double), 1,
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &element->snow.snow_canopy, sizeof(double), 1,
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
      }
      else {
        if (fscanf(init_state, " %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &element->snow.last_snow, &tmp_char, &element->snow.coverage,
            &element->snow.swq, &element->snow.surf_temp,
            &element->snow.surf_water, &element->snow.pack_temp,
            &element->snow.pack_water, &element->snow.density,
            &element->snow.coldcontent, &element->snow.snow_canopy) == EOF)
          nrerror("End of model state file found unexpectedly");
        element->snow.MELTING = (char) tmp_char;
      }
      if (element->snow.density > 0.)
        element->snow.depth = 1000. * element->snow.swq / element->snow.density;
      
      /* Read soil thermal node temperatures */
      for (int nidx = 0; nidx < state->options.Nnode; nidx++) {
        if (state->options.BINARY_STATE_FILE) {
          if (fread(&element->energy.T[nidx], sizeof(double), 1,
              init_state) != 1)
            nrerror("End of model state file found unexpectedly");
        } else {
          if (fscanf(init_state, " %lf",
              &element->energy.T[nidx]) == EOF)
            nrerror("End of model state file found unexpectedly");
        }
      }
    }
  }
  if ( state->options.LAKES ) {
    if ( state->options.BINARY_STATE_FILE ) {
      // Read both wet and dry fractions if using distributed precipitation
      for (int dist = 0; dist < Ndist; dist ++ ) {
	
	/* Read total soil moisture */
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
          if (fread(&prcp->lake_var.soil.layer[lidx].moist, sizeof(double), 1,
              init_state) != 1)
            nrerror("End of model state file found unexpectedly");
        }

        /* Read average ice content */
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
#if SPATIAL_FROST
          for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
            if ( fread( &prcp->lake_var.soil.layer[lidx].soil_ice[frost_area], sizeof(double), 1, init_state ) != 1 )
            nrerror("End of model state file found unexpectedly");
          }
#else
          if (fread(&prcp->lake_var.soil.layer[lidx].soil_ice, sizeof(double),
              1, init_state) != 1)
            nrerror("End of model state file found unexpectedly");
#endif // SPATIAL_FROST
        }

      }
      
      /* Read snow data */
      if ( fread( &prcp->lake_var.snow.last_snow, sizeof(int), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.snow.MELTING, sizeof(char), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.snow.coverage, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.snow.swq, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.snow.surf_temp, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.snow.surf_water, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.snow.pack_temp, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.snow.pack_water, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.snow.density, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.snow.coldcontent, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.snow.snow_canopy, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if(prcp->lake_var.snow.density > 0.)
	prcp->lake_var.snow.depth = 1000. * prcp->lake_var.snow.swq / prcp->lake_var.snow.density;
      
      /* Read soil thermal node temperatures */
      for (int nidx = 0; nidx < state->options.Nnode; nidx++ ) {
	if ( fread( &prcp->lake_var.energy.T[nidx], sizeof(double), 1, init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
      }

      /* Read lake-specific variables */
      if ( fread( &prcp->lake_var.activenod, sizeof(int), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.dz, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.surfdz, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.ldepth, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      for ( node = 0; node <= prcp->lake_var.activenod; node++ ) {
        if ( fread( &prcp->lake_var.surface[node], sizeof(double), 1, init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
      }
      if ( fread( &prcp->lake_var.sarea, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.volume, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      for ( node = 0; node < prcp->lake_var.activenod; node++ ) {
        if ( fread( &prcp->lake_var.temp[node], sizeof(double), 1, init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
      }
      if ( fread( &prcp->lake_var.tempavg, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.areai, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.new_ice_area, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.ice_water_eq, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.hice, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.tempi, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.swe, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.surf_temp, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.pack_temp, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.coldcontent, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.surf_water, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.pack_water, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.SAlbedo, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &prcp->lake_var.sdepth, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      
    }
    else {
      // Read both wet and dry fractions if using distributed precipitation
      for (int dist = 0; dist < Ndist; dist ++ ) {
	
	/* Read total soil moisture */
	for (int lidx = 0; lidx < state->options.Nlayer; lidx++ ) {
	  if ( fscanf(init_state," %lf", &prcp->lake_var.soil.layer[lidx].moist) == EOF )
	    nrerror("End of model state file found unexpectedly");
	}
	
        /* Read average ice content */
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++ ) {
#if SPATIAL_FROST
	  for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
	    if ( fscanf(init_state," %lf", &prcp->lake_var.soil.layer[lidx].soil_ice[frost_area]) == EOF )
	      nrerror("End of model state file found unexpectedly");
	  }
#else
	  if ( fscanf(init_state," %lf", &prcp->lake_var.soil.layer[lidx].soil_ice) == EOF )
	    nrerror("End of model state file found unexpectedly");
#endif // SPATIAL_FROST
	}
	
      }
      
      /* Read snow data */
      if ( fscanf(init_state," %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		  &prcp->lake_var.snow.last_snow, &tmp_char,
		  &prcp->lake_var.snow.coverage, &prcp->lake_var.snow.swq,
		  &prcp->lake_var.snow.surf_temp, &prcp->lake_var.snow.surf_water,
		  &prcp->lake_var.snow.pack_temp, &prcp->lake_var.snow.pack_water,
		  &prcp->lake_var.snow.density, &prcp->lake_var.snow.coldcontent,
		  &prcp->lake_var.snow.snow_canopy)
	   == EOF ) 
        nrerror("End of model state file found unexpectedly");
      prcp->lake_var.snow.MELTING = (char)tmp_char;
      if(prcp->lake_var.snow.density > 0.)
	prcp->lake_var.snow.depth = 1000. * prcp->lake_var.snow.swq / prcp->lake_var.snow.density;
      
      /* Read soil thermal node temperatures */
      for (int nidx = 0; nidx < state->options.Nnode; nidx++ ) {
	if ( fscanf(init_state," %lf", &prcp->lake_var.energy.T[nidx]) == EOF )
	  nrerror("End of model state file found unexpectedly");
      }

      /* Read lake-specific variables */
      if ( fscanf(init_state," %d", &prcp->lake_var.activenod) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.dz) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.surfdz) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.ldepth) == EOF )
	nrerror("End of model state file found unexpectedly");
      for ( node = 0; node <= prcp->lake_var.activenod; node++ ) {
        if ( fscanf(init_state," %lf", &prcp->lake_var.surface[node]) == EOF )
	  nrerror("End of model state file found unexpectedly");
      }
      if ( fscanf(init_state," %lf", &prcp->lake_var.sarea) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.volume) == EOF )
	nrerror("End of model state file found unexpectedly");
      for ( node = 0; node < prcp->lake_var.activenod; node++ ) {
        if ( fscanf(init_state," %lf", &prcp->lake_var.temp[node]) == EOF )
	  nrerror("End of model state file found unexpectedly");
      }
      if ( fscanf(init_state," %lf", &prcp->lake_var.tempavg) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.areai) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.new_ice_area) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.ice_water_eq) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.hice) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.tempi) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.swe) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.surf_temp) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.pack_temp) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.coldcontent) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.surf_water) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.pack_water) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.SAlbedo) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &prcp->lake_var.sdepth) == EOF )
	nrerror("End of model state file found unexpectedly");
      
    }
  }

}
