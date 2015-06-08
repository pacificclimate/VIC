#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include <string.h>

static char vcid[] = "$Id$";

int initialize_model_state(cell_info_struct* cell,
			   dmy_struct           dmy,
			   filep_struct         filep,
			   int                  Ndist,
			   const char*          initStateFilename,
			   const ProgramState  *state)
/**********************************************************************
  initialize_model_state      Keith Cherkauer	    April 17, 2000

  This routine initializes the model state (energy balance, water balance,
  and snow components).  If a state file is provided to the model than its
  contents are checked to see if it agrees with the current simulation
  set-up, if so it is used to initialize the model state.  If no state
  file is provided the model initializes all variables with defaults and
  the user should expect to throw out the beginning of the simulation 
  period as model start-up.

  UNITS: (m, s, kg, C, moisture in mm) unless otherwise specified

  Modifications:
  4-17-00 Modified from initialize_energy_bal.c and initialize_snow.c
          to provide a single controlling routine for initializing the
          model state.
  9-00    Fixed bug where initialization of soil node temperatures 
          and moistures was within two vegetation loops, thus only
          the first vegetation type was properly initialized.     KAC
  2-19-03 Modified to initialize soil and vegetation parameters for
          the dry grid cell fraction, if distributed precipitation
          is activated.                                           KAC
  11-18-02 Modified to initialize lake and wetland algorithms 
          variables.                                              LCB
  2-10-03 Fixed looping problem with initialization of soil moisture. KAC
  3-12-03 Modified so that soil layer ice content is only calculated 
          when frozen soil is implemented and active in the current 
          grid cell.                                                KAC
  04-10-03 Modified to read storm parameters from model state file.  KAC
  04-25-03 Modified to work with vegetation type specific storm 
           parameters.                                              KAC
  07-May-04 Initialize cell->soil_con.dz_node[Nnodes] to 0.0, since it is
	    accessed in set_node_parameters().				TJB
  01-Nov-04 Added support for state files containing SPATIAL_FROST
	    and LAKE_MODEL state variables.				TJB
  2006-Apr-21 Replaced Cv (uninitialized) with lake_con.Cl[0] in
	      surfstor calculation.					TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
              save_data structure to track changes in moisture storage
              over each time step; this needs initialization here.	TJB
  2006-Oct-10 Added snow[veg][band].snow_canopy to save_data.swe.	TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      This included removing the unused init_snow file.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2007-Apr-24 Added EXP_TRANS option.					JCA
  2007-Apr-24 Zsum_node loaded into soil_con structure for later use
	      without having to recalculate.				JCA
  2007-Aug-09 Added features for EXCESS_ICE option.			JCA
  2007-Aug-21 Return value of ErrorFlag if error in
	      distribute_node_moisture_properties.			JCA
  2007-Sep-18 Check for soil moist exceeding max moist moved from
	      read_initial_model_state to here.				JCA
  2007-Oct-24 Modified initialize_lake() to return ErrorFlag.		TJB
  2008-Mar-01 Reinserted missing logic for QUICK_FS in calls to
	      distribute_node_moisture_properties() and
	      estimate_layer_ice_content().				TJB
  2009-Feb-09 Removed dz_node from call to
	      distribute_node_moisture_properties.			KAC via TJB
  2009-Feb-09 Removed dz_node from call to find_0_degree_front.		KAC via TJB
  2009-Mar-15 Modified to not call estimate_layer_ice_content() if
	      not modeling frozen soil.					KAC via TJB
  2009-Mar-16 Added resid_moist to argument list of
	      estimate_layer_ice_content().  This allows computation
	      of min_liq, the minimum allowable liquid water content
	      in each layer as a function of temperature.		TJB
  2009-Jun-09 Modified to use extension of veg_lib structure to contain
	      bare soil information.					TJB
  2009-Jul-26 Added initial estimate of incoming longwave at surface
	      (LongUnderOut) for use in canopy snow T iteration.	TJB
  2009-Jul-31 Removed extra lake/wetland veg tile.			TJB
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.		TJB
  2009-Sep-19 Made initialization of Tfoliage more accurate for snow bands.	TJB
  2009-Sep-28 Added initialization of energy structure.			TJB
  2009-Nov-15 Added check to ensure that depth of first thermal node
	      is <= depth of first soil layer.				TJB
  2009-Dec-11 Removed initialization of save_data structure, since this
	      is now performed by the initial call to put_data().	TJB
  2009-Dec-11 Removed min_liq and state->options.MIN_LIQ.			TJB
  2010-Nov-11 Updated call to initialize_lake() to accommodate new
	      skip_hydro flag.						TJB
  2011-Mar-01 Updated calls to initialize_soil() and initialize_lake()
	      to accommodate new arguments.  Added more detailed validation
	      of soil moisture.						TJB
  2011-Mar-05 Added validation of initial soil moisture, ice, and snow
	      variables to make sure they are self-consistent.		TJB
  2011-May-31 Removed state->options.GRND_FLUX.  Now soil temperatures and
	      ground flux are always computed.				TJB
  2011-Jun-03 Added state->options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.				TJB
  2011-Jul-05 Changed logic initializing soil temperatures so that
	      type of initialization depends solely on state->options.QUICK_FLUX;
	      state->options.Nnodes is no longer automatically reset here.	TJB
**********************************************************************/
{
#if QUICK_FS
  extern const double temps[];
#endif

  char     ErrStr[MAXSTRING];
  char     FIRST_VEG;
  double   tmp_moist[MAX_LAYERS];
  double   tmp_runoff;
  int      frost_area;
  int      ErrorFlag;
  double   Zsum, dp;
  double   tmpdp, tmpadj, Bexp;
  double   Tair;
  double   tmp;
  double  *M;
  const int NUM_HRU = cell->prcp.hruList.size();
  double   moist[NUM_HRU][MAX_LAYERS];
#if SPATIAL_FROST
  double   ice[NUM_HRU][MAX_LAYERS][FROST_SUBAREAS];
#else
  double   ice[NUM_HRU][MAX_LAYERS];
#endif // SPATIAL_FROST
  double   Aufwc, Bufwc;
  double   Clake;
  double   precipitation_mu;
  double   surf_swq;
  double   pack_swq;
  double   sum_mindepth = 0, sum_depth_pre = 0, sum_depth_post = 0, tmp_mindepth = 0; //Excess Ice option variables

  double surf_temp = cell->atmos[0].air_temp[state->NR];

  // Initialize soil depths
  dp = cell->soil_con.dp;

  FIRST_VEG = TRUE;

  // increase initial soil surface temperature if air is very cold
  Tair = surf_temp;
  if ( surf_temp < -1. ) surf_temp = -1.;
  
  // initialize storm parameters to start a new simulation
  for (std::vector<HRU>::iterator hru = cell->prcp.hruList.begin(); hru != cell->prcp.hruList.end(); ++hru) {
    hru->init_STILL_STORM = 0;
    hru->init_DRY_TIME = INVALID_INT;
  }
  
  /********************************************
    Initialize all snow pack variables 
    - some may be reset if state file present
  ********************************************/

  initialize_snow(cell->prcp.hruList);

  /********************************************
    Initialize all soil layer variables 
    - some may be reset if state file present
  ********************************************/

  initialize_soil(cell->prcp.hruList, WET, &cell->soil_con, state);
  if ( state->options.DIST_PRCP )
    initialize_soil(cell->prcp.hruList, DRY, &cell->soil_con, state);

  /********************************************
    Initialize all vegetation variables 
    - some may be reset if state file present
  ********************************************/

  initialize_veg(cell->prcp.hruList, WET);
  if ( state->options.DIST_PRCP )
    initialize_veg(cell->prcp.hruList, DRY);

  /********************************************
    Initialize all lake variables 
  ********************************************/

  if ( state->options.LAKES && cell->lake_con.Cl[0] > 0) {
    for (std::vector<HRU>::iterator hru = cell->prcp.hruList.begin(); hru != cell->prcp.hruList.end(); ++hru) {
      if (hru->veg_con.vegClass == cell->lake_con.lake_idx) {
        hru->veg_con.LAKE = 1;
        ErrorFlag = initialize_lake(&cell->prcp.lake_var, cell->lake_con, &cell->soil_con, &(hru->cell[WET]), surf_temp, 0);
        if (ErrorFlag == ERROR) return(ErrorFlag);
      }
    }
  }

  /********************************************
    Initialize all spatial frost variables 
  ********************************************/

#if SPATIAL_FROST
  for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
    if ( FROST_SUBAREAS == 1 ) cell->soil_con.frost_fract[frost_area] = 1.;
    else if (FROST_SUBAREAS == 2 ) cell->soil_con.frost_fract[frost_area] = 0.5;
    else {
      cell->soil_con.frost_fract[frost_area] = 1. / (FROST_SUBAREAS - 1);
      if ( frost_area == 0 || frost_area == FROST_SUBAREAS-1 ) 
	cell->soil_con.frost_fract[frost_area] /= 2.;
    }
  }
#endif // SPATIAL_FROST

  /********************************************************
    Compute grid cell fractions for all subareas used in 
    spatial distribution of soil frost routines.
  ********************************************************/

#if QUICK_FS
#error // QUICK_FS is an untested code path. Continue at your own risk!
  if(state->options.FROZEN_SOIL) {

    /***********************************************************
      Prepare table of maximum unfrozen water content values
      - This linearizes the equation for maximum unfrozen water
        content, reducing computation time for the frozen soil
        model.
    ***********************************************************/

    for(lidx=0;lidx<state->options.Nlayer;lidx++) {
      for(ii=0;ii<QUICK_FS_TEMPS;ii++) {
	Aufwc = maximum_unfrozen_water(temps[ii], 1.0, 
				       cell->soil_con.bubble[lidx],
				       cell->soil_con.expt[lidx]);
	Bufwc = maximum_unfrozen_water(temps[ii+1], 1.0, 
				       cell->soil_con.bubble[lidx],
				       cell->soil_con.expt[lidx]);
	cell->soil_con.ufwc_table_layer[lidx][ii][0]
	  = linear_interp(0., temps[ii], temps[ii+1], Aufwc, Bufwc);
	cell->soil_con.ufwc_table_layer[lidx][ii][1]
	  = (Bufwc - Aufwc) / (temps[ii+1] - temps[ii]);
      }
    }
  }  
#endif // QUICK_FS



  // initialize miscellaneous energy balance terms
  for (std::vector<HRU>::iterator it = cell->prcp.hruList.begin(); it != cell->prcp.hruList.end(); ++it) {
	  /* Set fluxes to 0 */
	  it->energy.advected_sensible = 0.0;
	  it->energy.advection         = 0.0;
	  it->energy.AtmosError        = 0.0;
	  it->energy.AtmosLatent       = 0.0;
	  it->energy.AtmosLatentSub    = 0.0;
	  it->energy.AtmosSensible     = 0.0;
	  it->energy.canopy_advection  = 0.0;
	  it->energy.canopy_latent     = 0.0;
	  it->energy.canopy_latent_sub = 0.0;
	  it->energy.canopy_refreeze   = 0.0;
	  it->energy.canopy_sensible   = 0.0;
	  it->energy.deltaCC           = 0.0;
	  it->energy.deltaH            = 0.0;
	  it->energy.error             = 0.0;
	  it->energy.fusion            = 0.0;
	  it->energy.grnd_flux         = 0.0;
	  it->energy.latent            = 0.0;
	  it->energy.latent_sub        = 0.0;
	  it->energy.longwave          = 0.0;
	  it->energy.LongOverIn        = 0.0;
	  it->energy.LongUnderIn       = 0.0;
	  it->energy.LongUnderOut      = 0.0;
	  it->energy.melt_energy       = 0.0;
	  it->energy.NetLongAtmos      = 0.0;
	  it->energy.NetLongOver       = 0.0;
	  it->energy.NetLongUnder      = 0.0;
	  it->energy.NetShortAtmos     = 0.0;
	  it->energy.NetShortGrnd      = 0.0;
	  it->energy.NetShortOver      = 0.0;
	  it->energy.NetShortUnder     = 0.0;
	  it->energy.out_long_canopy   = 0.0;
	  it->energy.out_long_surface  = 0.0;
	  it->energy.refreeze_energy   = 0.0;
	  it->energy.sensible          = 0.0;
	  it->energy.shortwave         = 0.0;
	  it->energy.ShortOverIn       = 0.0;
	  it->energy.ShortUnderIn      = 0.0;
	  it->energy.snow_flux         = 0.0;
	  it->energy.AlbedoLake        = 0.0; //new
	  it->energy.AlbedoOver        = 0.0; //new
	  it->energy.fdepth[0]         = 0.0; //new
	  it->energy.fdepth[1]         = 0.0; //new
	  it->energy.fdepth[2]         = 0.0; //new
	  it->energy.tdepth[0]         = 0.0; //new
	  it->energy.tdepth[1]         = 0.0; //new
	  it->energy.tdepth[2]         = 0.0; //new
	  it->energy.unfrozen          = 0.0; //new
	  it->energy.glacier_flux      = 0.0; //new

	  /* Initial estimate of LongUnderOut for use by snow_intercept() */
	  tmp = it->energy.T[0] + KELVIN;
	  it->energy.LongUnderOut = STEFAN_B * tmp * tmp * tmp * tmp;
	  it->energy.Tfoliage     = Tair + cell->soil_con.Tfactor[it->bandIndex];

	  it->veg_var[0].canopyevap   = 0.0; //new
	  it->veg_var[1].canopyevap   = 0.0; //new

	  it->cell[0].layer[0].T      = 0.0; //new
	  it->cell[0].layer[1].T      = 0.0; //new
	  it->cell[0].layer[2].T      = 0.0; //new
	  it->cell[1].layer[0].T      = 0.0; //new
	  it->cell[1].layer[1].T      = 0.0; //new
	  it->cell[1].layer[2].T      = 0.0; //new

  }


  /************************************************************************
    CASE 1: Not using quick ground heat flux, and initial conditions files 
    provided
  ************************************************************************/

  if (state->options.INIT_STATE) {

#if EXCESS_ICE
#error // EXCESS_ICE is an untested code path. Continue at your own risk!
    sum_mindepth = 0;
    sum_depth_pre = 0;
    for( lidx = 0; lidx < state->options.Nlayer; lidx++ ) {
      tmp_mindepth = (float)(int)(cell->soil_con.min_depth[lidx] * 1000 + 0.5) / 1000;
      sum_mindepth += tmp_mindepth;
      sum_depth_pre += cell->soil_con.depth[lidx];
    }
#endif


    read_initial_model_state(initStateFilename, cell, cell->prcp.hruList.size(), Ndist, state);

#if EXCESS_ICE
    // calculate dynamic soil and veg properties if excess_ice is present
    sum_depth_post = 0;
    for( lidx = 0; lidx < state->options.Nlayer; lidx++ )
    sum_depth_post += cell->soil_con.depth[lidx];
    if( sum_depth_post != sum_depth_pre) {
      /*update soil_con properties*/
      for( lidx = 0; lidx < state->options.Nlayer; lidx++ ) {
        cell->soil_con.bulk_dens_min[lidx] *= (1.0-cell->soil_con.effective_porosity[lidx])*cell->soil_con.soil_density[lidx]/cell->soil_con.bulk_density[lidx];
        if (cell->soil_con.organic[lidx] > 0)
        cell->soil_con.bulk_dens_org[lidx] *= (1.0-cell->soil_con.effective_porosity[lidx])*cell->soil_con.soil_density[lidx]/cell->soil_con.bulk_density[lidx];
        cell->soil_con.bulk_density[lidx] = (1.0-cell->soil_con.effective_porosity[lidx])*cell->soil_con.soil_density[lidx];
        cell->soil_con.max_moist[lidx] = cell->soil_con.depth[lidx] * cell->soil_con.effective_porosity[lidx] * 1000.;
      } //loop for each soil layer      

      /********update remaining soil_con properties**********/
      /* update Maximum Infiltration for Upper Layers */
      if(state->options.Nlayer==2)
      cell->soil_con.max_infil = (1.0+cell->soil_con.b_infilt)*cell->soil_con.max_moist[0];
      else
      cell->soil_con.max_infil = (1.0+cell->soil_con.b_infilt)*(cell->soil_con.max_moist[0]+cell->soil_con.max_moist[1]);

      /* Soil Layer Critical and Wilting Point Moisture Contents */
      for(lidx=0;lidx<state->options.Nlayer;lidx++) { //soil layer
        cell->soil_con.Wcr[lidx] = cell->soil_con.Wcr_FRACT[lidx] * cell->soil_con.max_moist[lidx];
        cell->soil_con.Wpwp[lidx] = cell->soil_con.Wpwp_FRACT[lidx] * cell->soil_con.max_moist[lidx];
        if(cell->soil_con.Wpwp[lidx] > cell->soil_con.Wcr[lidx]) {
          sprintf(ErrStr,"Updated wilting point moisture (%f mm) is greater than updated critical point moisture (%f mm) for layer %d.\n\tIn the soil parameter file, Wpwp_FRACT MUST be <= Wcr_FRACT.\n",
              cell->soil_con.Wpwp[lidx], cell->soil_con.Wcr[lidx], lidx);
          nrerror(ErrStr);
        }
        if(cell->soil_con.Wpwp[lidx] < cell->soil_con.resid_moist[lidx] * cell->soil_con.depth[lidx] * 1000.) {
          sprintf(ErrStr,"Updated wilting point moisture (%f mm) is less than updated residual moisture (%f mm) for layer %d.\n\tIn the soil parameter file, Wpwp_FRACT MUST be >= resid_moist / (1.0 - bulk_density/soil_density).\n",
              cell->soil_con.Wpwp[lidx], cell->soil_con.resid_moist[lidx] * cell->soil_con.depth[lidx] * 1000., lidx);
          nrerror(ErrStr);
        }
      }

      /* If BASEFLOW = NIJSSEN2001 then convert ARNO baseflow
       parameters d1, d2, d3, and d4 to Ds, Dsmax, Ws, and c */
      if(state->options.BASEFLOW == NIJSSEN2001) {
        lidx = state->options.Nlayer-1;
        cell->soil_con.Dsmax = cell->soil_con.Dsmax_orig *
        pow((double)(1./(cell->soil_con.max_moist[lidx]-cell->soil_con.Ws_orig)), -cell->soil_con.c) +
        cell->soil_con.Ds_orig * cell->soil_con.max_moist[lidx];
        cell->soil_con.Ds = cell->soil_con.Ds_orig * cell->soil_con.Ws_orig / cell->soil_con.Dsmax_orig;
        cell->soil_con.Ws = cell->soil_con.Ws_orig/cell->soil_con.max_moist[lidx];
      }

      /*********** update root fractions ***************/
      calc_root_fractions(cell->prcp.hruList, soil_con);

#if VERBOSE
      /* write changes to screen */
      fprintf(stderr,"Soil properties initialized from state file:\n");
      for(lidx=0;lidx<state->options.Nlayer;lidx++) { //soil layer
        fprintf(stderr,"\tFor layer %d:\n",lidx+1);
        fprintf(stderr,"\t\tDepth of soil layer = %.2f m.\n",cell->soil_con.depth[lidx]);
        fprintf(stderr,"\t\tEffective porosity = %.2f.\n",cell->soil_con.effective_porosity[lidx]);
        fprintf(stderr,"\t\tBulk density = %.2f kg/m^3.\n",cell->soil_con.bulk_density[lidx]);
      }
      fprintf(stderr,"\tDamping depth = %.2f m.\n",cell->soil_con.dp);
      if(sum_depth_post == sum_mindepth)
      fprintf(stderr,"\tExcess ice is no longer present in the soil column.\n");
#endif //VERBOSE
    } //updated initial conditions due to state file
#endif //EXCESS_ICE

    /******Check that soil moisture does not exceed maximum allowed************/
    for (std::vector<HRU>::iterator it = cell->prcp.hruList.begin();
        it != cell->prcp.hruList.end(); ++it) {

      for (int dist = 0; dist < Ndist; dist++) {
        hru_data_struct& cellRef = it->cell[dist];
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {

          if (cellRef.layer[lidx].moist > cell->soil_con.max_moist[lidx]) {
            fprintf(stderr,
                "WARNING: Initial soil moisture (%f mm) exceeds maximum (%f mm) in layer %d for veg tile %d and snow band %d.  Resetting to maximum.\n",
                cellRef.layer[lidx].moist, cell->soil_con.max_moist[lidx], lidx,
                it->veg_con.vegIndex, it->bandIndex);
#if SPATIAL_FROST
            for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++)
            cellRef.layer[lidx].soil_ice[frost_area] *= cell->soil_con.max_moist[lidx]/cellRef.layer[lidx].moist;
#else
            cellRef.layer[lidx].soil_ice *= cell->soil_con.max_moist[lidx]
                / cellRef.layer[lidx].moist;
#endif
            cellRef.layer[lidx].moist = cell->soil_con.max_moist[lidx];
            tmp_moist[lidx] = cellRef.layer[lidx].moist;
          }
#if SPATIAL_FROST
          for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++) {
            if (cellRef.layer[lidx].soil_ice[frost_area] > cellRef.layer[lidx].moist)
            cellRef.layer[lidx].soil_ice[frost_area] = cellRef.layer[lidx].moist;
          }
#else
          if (cellRef.layer[lidx].soil_ice > cellRef.layer[lidx].moist)
            cellRef.layer[lidx].soil_ice = cellRef.layer[lidx].moist;
#endif
        }
        compute_runoff_and_asat(&cell->soil_con, tmp_moist, 0, &(cellRef.asat),
            &tmp_runoff, state);
      }

      // Override possible bad values of soil moisture under lake coming from state file
      // (ideally we wouldn't store these in the state file in the first place)
      if (state->options.LAKES && it->veg_con.vegClass == cell->lake_con.lake_idx) {
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
          cell->prcp.lake_var.soil.layer[lidx].moist =
              cell->soil_con.max_moist[lidx];
#if SPATIAL_FROST
          for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++) {
            if (cell->prcp.lake_var->soil.layer[lidx].soil_ice[frost_area] > cell->prcp.lake_var->soil.layer[lidx].moist)
            cell->prcp.lake_var->soil.layer[lidx].soil_ice[frost_area] = cell->prcp.lake_var->soil.layer[lidx].moist;
          }
#else
          if (cell->prcp.lake_var.soil.layer[lidx].soil_ice
              > cell->prcp.lake_var.soil.layer[lidx].moist)
            cell->prcp.lake_var.soil.layer[lidx].soil_ice =
                cell->prcp.lake_var.soil.layer[lidx].moist;
#endif
        }
      }
    }


    /****** initialize moist and ice ************/
    int hruIndex = 0;
    for (std::vector<HRU>::iterator it = cell->prcp.hruList.begin(); it != cell->prcp.hruList.end(); ++it, ++hruIndex) {

      // Initialize soil for existing vegetation types
      double Cv = it->veg_con.Cv;

      if (Cv > 0) {
        hru_data_struct& cellRef = it->cell[WET];
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
          moist[hruIndex][lidx] = cellRef.layer[lidx].moist;

#if SPATIAL_FROST
          for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
          ice[hruIndex][lidx][frost_area] = cellRef.layer[lidx].soil_ice[frost_area];
#else
          ice[hruIndex][lidx] = cellRef.layer[lidx].soil_ice;
#endif
        }
      }
    }

    /******Check that snow pack terms are self-consistent************/
    for (std::vector<HRU>::iterator it = cell->prcp.hruList.begin(); it != cell->prcp.hruList.end(); ++it) {
      if (it->snow.swq > MAX_SURFACE_SWE) {
        pack_swq = it->snow.swq - MAX_SURFACE_SWE;
        surf_swq = MAX_SURFACE_SWE;
      } else {
        pack_swq = 0;
        surf_swq = it->snow.swq;
        it->snow.pack_temp = 0;
      }
      if (it->snow.surf_water > LIQUID_WATER_CAPACITY * surf_swq) {
        it->snow.pack_water += it->snow.surf_water
            - (LIQUID_WATER_CAPACITY * surf_swq);
        it->snow.surf_water = LIQUID_WATER_CAPACITY * surf_swq;
      }
      if (it->snow.pack_water > LIQUID_WATER_CAPACITY * pack_swq) {
        it->snow.pack_water = LIQUID_WATER_CAPACITY * pack_swq;
      }
    }

  }
  
  /************************************************************************
    CASE 2: Initialize soil if using quick heat flux, and no initial
    soil properties file given
  ************************************************************************/
    
  else if(state->options.QUICK_FLUX) {

    /* Initialize soil node thicknesses */
    cell->soil_con.dz_node[0]   = cell->soil_con.depth[0];
    cell->soil_con.dz_node[1]   = cell->soil_con.depth[0];
    cell->soil_con.dz_node[2]   = 2. * (dp - 1.5 * cell->soil_con.depth[0]);
    cell->soil_con.Zsum_node[0] = 0;
    cell->soil_con.Zsum_node[1] = cell->soil_con.depth[0];
    cell->soil_con.Zsum_node[2] = dp;

    int hruIndex = 0;
    for (std::vector<HRU>::iterator it = cell->prcp.hruList.begin(); it != cell->prcp.hruList.end(); ++it, ++hruIndex) {
      // Initialize soil for existing vegetation types
      double Cv = it->veg_con.Cv;

      if (Cv > 0) {
        /* Initialize soil node temperatures */
        it->energy.T[0] = surf_temp;
        it->energy.T[1] = surf_temp;
//	    it->energy.T[2] = cell->soil_con.avg_temp + (surf_temp-cell->soil_con.avg_temp)*exp(-(cell->soil_con.Zsum_node[2]-cell->soil_con.Zsum_node[1])/dp);
        it->energy.T[2] = cell->soil_con.avg_temp;

        /* Initialize soil layer thicknesses */
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
          moist[hruIndex][lidx] = it->cell[WET].layer[lidx].moist;
#if SPATIAL_FROST
          for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
          ice[hruIndex][lidx][frost_area] = 0.;
#else
          ice[hruIndex][lidx] = 0.;
#endif
        }
      }
    }
  }

  /*****************************************************************
    CASE 3: Initialize Energy Balance Variables if not using quick
    ground heat flux, and no Initial Condition File Given 
  *****************************************************************/
  else if(!state->options.QUICK_FLUX) {
    int hruIndex = 0;
    for (std::vector<HRU>::iterator it = cell->prcp.hruList.begin(); it != cell->prcp.hruList.end(); ++it, ++hruIndex) {
      // Initialize soil for existing vegetation types
      double Cv = it->veg_con.Cv;

      if (Cv > 0) {

        if (!state->options.EXP_TRANS) {
          /* Initialize soil node temperatures and thicknesses
           Nodes set at surface, the depth of the first layer,
           twice the depth of the first layer, and at the
           damping depth.  Extra nodes are placed equal distance
           between the damping depth and twice the depth of the
           first layer. */

          it->energy.T[0] = surf_temp;
          cell->soil_con.dz_node[0] = cell->soil_con.depth[0];
          cell->soil_con.dz_node[1] = cell->soil_con.depth[0];
          cell->soil_con.dz_node[2] = cell->soil_con.depth[0];
          it->energy.T[state->options.Nnode - 1] = cell->soil_con.avg_temp;
          it->energy.T[1] = exp_interp(cell->soil_con.depth[0], 0., dp,
              surf_temp, cell->soil_con.avg_temp);
          it->energy.T[2] = exp_interp(2. * cell->soil_con.depth[0], 0.,
              dp, surf_temp, cell->soil_con.avg_temp);

          cell->soil_con.Zsum_node[0] = 0;
          cell->soil_con.Zsum_node[1] = cell->soil_con.depth[0];
          Zsum = 2. * cell->soil_con.depth[0];
          cell->soil_con.Zsum_node[2] = Zsum;
          tmpdp = dp - cell->soil_con.depth[0] * 2.5;
          tmpadj = 3.5;
          for (int index = 3; index < state->options.Nnode - 1; index++) {
            if (FIRST_VEG) {
              cell->soil_con.dz_node[index] = tmpdp / (((double) state->options.Nnode - tmpadj));
            }
            Zsum += (cell->soil_con.dz_node[index] + cell->soil_con.dz_node[index - 1])
                / 2.;
            cell->soil_con.Zsum_node[index] = Zsum;
            it->energy.T[index] = exp_interp(Zsum, 0.,
                cell->soil_con.dp, surf_temp, cell->soil_con.avg_temp);
          }
          if (FIRST_VEG) {
            FIRST_VEG = FALSE;
            cell->soil_con.dz_node[state->options.Nnode - 1] = (dp - Zsum
                - cell->soil_con.dz_node[state->options.Nnode - 2] / 2.) * 2.;
            Zsum += (cell->soil_con.dz_node[state->options.Nnode - 2]
                + cell->soil_con.dz_node[state->options.Nnode - 1]) / 2.;
            cell->soil_con.Zsum_node[state->options.Nnode - 1] = Zsum;
            if ((int) (Zsum * 1000 + 0.5) != (int) (dp * 1000 + 0.5)) {
              sprintf(ErrStr,
                  "Sum of thermal node thicknesses (%f) in initialize_model_state do not equal dp (%f), check initialization procedure",
                  Zsum, dp);
              nrerror(ErrStr);
            }
          }
        } else { /* exponential grid transformation, EXP_TRANS = TRUE*/

          /*calculate exponential function parameter */
          if (FIRST_VEG) {
            Bexp = log(dp + 1.) / (double) (state->options.Nnode - 1); //to force Zsum=dp at bottom node
            for (int index = 0; index <= state->options.Nnode - 1; index++)
              cell->soil_con.Zsum_node[index] = exp(Bexp * index) - 1.;
            if (cell->soil_con.Zsum_node[0] > cell->soil_con.depth[0]) {
              sprintf(ErrStr,
                  "Depth of first thermal node (%f) in initialize_model_state is greater than depth of first soil layer (%f); increase the number of nodes or decrease the thermal damping depth dp (%f)",
                  cell->soil_con.Zsum_node[0], cell->soil_con.depth[0], dp);
              nrerror(ErrStr);
            }
          }

          //top node
          int index = 0;
          if (FIRST_VEG)
            cell->soil_con.dz_node[index] = cell->soil_con.Zsum_node[index + 1]
                - cell->soil_con.Zsum_node[index];
          it->energy.T[index] = surf_temp;
          //middle nodes
          for (index = 1; index < state->options.Nnode - 1; index++) {
            if (FIRST_VEG) {
              cell->soil_con.dz_node[index] =
                  (cell->soil_con.Zsum_node[index + 1] - cell->soil_con.Zsum_node[index])
                      / 2.
                      + (cell->soil_con.Zsum_node[index]
                          - cell->soil_con.Zsum_node[index - 1]) / 2.;
            }
            it->energy.T[index] = exp_interp(
                cell->soil_con.Zsum_node[index], 0., cell->soil_con.dp, surf_temp,
                cell->soil_con.avg_temp);
          }
          //bottom node
          index = state->options.Nnode - 1;
          if (FIRST_VEG)
            cell->soil_con.dz_node[index] = cell->soil_con.Zsum_node[index]
                - cell->soil_con.Zsum_node[index - 1];
          it->energy.T[index] = cell->soil_con.avg_temp;
        }

        //initialize moisture and ice for each soil layer
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
          moist[hruIndex][lidx] = it->cell[WET].layer[lidx].moist;
#if SPATIAL_FROST
          for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
          ice[hruIndex][lidx][frost_area] = 0.;
#else
          ice[hruIndex][lidx] = 0.;
#endif
        }
      }
    }
  }

  /*********************************
    CASE 4: Unknown option
  *********************************/
  else {
    for (std::vector<HRU>::iterator hru = cell->prcp.hruList.begin(); hru != cell->prcp.hruList.end(); ++hru) {
      // Initialize soil for existing vegetation types

      if (hru->veg_con.Cv > 0) {
        for (int index = 0; index < state->options.Nlayer; index++) {
          cell->soil_con.dz_node[index] = 1.;
        }
      }
    }
  }

  /********************************************
    Initialize subsidence 
  ********************************************/

#if EXCESS_ICE
  for (int lidx = 0; lidx < state->options.Nlayer; lidx++ )
    cell->soil_con.subsidence[lidx] = 0.0;
    
#endif // EXCESS_ICE

  /******************************************
    Initialize soil thermal node properties 
  ******************************************/

  FIRST_VEG = TRUE;
  int hruIndex = 0;
  for (std::vector<HRU>::iterator it = cell->prcp.hruList.begin(); it != cell->prcp.hruList.end(); ++it, ++hruIndex) {
    // Initialize soil for existing vegetation types
    double Cv = it->veg_con.Cv;

    if (Cv > 0) {
      // Initialize soil for existing snow elevation bands
      if (cell->soil_con.AreaFract[it->bandIndex] > 0.) {

        /** Set soil properties for all soil nodes **/
        if (FIRST_VEG) {
          FIRST_VEG = FALSE;
          set_node_parameters(cell->soil_con.dz_node, cell->soil_con.Zsum_node,
              cell->soil_con.max_moist_node, cell->soil_con.expt_node,
              cell->soil_con.bubble_node, cell->soil_con.alpha,
              cell->soil_con.beta, cell->soil_con.gamma, cell->soil_con.depth,
              cell->soil_con.max_moist, cell->soil_con.expt,
              cell->soil_con.bubble, cell->soil_con.quartz,
              cell->soil_con.ufwc_table_node, cell->soil_con.porosity,
              cell->soil_con.effective_porosity, cell->soil_con.porosity_node,
              cell->soil_con.effective_porosity_node, state->options.Nnode,
              state->options.Nlayer, cell->soil_con.FS_ACTIVE, state);
        }

        /* set soil moisture properties for all soil thermal nodes */
        ErrorFlag = distribute_node_moisture_properties(&it->energy, &cell->soil_con, moist[hruIndex], state);
        if (ErrorFlag == ERROR)
          return (ErrorFlag);

        /* initialize layer moistures and ice contents */
        for (int curDist = 0; curDist < Ndist; curDist++) {
          hru_data_struct& cellRef = it->cell[curDist];
          for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
            cellRef.layer[lidx].moist = moist[hruIndex][lidx];
#if SPATIAL_FROST
            for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )

            cellRef.layer[lidx].soil_ice[frost_area] = ice[hruIndex][lidx][frost_area];
#else
            cellRef.layer[lidx].soil_ice = ice[hruIndex][lidx];
#endif
          }
          if (state->options.QUICK_FLUX) {
            ErrorFlag = estimate_layer_ice_content_quick_flux(cellRef.layer, it->energy.T[0], it->energy.T[1], &cell->soil_con, state);
          } else {
            ErrorFlag = estimate_layer_ice_content(cellRef.layer, it->energy.T, state->options.Nnode, state->options.Nlayer, &cell->soil_con, state);
          }
        }

        /* Find freezing and thawing front depths */
        if (!state->options.QUICK_FLUX && cell->soil_con.FS_ACTIVE)
          find_0_degree_fronts(&it->energy, cell->soil_con.Zsum_node, it->energy.T, state->options.Nnode);
      }
    }
    else if(Cv <= 0 || cell->soil_con.AreaFract[it->bandIndex] <= 0.) //new
    {
    	for(int nidx=0; nidx < state->options.Nnode; nidx++) //new
    	{
    		it->energy.moist[nidx] = 0; //new
    		it->energy.kappa_node[nidx] = 0; //new
    		it->energy.ice_content[nidx] = 0; //new
    		it->energy.Cs_node[nidx] = 0; //new
    	}

    }


  }	

  // initialize miscellaneous energy balance terms
//  for (std::vector<HRU>::iterator it = cell->prcp.hruList.begin(); it != cell->prcp.hruList.end(); ++it) {
//      /* Set fluxes to 0 */
//      it->energy.advected_sensible = 0.0;
//      it->energy.advection         = 0.0;
//      it->energy.AtmosError        = 0.0;
//      it->energy.AtmosLatent       = 0.0;
//      it->energy.AtmosLatentSub    = 0.0;
//      it->energy.AtmosSensible     = 0.0;
//      it->energy.canopy_advection  = 0.0;
//      it->energy.canopy_latent     = 0.0;
//      it->energy.canopy_latent_sub = 0.0;
//      it->energy.canopy_refreeze   = 0.0;
//      it->energy.canopy_sensible   = 0.0;
//      it->energy.deltaCC           = 0.0;
//      it->energy.deltaH            = 0.0;
//      it->energy.error             = 0.0;
//      it->energy.fusion            = 0.0;
//      it->energy.grnd_flux         = 0.0;
//      it->energy.latent            = 0.0;
//      it->energy.latent_sub        = 0.0;
//      it->energy.longwave          = 0.0;
//      it->energy.LongOverIn        = 0.0;
//      it->energy.LongUnderIn       = 0.0;
//      it->energy.LongUnderOut      = 0.0;
//      it->energy.melt_energy       = 0.0;
//      it->energy.NetLongAtmos      = 0.0;
//      it->energy.NetLongOver       = 0.0;
//      it->energy.NetLongUnder      = 0.0;
//      it->energy.NetShortAtmos     = 0.0;
//      it->energy.NetShortGrnd      = 0.0;
//      it->energy.NetShortOver      = 0.0;
//      it->energy.NetShortUnder     = 0.0;
//      it->energy.out_long_canopy   = 0.0;
//      it->energy.out_long_surface  = 0.0;
//      it->energy.refreeze_energy   = 0.0;
//      it->energy.sensible          = 0.0;
//      it->energy.shortwave         = 0.0;
//      it->energy.ShortOverIn       = 0.0;
//      it->energy.ShortUnderIn      = 0.0;
//      it->energy.snow_flux         = 0.0;
//      /* Initial estimate of LongUnderOut for use by snow_intercept() */
//      tmp = it->energy.T[0] + KELVIN;
//      it->energy.LongUnderOut = STEFAN_B * tmp * tmp * tmp * tmp;
//      it->energy.Tfoliage     = Tair + cell->soil_con.Tfactor[it->bandIndex];
//
//  }

  // initialize Tfallback counters
  for (std::vector<HRU>::iterator it = cell->prcp.hruList.begin(); it != cell->prcp.hruList.end(); ++it) {
    it->energy.Tfoliage_fbcount = 0;
    it->energy.Tcanopy_fbcount = 0;
    it->energy.Tsurf_fbcount = 0;
    for (int index = 0; index < state->options.Nnode - 1; index++) {
      it->energy.T_fbcount[index] = 0;
    }
  }

  return(0);
}








int update_thermal_nodes(dist_prcp_struct    *prcp,
			  int                  Nveg,
			  int                  Nnodes,
			  int                  Ndist,
			  soil_con_struct     *soil_con,
			  veg_con_struct      *veg_con,
			  const ProgramState  *state)
/**********************************************************************
  update_thermal_nodes           Jennifer Adam        August 16, 2007

  This routine is run after subsidence occurs (used only for EXCESS_ICE option).
  This routine updates the node depths and interpolates the current
  node temperatures to the new depths, then recalculates the nodal
  thermal properties.  Much of this routine is taken directly from
  initialize_model_state.

  Modifications:
  2009-Feb-09 Removed dz_node from call to
	      distribute_node_moisture_properties.			KAC via TJB
  2009-Feb-09 Removed dz_node from call to find_0_degree_front.		KAC via TJB
**********************************************************************/
{
  char     ErrStr[MAXSTRING];
  char     FIRST_VEG;
  int      ErrorFlag;
  double   Zsum;
  double   tmpdp, tmpadj, Bexp;

  double Tnode_prior[MAX_NODES];
  double Zsum_prior[MAX_NODES];
  

  FIRST_VEG = TRUE;

  /*****************************************************************
    Update soil thermal node depths, thicknesses, and temperatures.
    CASE 3: Initialize Energy Balance Variables if not using quick
    ground heat flux, and no Initial Condition File Given 
    (Currently this is the only case that works with EXCESS_ICE.)
  *****************************************************************/

  /*****************************************************************
    Update soil thermal node depths and thicknesses.
  *****************************************************************/
  //set previous Zsum
  for (int index = 0; index < Nnodes; index++ )
    Zsum_prior[index] = soil_con->Zsum_node[index];

  if(!state->options.EXP_TRANS){
    /* Nodes set at surface, the depth of the first layer,
       twice the depth of the first layer, and at the
       damping depth.  Extra nodes are placed equal distance
       between the damping depth and twice the depth of the
       first layer. */
    
    soil_con->dz_node[0] = soil_con->depth[0];
    soil_con->dz_node[1] = soil_con->depth[0];
    soil_con->dz_node[2] = soil_con->depth[0];	  
    soil_con->Zsum_node[0] = 0;
    soil_con->Zsum_node[1] = soil_con[0].depth[0];
    Zsum   = 2. * soil_con[0].depth[0];
    soil_con->Zsum_node[2] = Zsum;
    tmpdp  = soil_con->dp - soil_con[0].depth[0] * 2.5;
    tmpadj = 3.5;
    for (int index = 3; index < Nnodes-1; index++ ) {
      soil_con->dz_node[index] = tmpdp/(((double)Nnodes-tmpadj));
      Zsum += (soil_con->dz_node[index]
	       +soil_con->dz_node[index-1])/2.;
      soil_con->Zsum_node[index] = Zsum;
    }
    soil_con->dz_node[Nnodes-1] = (soil_con->dp - Zsum 
				   - soil_con->dz_node[Nnodes-2] 
				   / 2. ) * 2.;
    Zsum += (soil_con->dz_node[Nnodes-2]
	     +soil_con->dz_node[Nnodes-1])/2.;
    soil_con->Zsum_node[Nnodes-1] = Zsum;
    if((int)(Zsum*1000+0.5) != (int)(soil_con->dp*1000+0.5)) {
      sprintf(ErrStr,"Sum of thermal node thicknesses (%f) in initialize_model_state do not equal dp (%f), check initialization procedure",Zsum,soil_con->dp);
      nrerror(ErrStr);
    }
  }
  else{ /* exponential grid transformation, EXP_TRANS = TRUE*/
    
    /*calculate exponential function parameter */
    Bexp = log(soil_con->dp+1.)/(double)(Nnodes-1); //to force Zsum=dp at bottom node
    for (int index = 0; index <= Nnodes-1; index++ )
      soil_con->Zsum_node[index] = exp(Bexp*index)-1.;
    
    //top node	  
    int index=0;
    soil_con->dz_node[index] = soil_con->Zsum_node[index+1]-soil_con->Zsum_node[index];
    //middle nodes
    for ( index = 1; index < Nnodes-1; index++ ) {
      soil_con->dz_node[index] = (soil_con->Zsum_node[index+1]-soil_con->Zsum_node[index])/2.+(soil_con->Zsum_node[index]-soil_con->Zsum_node[index-1])/2.;
    }
    //bottom node
    index=Nnodes-1;
    soil_con->dz_node[index] = soil_con->Zsum_node[index]-soil_con->Zsum_node[index-1];
  }
#if VERBOSE
  fprintf(stderr,"More updated parameters in soil_con: dz_node and Zsum_node.\n");
#endif

  /******************************************
    Update soil thermal node temperatures via linear interpolation.
  ******************************************/
  for (std::vector<HRU>::iterator it = prcp->hruList.begin(); it != prcp->hruList.end(); ++it) {
    if (veg_con[it->veg_con.vegIndex].Cv > 0) {
      if (soil_con->AreaFract[it->bandIndex] > 0.) {
        //set previous temperatures
        for (int index = 0; index < Nnodes; index++)
          Tnode_prior[index] = it->energy.T[index];
        //top node: no need to update surface temperature
        //remaining nodes
        for (int index = 1; index < Nnodes; index++) {
          it->energy.T[index] = linear_interp(soil_con->Zsum_node[index],
              Zsum_prior[index - 1], Zsum_prior[index], Tnode_prior[index - 1],
              Tnode_prior[index]);
        }	  //node
      }
    }
  }

  /******************************************
    Update soil thermal node properties 
  ******************************************/  
  FIRST_VEG = TRUE;
  for (std::vector<HRU>::iterator it = prcp->hruList.begin(); it != prcp->hruList.end(); ++it) {
    if (veg_con[it->veg_con.vegIndex].Cv > 0) {
      // Initialize soil for existing snow elevation bands
      if (soil_con->AreaFract[it->bandIndex] > 0.) {
        /** Set soil properties for all soil nodes **/
        if (FIRST_VEG) {
          FIRST_VEG = FALSE;
          set_node_parameters(soil_con->dz_node, soil_con->Zsum_node,
              soil_con->max_moist_node, soil_con->expt_node,
              soil_con->bubble_node, soil_con->alpha, soil_con->beta,
              soil_con->gamma, soil_con->depth, soil_con->max_moist,
              soil_con->expt, soil_con->bubble, soil_con->quartz,
              soil_con->ufwc_table_node, soil_con->porosity,
              soil_con->effective_porosity, soil_con->porosity_node,
              soil_con->effective_porosity_node, Nnodes, state->options.Nlayer,
              soil_con->FS_ACTIVE, state);
        }

        double moist[MAX_LAYERS];
        for (int lidx = 0; lidx < state->options.Nlayer; lidx++)
          moist[lidx] = it->cell[WET].layer[lidx].moist;

        /* set soil moisture properties for all soil thermal nodes */
        if (!(state->options.LAKES && veg_con->LAKE != 0)) {
          ErrorFlag = distribute_node_moisture_properties(&it->energy, soil_con,
              moist, state);
          if (ErrorFlag == ERROR)
            return (ErrorFlag);
        }

        /* initialize layer moistures and ice contents */
        for (int curDist = 0; curDist < Ndist; curDist++) {
          hru_data_struct& cellRef = it->cell[curDist];
          if (!(state->options.LAKES && veg_con->LAKE != 0)) {
            if (state->options.QUICK_FLUX) {
              ErrorFlag = estimate_layer_ice_content_quick_flux(cellRef.layer,
                  it->energy.T[0], it->energy.T[1], soil_con, state);
            } else {
              ErrorFlag = estimate_layer_ice_content(cellRef.layer,
                  it->energy.T, Nnodes, state->options.Nlayer, soil_con, state);
            }
          }
        }

        /* Find freezing and thawing front depths */
        if (!state->options.QUICK_FLUX && soil_con->FS_ACTIVE)
          if (!(state->options.LAKES && veg_con->LAKE != 0))
            find_0_degree_fronts(&it->energy, soil_con->Zsum_node, it->energy.T,
                Nnodes);
      }
    }
  }

  return (0);
}
