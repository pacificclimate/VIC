#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

int  full_energy(char                 NEWCELL,
                 int                  time_step_record,
                 atmos_data_struct   *atmos,
                 dist_prcp_struct    *prcp,
                 const dmy_struct    *dmy,
                 lake_con_struct     *lake_con,
                 const soil_con_struct     *soil_con,
                 const veg_con_struct      *veg_con,
                 WriteDebug          *writeDebug,
                 const ProgramState  *state)
/**********************************************************************
	full_energy	Keith Cherkauer		January 8, 1997

  This subroutine controls the model core, it solves both the energy
  and water balance models, as well as frozen soils.  

  modifications:
  07-98 restructured to fix problems with distributed precipitation, 
        and to add the ability to solve the snow model at different 
	elevation bands within a single grid cell.                 KAC
  01-19-00 modified to work with the new atmosphere data structure 
           implemented when the radiation forcing routines were 
	   updated.  Also modified to use the new simplified
	   soil moisture storage for the frozen soil algorithm.    KAC
  12-01-00 modified to include the lakes and wetlands algorithm.   KAC
  11-18-02 Modified to handle blowing snow.  Also added debugging
           output for lake model.                                  LCB
  05-27-03 Updated passing of veg_con parameters for blowing snow
           to surface_fluxes.  Original did not account for the fact
           that veg_con is not allocated for veg = Nveg (bare soil)
           case.  This eliminates a memory error.                  KAC
  28-Sep-04 Added aero_resist_used to store the aerodynamic resistance
	    used in flux calculations.				TJB
  2006-Sep-23 Implemented flexible output configuration; now computation
	      of soil wetness and root zone soil moisture happens here.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.					TJB
  2007-Apr-04 Modified to handle grid cell errors by returning to the
	      main subroutine, rather than ending the simulation.		GCT/KAC
  2007-May-01 Added case of SPATIAL_FROST = TRUE in modifications
	      from 2006-Sep-23.							GCT
  2007-Aug-10 Added features for EXCESS_ICE option.
	      Including calculating subsidence for each layer and
	      updating soil depth, effective porosity,
	      bulk density, and soil moisture and fluxes by calling
	      runoff function if subsidence occurs.				JCA
  2007-Sep-07 No longer resets ice content to previous time-step ice content if
	      subsidence has occurred.						JCA
  2007-Sep-19 Added MAX_SUBSIDENCE parameter to EXCESS_ICE option.		JCA
  2007-Sep-19 Fixed bug in subsidence calculation.				JCA
  2007-Nov-06 Added veg_con to parameter list of lakemain().  Replaced
	      lake.fraci with lake.areai.					LCB via TJB
  2008-Jan-23 Changed ice0 from a scalar to an array.  Previously,
	      when state->options.SNOW_BAND > 1, the value of ice0 computed
	      for earlier bands was always overwritten by the value
	      of ice0 computed for the final band (even if the final
	      band had 0 area).							JS via TJB
  2008-May-05 Changed moist from a scalar to an array (moist0).  Previously,
	      when state->options.SNOW_BAND > 1, the value of moist computed
	      for earlier bands was always overwritten by the value
	      of moist computed for the final band (even if the final
	      band had 0 area).							KAC via TJB
  2009-Jan-16 Modified aero_resist_used and Ra_used to become arrays of
	      two elements (surface and overstory); added
	      state->options.AERO_RESIST_CANSNOW.					TJB
  2009-May-17 Added asat to cell_data.						TJB
  2009-Jun-09 Modified to use extension of state->veg_lib structure to contain
	      bare soil information.						TJB
  2009-Jun-09 Modified to compute aero_resist for all potential evap
	      landcover types.							TJB
  2009-Jun-09 Cell_data structure now only stores final aero_resist
	      values (called "aero_resist").  Preliminary uncorrected
	      aerodynamic resistances for current vegetation and various
	      reference land cover types for use in potential evap
	      calculations is stored in temporary array aero_resist.		TJB
  2009-Jun-26 Simplified argument list of runoff() by passing all cell_data
	      variables via a single reference to the cell data structure.	TJB
  2009-Jul-22 Fixed error in assignment of cell.aero_resist.			TJB
  2009-Jul-31 Wetland portion of lake/wetland tile is now processed in
	      full_energy() instead of wetland_energy().  Lake funcions are
	      now called directly from full_energy instead of lakemain().	TJB
  2009-Sep-28 Moved lake_snow and lake_energy into lake_var structure.
	      Removed calls to initialize_prcp and update_prcp.			TJB
  2009-Sep-30 Miscellaneous fixes for lake model.				TJB
  2009-Oct-05 Miscellaneous fixes for lake model, including updating/
	      rescaling of lake and wetland storages and fluxes to account
	      for changes in lake area.						TJB
  2009-Nov-09 Changed definition of lake->sarea to include ice extent; other
	      changes to handle case when lake fraction goes to 0.		LCB via TJB
  2010-Mar-31 Added runoff_in.							TJB
  2010-Sep-24 Changed atmos.runoff_in to atmos.channel_in.  Added
	      lake_var.channel_in to store it.					TJB
  2010-Nov-02 Changed units of lake_var moisture fluxes to volume (m3).		TJB
  2010-Nov-26 Changed argument list of water_balance().				TJB
  2011-May-31 Prepare_full_energy() is now always called.			TJB
  2011-Jun-03 Added state->options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.					TJB
  2012-Jan-01 Modified condition for determining whether to simulate lakes
	      to check whether lake_idx >= 0.					TJB

**********************************************************************/
{
  char                   overstory;
  int                    Ndist;
  int                    Nveg;
  int                    veg_class;
  int                    Nbands;
  int                    ErrorFlag;
  int                    frost_area;
  double                 out_prec[2*MAX_BANDS];
  double                 out_rain[2*MAX_BANDS];
  double                 out_snow[2*MAX_BANDS];
  double                 out_short=0;
  double                 dp;
  double                 ice0[MAX_BANDS];
  double                 moist0[MAX_BANDS];
  double                 surf_atten;
  double                 Tend_surf;
  double                 Tend_grnd;
  double                 wind_h;
  double                 height;
  double                 displacement[3];
  double                 roughness[3];
  double                 ref_height[3];
  double               **aero_resist;
  double                 Cv;
  double                 latent_heat_Le;
  double                 Melt[2*MAX_BANDS];
  double                 bare_albedo;
  double                 snow_inflow[MAX_BANDS];
  double                 rainonly;
  double                 sum_runoff;
  double                 sum_baseflow;
  double                 tmp_wind[3];
  double                 tmp_mu;
  double                 tmp_total_moist;
  double                 gauge_correction[2];
  float 	               lag_one;
  float 	               sigma_slope;
  float  	               fetch;
  int                    pet_veg_class;
  double                 lakefrac;
  double                 fraci;
  double                 wetland_runoff;
  double                 wetland_baseflow;
  double                 oldsnow;
  double                 snowprec;
  double                 rainprec;

  //Excess ice variables
  int                    SubsidenceUpdate = 0;
  int                    index = 0;
  char                   ErrStr[MAXSTRING];
  double                 max_ice_layer = 0; //mm/mm
  double                 ave_ice_fract = 0; //mm/mm
  double                 ave_ice = 0, tmp_ice = 0; //mm
  double                 ice_layer = 0; //mm
  double                 subsidence[MAX_LAYERS]; //mm
  double                 total_subsidence = 0; //m
  double                 tmp_subsidence = 0; //mm
  double                 total_meltwater = 0; //mm
  double                 tmp_depth = 0, tmp_depth_prior = 0; //m
  double                 ppt[2]; 
  double                 moist_prior[2][MAX_VEG][MAX_BANDS][MAX_LAYERS]; //mm
  double                 evap_prior[2][MAX_VEG][MAX_BANDS][MAX_LAYERS]; //mm

  /* Allocate aero_resist array */
  aero_resist = (double**) calloc(N_PET_TYPES + 1, sizeof(double*));
  for (int p = 0; p < N_PET_TYPES + 1; p++) {
    aero_resist[p] = (double*) calloc(3, sizeof(double));
  }

  /* set variables for distributed precipitation */
  if (state->options.DIST_PRCP)
    Ndist = 2;
  else
    Ndist = 1;
  Nbands = state->options.SNOW_BAND;

  /* Set number of vegetation types */
  Nveg = veg_con[0].vegetat_type_num;

  /** Set Damping Depth **/
  dp = soil_con->dp;

  /* Compute gauge undercatch correction factors 
   - this assumes that the gauge is free of vegetation effects, so gauge
   correction is constant for the entire grid cell */
  if (state->options.CORRPREC && atmos->prec[state->NR] > 0)
    correct_precip(gauge_correction, atmos->wind[state->NR], state->global_param.wind_h,
        soil_con->rough, soil_con->snow_rough);
  else {
    gauge_correction[0] = 1;
    gauge_correction[1] = 1;
  }
  atmos->out_prec = 0;
  atmos->out_rain = 0;
  atmos->out_snow = 0;

  /* initialize prior moist and ice for subsidence calculations */
#if EXCESS_ICE
  for (std::vector<HRU>::iterator it = prcp->hruList.begin(); it != prcp->hruList.end(); ++it) {
    for ( dist = 0; dist < Ndist; dist++ ) {
      for(lidx=0;lidx<state->options.Nlayer;lidx++) {
        moist_prior[dist][it->vegIndex][it->bandIndex][lidx] = it->cell[dist].layer[lidx].moist;
        evap_prior[dist][it->vegIndex][it->bandIndex][lidx] = 0; //initialize
      }
    }
  }
#endif //EXCESS_ICE

  /**************************************************
   Solve Energy and/or Water Balance for Each
   Vegetation Type
   **************************************************/
  for (std::vector<HRU>::iterator it = prcp->hruList.begin(); it != prcp->hruList.end(); ++it) {

    /** Solve Veg Type only if Coverage Greater than 0% **/
    if (veg_con[it->vegIndex].Cv > 0.0) {
      Cv = veg_con[it->vegIndex].Cv;
      Nbands = state->options.SNOW_BAND;

      /** Lake-specific processing **/
      if (veg_con[it->vegIndex].LAKE) {

        /* Update areai to equal new ice area from previous time step. */
        prcp->lake_var.areai = prcp->lake_var.new_ice_area;

        /* Compute lake fraction and ice-covered fraction */
        if (prcp->lake_var.areai < 0)
          prcp->lake_var.areai = 0;
        if (prcp->lake_var.sarea > 0) {
          fraci = prcp->lake_var.areai / prcp->lake_var.sarea;
          if (fraci > 1.0)
            fraci = 1.0;
        } else
          fraci = 0.0;
        lakefrac = prcp->lake_var.sarea / lake_con->basin[0];

        Nbands = 1;
        Cv *= (1 - lakefrac);

        if (Cv == 0)
          continue;

      }

      /**************************************************
       Initialize Model Parameters
       **************************************************/

      if (soil_con->AreaFract[it->bandIndex] > 0) {

        /* Initialize prcp->energy balance variables */
        it->energy.shortwave = 0;
        it->energy.longwave = 0.;

        /* Initialize snow variables */
        it->snow.vapor_flux = 0.;
        it->snow.canopy_vapor_flux = 0.;
        snow_inflow[it->bandIndex] = 0.;
        Melt[it->bandIndex * 2] = 0.;

      }

      if (it->bandIndex == 0) {
        /* Initialize precipitation storage */
        for (int j = 0; j < 2 * MAX_BANDS; j++) {
          out_prec[j] = 0;
          out_rain[j] = 0;
          out_snow[j] = 0;
        }

        /** Define vegetation class number **/
        veg_class = veg_con[it->vegIndex].veg_class;

        /** Assign wind_h **/
        /** Note: this is ignored below **/
        wind_h = state->veg_lib[veg_class].wind_h;

        /** Compute Surface Attenuation due to Vegetation Coverage **/
        surf_atten =
            exp(
                -state->veg_lib[veg_class].rad_atten
                    * state->veg_lib[veg_class].LAI[dmy[time_step_record].month
                        - 1]);

        /* Initialize soil thermal properties for the top two layers */
        prepare_full_energy(it->vegIndex, Nveg, state->options.Nnode, prcp, soil_con,
            moist0, ice0, state);

        /** Compute Bare (free of snow) Albedo **/
        bare_albedo =
            state->veg_lib[veg_class].albedo[dmy[time_step_record].month - 1];

        /*************************************
         Compute the aerodynamic resistance
         for current veg cover and various
         types of potential evap
         *************************************/

        /* Loop over types of potential evap, plus current veg */
        /* Current veg will be last */
        for (int p = 0; p < N_PET_TYPES + 1; p++) {

          /* Initialize wind speeds */
          tmp_wind[0] = atmos->wind[state->NR];
          tmp_wind[1] = INVALID;
          tmp_wind[2] = INVALID;

          /* Set surface descriptive variables */
          if (p < N_PET_TYPES_NON_NAT) {
            pet_veg_class = state->veg_lib[0].NVegLibTypes + p;
          } else {
            pet_veg_class = veg_class;
          }
          displacement[0] =
              state->veg_lib[pet_veg_class].displacement[dmy[time_step_record].month
                  - 1];
          roughness[0] =
              state->veg_lib[pet_veg_class].roughness[dmy[time_step_record].month
                  - 1];
          overstory = state->veg_lib[pet_veg_class].overstory;
          if (p >= N_PET_TYPES_NON_NAT)
            if (roughness[0] == 0)
              roughness[0] = soil_con->rough;

          /* Estimate vegetation height */
          height = calc_veg_height(displacement[0]);

          /* Estimate reference height */
          if (displacement[0] < wind_h)
            ref_height[0] = wind_h;
          else
            ref_height[0] = displacement[0] + wind_h + roughness[0];

          /* Compute aerodynamic resistance over various surface types */
          /* Do this not only for current veg but also all types of PET */
          int ErrorFlag = CalcAerodynamic(overstory, height,
              state->veg_lib[pet_veg_class].trunk_ratio, soil_con->snow_rough,
              soil_con->rough, state->veg_lib[pet_veg_class].wind_atten,
              aero_resist[p], tmp_wind, displacement, ref_height, roughness);
          if (ErrorFlag == ERROR)
            return (ERROR);

        }

      }

      /* Initialize final aerodynamic resistance values */
      if (soil_con->AreaFract[it->bandIndex] > 0) {
        it->cell[WET].aero_resist[0] = aero_resist[N_PET_TYPES][0];
        it->cell[WET].aero_resist[1] = aero_resist[N_PET_TYPES][1];
      }


      /**************************************************
       Store Water Balance Terms for Debugging
       **************************************************/

#if LINK_DEBUG
      if(state->debug.DEBUG || state->debug.PRT_MOIST || state->debug.PRT_BALANCE) {
        /** Compute current total moisture for water balance check **/
        if (it->bandIndex == 0) {
          store_moisture_for_debug(it->vegIndex, Nveg, prcp->mu, prcp->hruList, soil_con, state);
        }
        if (state->debug.PRT_BALANCE) {
          for (int j = 0; j < Ndist; j++) {
            if (soil_con->AreaFract[it->bandIndex] > 0) {
              for (int i = 0; i < state->options.Nlayer + 3; i++) {
                state->debug.inflow[j][it->bandIndex][i] = 0;
                state->debug.outflow[j][it->bandIndex][i] = 0;
              }

            }
          }
        }
      }
#endif // LINK_DEBUG
      /******************************
       Solve ground surface fluxes
       ******************************/

      if (soil_con->AreaFract[it->bandIndex] > 0) {

        lag_one = veg_con[it->vegIndex].lag_one;
        sigma_slope = veg_con[it->vegIndex].sigma_slope;
        fetch = veg_con[it->vegIndex].fetch;

        /* Initialize pot_evap */
        for (int p = 0; p < N_PET_TYPES; p++)
          it->cell[WET].pot_evap[p] = 0;

        ErrorFlag = surface_fluxes(overstory, bare_albedo, height, ice0[it->bandIndex],
            moist0[it->bandIndex], SubsidenceUpdate, evap_prior[DRY][it->vegIndex][it->bandIndex],
            evap_prior[WET][it->vegIndex][it->bandIndex], prcp->mu[it->vegIndex], surf_atten,
            &(Melt[it->bandIndex * 2]), &latent_heat_Le, aero_resist, displacement,
            gauge_correction, &out_prec[it->bandIndex * 2], &out_rain[it->bandIndex * 2],
            &out_snow[it->bandIndex * 2], ref_height, roughness, &snow_inflow[it->bandIndex],
            tmp_wind, veg_con[it->vegIndex].root, Nbands, Ndist, state->options.Nlayer,
            Nveg, it->bandIndex, dp, it->vegIndex, time_step_record, veg_class, atmos, dmy,
            &(it->energy), &(it->cell[DRY]),
            &(it->cell[WET]), &(it->snow), soil_con,
            &(it->veg_var[DRY]),
            &(it->veg_var[WET]), lag_one, sigma_slope, fetch,
            state);

        if (ErrorFlag == ERROR)
          return (ERROR);

        atmos->out_prec += out_prec[it->bandIndex * 2] * Cv * soil_con->AreaFract[it->bandIndex];
        atmos->out_rain += out_rain[it->bandIndex * 2] * Cv * soil_con->AreaFract[it->bandIndex];
        atmos->out_snow += out_snow[it->bandIndex * 2] * Cv * soil_con->AreaFract[it->bandIndex];

        /********************************************************
         Compute soil wetness and root zone soil moisture
         ********************************************************/
        // Loop through distributed precipitation fractions
        for (int dist = 0; dist < Ndist; dist++) {
          cell_data_struct& cellRef = it->cell[dist];
          cellRef.rootmoist = 0;
          cellRef.wetness = 0;
          for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
            if (veg_con->root[lidx] > 0) {
              cellRef.rootmoist += cellRef.layer[lidx].moist;
            }
#if EXCESS_ICE
            cellRef.wetness += (cellRef.layer[lidx].moist - soil_con->Wpwp[lidx])/(soil_con->effective_porosity[lidx]*soil_con->depth[lidx]*1000 - soil_con->Wpwp[lidx]);
#else
            cellRef.wetness +=
                (cellRef.layer[lidx].moist - soil_con->Wpwp[lidx])
                    / (soil_con->porosity[lidx] * soil_con->depth[lidx] * 1000
                        - soil_con->Wpwp[lidx]);
#endif
          }
          cellRef.wetness /= state->options.Nlayer;
        }
      }

      /****************************
       Controls Debugging Output
       ****************************/
#if LINK_DEBUG

      for(int j = 0; j < Ndist; j++) {
        if(j == 0)
        tmp_mu = prcp->mu[it->vegIndex];
        else
        tmp_mu = 1. - prcp->mu[it->vegIndex];
        /** for debugging water balance: [0] = vegetation,
         [1] = ground prcp->snow, [2..Nlayer+1] = soil layers **/
        if (state->debug.PRT_BALANCE) {
          if (soil_con->AreaFract[it->bandIndex] > 0) {
            state->debug.inflow[j][it->bandIndex][state->options.Nlayer + 2] +=
                out_prec[j + it->bandIndex * 2] * soil_con->Pfactor[it->bandIndex];
            state->debug.inflow[j][it->bandIndex][0] = 0.;
            state->debug.inflow[j][it->bandIndex][1] = 0.;
            state->debug.outflow[j][it->bandIndex][0] = 0.;
            state->debug.outflow[j][it->bandIndex][1] = 0.;
            state->debug.inflow[j][it->bandIndex][0] += out_prec[j + it->bandIndex * 2]
                * soil_con->Pfactor[it->bandIndex];
            state->debug.outflow[j][it->bandIndex][0] +=
                it->veg_var[j].throughfall;
            if (j == 0)
              state->debug.inflow[j][it->bandIndex][1] += snow_inflow[it->bandIndex];
            state->debug.outflow[j][it->bandIndex][1] += Melt[it->bandIndex * 2 + j];
          }

        }

        //TODO: convert debug struct into local per element information!
        //writeDebug->write_debug(atmos, soil_con, prcp->cell[j][iveg], prcp->energy[iveg],
        //prcp->snow[iveg], prcp->veg_var[j][iveg], &(dmy[time_step_record]), out_short,
        //tmp_mu, Nveg, iveg, time_step_record, j, NEWCELL, state);
      }
#endif // LINK_DEBUG
    } /** end current vegetation type **/
  } /** end of vegetation loop **/

  for (int p = 0; p < N_PET_TYPES + 1; p++) {
    free((char *) aero_resist[p]);
  }
  free((char *) aero_resist);

  /****************************
   Calculate Subsidence
   ****************************/
#if EXCESS_ICE
  total_subsidence = 0;
  total_meltwater = 0; //for lake model only
  for(lidx=0;lidx<state->options.Nlayer;lidx++) { //soil layer
    subsidence[lidx] = 0;
    if(soil_con->effective_porosity[lidx]>soil_con->porosity[lidx]) {

      /* find average ice/porosity fraction and sub-grid with greatest ice/porosity fraction */
      ave_ice = 0;
      max_ice_layer = 0;
      for(iveg = 0; iveg <= Nveg; iveg++) { //iveg
        if (veg_con[iveg].Cv > 0.) {
          Cv = veg_con[iveg].Cv;
          Nbands = state->options.SNOW_BAND;
          if (veg_con[iveg].LAKE) {
            Cv *= (1-lakefrac);
            Nbands = 1;
          }
          for(band = 0; band < Nbands; band++) { //band
            cell_data_struct& cellRef = prcp->cell[dist][iveg][band];
            if(soil_con->AreaFract[band] > 0) {
              for ( dist = 0; dist < Ndist; dist++ ) { // wet/dry
                if(dist==0)
                tmp_mu = prcp->mu[iveg];
                else
                tmp_mu = 1. - prcp->mu[iveg];
#if SPATIAL_FROST
                tmp_ice = 0;
                for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) { //frost area
                  tmp_ice += (cellRef.layer[lidx].ice[frost_area]
                      * soil_con->frost_fract[frost_area]);
                  ice_layer = cellRef.layer[lidx].ice[frost_area];
                  if(ice_layer>=max_ice_layer)
                  max_ice_layer = ice_layer;
                } // frost area
#else //SPATIAL_FROST
                tmp_ice = cellRef.layer[lidx].soil_ice;
                ice_layer = cellRef.layer[lidx].soil_ice;
                if(ice_layer>=max_ice_layer)
                max_ice_layer = ice_layer;
#endif //SPATIAL_FROST
                ave_ice += tmp_ice * Cv * tmp_mu * soil_con->AreaFract[band];
              } // wet/dry
            }
          } //band
        }
      } //iveg
      ave_ice_fract = ave_ice/soil_con->max_moist[lidx];

      /*check to see if threshold is exceeded by average ice/porosity fraction*/
      if(ave_ice_fract <= ICE_AT_SUBSIDENCE) {
        SubsidenceUpdate = 1;

        /*calculate subsidence based on maximum ice content in layer*/
        /*constrain subsidence by MAX_SUBSIDENCE*/
        tmp_depth_prior = soil_con->depth[lidx]; //m
        tmp_subsidence = (1000.*tmp_depth_prior - max_ice_layer);//mm
        if(tmp_subsidence > MAX_SUBSIDENCE)
        tmp_subsidence = MAX_SUBSIDENCE;
        tmp_depth = tmp_depth_prior - tmp_subsidence/1000.;//m
        if(tmp_depth <= soil_con->min_depth[lidx])
        tmp_depth = soil_con->min_depth[lidx];
        soil_con->depth[lidx] = (float)(int)(tmp_depth * 1000 + 0.5) / 1000;//m
        subsidence[lidx] = (tmp_depth_prior - soil_con->depth[lidx])*1000.;//mm
        total_subsidence += (tmp_depth_prior - soil_con->depth[lidx]);//m

        if(subsidence[lidx] > 0 ) {
#if VERBOSE
          fprintf(stderr,"Subsidence of %.3f m in layer %d:\n",subsidence[lidx]/1000.,lidx+1);
          fprintf(stderr,"\t\tOccurred for record=%d: year=%d, month=%d, day=%d, hour=%d\n",rec,dmy[time_step_record].year,dmy[time_step_record].month,dmy[time_step_record].day, dmy[time_step_record].hour);
          fprintf(stderr,"\t\tDepth of soil layer decreased from %.3f m to %.3f m.\n",tmp_depth_prior,soil_con->depth[lidx]);
#endif	  

          /*update soil_con properties*/
#if VERBOSE
          fprintf(stderr,"\t\tEffective porosity decreased from %.3f to %.3f.\n",soil_con->effective_porosity[lidx],1.0-(1.0-soil_con->effective_porosity[lidx])*tmp_depth_prior/soil_con->depth[lidx]);
#endif
          soil_con->effective_porosity[lidx]=1.0-(1.0-soil_con->effective_porosity[lidx])*tmp_depth_prior/soil_con->depth[lidx];
          if(tmp_depth <= soil_con->min_depth[lidx])
          soil_con->effective_porosity[lidx]=soil_con->porosity[lidx];
#if VERBOSE
          fprintf(stderr,"\t\tBulk density increased from %.2f kg/m^3 to %.2f kg/m^3.\n",soil_con->bulk_density[lidx],(1.0-soil_con->effective_porosity[lidx])*soil_con->soil_density[lidx]);
#endif
          soil_con->bulk_dens_min[lidx] *= (1.0-soil_con->effective_porosity[lidx])*soil_con->soil_density[lidx]/soil_con->bulk_density[lidx];
          if (soil_con->organic[lidx] > 0)
          soil_con->bulk_dens_org[lidx] *= (1.0-soil_con->effective_porosity[lidx])*soil_con->soil_density[lidx]/soil_con->bulk_density[lidx];
          soil_con->bulk_density[lidx] = (1.0-soil_con->effective_porosity[lidx])*soil_con->soil_density[lidx]; //adjust bulk density
          total_meltwater += soil_con->max_moist[lidx] - soil_con->depth[lidx] * soil_con->effective_porosity[lidx] * 1000.;//for lake model (uses prior max_moist,
                                                                                                                            //so must come before new max_moist calculation
          soil_con->max_moist[lidx] = soil_con->depth[lidx] * soil_con->effective_porosity[lidx] * 1000.;

        }	                                                   //subsidence occurs
      }	                                                    //threshold exceeded

    }	                                                       //excess ice exists
  } //loop for each soil layer
  if(total_subsidence>0) {

    /********update remaining soil_con properties**********/
#if VERBOSE
    fprintf(stderr,"Damping depth decreased from %.3f m to %.3f m.\n",soil_con->dp,soil_con->dp-total_subsidence);
#endif
    soil_con->dp -= total_subsidence;  //adjust damping depth

#if VERBOSE
    fprintf(stderr,"More updated parameters in soil_con: max_infil, Wcr, and Wpwp\n");
#endif  
    /* update Maximum Infiltration for Upper Layers */
    if(state->options.Nlayer==2)
    soil_con->max_infil = (1.0+soil_con->b_infilt)*soil_con->max_moist[0];
    else
    soil_con->max_infil = (1.0+soil_con->b_infilt)*(soil_con->max_moist[0]+soil_con->max_moist[1]);

    /* Soil Layer Critical and Wilting Point Moisture Contents */
    for(lidx=0;lidx<state->options.Nlayer;lidx++) {  //soil layer
      soil_con->Wcr[lidx] = soil_con->Wcr_FRACT[lidx] * soil_con->max_moist[lidx];
      soil_con->Wpwp[lidx] = soil_con->Wpwp_FRACT[lidx] * soil_con->max_moist[lidx];
      if(soil_con->Wpwp[lidx] > soil_con->Wcr[lidx]) {
        sprintf(ErrStr,"Updated wilting point moisture (%f mm) is greater than updated critical point moisture (%f mm) for layer %d.\n\tIn the soil parameter file, Wpwp_FRACT MUST be <= Wcr_FRACT.\n",
            soil_con->Wpwp[lidx], soil_con->Wcr[lidx], lidx);
        nrerror(ErrStr);
      }
      if(soil_con->Wpwp[lidx] < soil_con->resid_moist[lidx] * soil_con->depth[lidx] * 1000.) {
        sprintf(ErrStr,"Updated wilting point moisture (%f mm) is less than updated residual moisture (%f mm) for layer %d.\n\tIn the soil parameter file, Wpwp_FRACT MUST be >= resid_moist / (1.0 - bulk_density/soil_density).\n",
            soil_con->Wpwp[lidx], soil_con->resid_moist[lidx] * soil_con->depth[lidx] * 1000., lidx);
        nrerror(ErrStr);
      }
    }

    /* If BASEFLOW = NIJSSEN2001 then convert ARNO baseflow
     parameters d1, d2, d3, and d4 to Ds, Dsmax, Ws, and c */
    if(state->options.BASEFLOW == NIJSSEN2001) {
      lidx = state->options.Nlayer-1;
      soil_con->Dsmax = soil_con->Dsmax_orig *
      pow((double)(1./(soil_con->max_moist[lidx]-soil_con->Ws_orig)), -soil_con->c) +
      soil_con->Ds_orig * soil_con->max_moist[lidx];
      soil_con->Ds = soil_con->Ds_orig * soil_con->Ws_orig / soil_con->Dsmax_orig;
      soil_con->Ws = soil_con->Ws_orig/soil_con->max_moist[lidx];
#if VERBOSE
      fprintf(stderr,"More updated parameters in soil_con: Dsmax, Ds, Ws\n");
#endif    
    }

    /*********** update root fractions ***************/
    fprintf(stderr,"Updated parameter in veg_con: root\n");
    calc_root_fractions(veg_con, soil_con);

    /**********redistribute soil moisture (call runoff function)*************/
    /* If subsidence occurs, recalculate runoff, baseflow, and soil moisture
     using soil moisture values from previous time-step; i.e.
     as if prior runoff call did not occur.*/
    for(iveg = 0; iveg <= Nveg; iveg++) {
      if (veg_con[iveg].Cv > 0.) {
        Nbands = state->options.SNOW_BAND;
        if (veg_con[iveg].LAKE) {
          Nbands = 1;
        }
        for ( band = 0; band < Nbands; band++ ) {
          for ( dist = 0; dist < Ndist; dist++ ) {
            for(lidx=0;lidx<state->options.Nlayer;lidx++) {
              prcp->cell[dist][iveg][band].layer[lidx].moist = moist_prior[dist][iveg][band][lidx];
              prcp->cell[dist][iveg][band].layer[lidx].evap = evap_prior[dist][iveg][band][lidx];
            }
          }
        }
      }
    }
    for(iveg = 0; iveg <= Nveg; iveg++) {
      if (veg_con[iveg].Cv > 0.) {
        Nbands = state->options.SNOW_BAND;
        if (veg_con[iveg].LAKE) {
          Nbands = 1;
        }
        for ( band = 0; band < Nbands; band++ ) {
          if( soil_con->AreaFract[band] > 0 ) {

            //set inflow for runoff call
            ppt[WET]=prcp->cell[WET][iveg][band].inflow;
            ppt[DRY]=prcp->cell[DRY][iveg][band].inflow;

            ErrorFlag = runoff(&(prcp->cell[WET][iveg][band]), &(prcp->cell[DRY][iveg][band]), &(prcp->energy[iveg][band]),
                soil_con, ppt,
                SubsidenceUpdate,
                soil_con->frost_fract,
                prcp->mu[iveg], state->global_param.dt, state->options.Nnode, band, time_step_record, iveg);
            if ( ErrorFlag == ERROR ) return ( ERROR );

          }
        }  //band
      }
    }  //veg

    /**********interpolate nodal temperatures to new depths and recalculate thermal properties***********/
    ErrorFlag = update_thermal_nodes(prcp, Nveg, state->options.Nnode, Ndist, soil_con, veg_con);
    if ( ErrorFlag == ERROR ) return ( ERROR );

  }  //subsidence occurs

  /********************************************
   Save subsidence for output
   ********************************************/
  for(lidx=0;lidx<state->options.Nlayer;lidx++)
  soil_con->subsidence[lidx] = subsidence[lidx];

#endif //EXCESS_ICE
  /****************************
   Run Lake Model
   ****************************/

  /** Compute total runoff and baseflow for all vegetation types
   within each snowband. **/
  if (state->options.LAKES && lake_con->lake_idx >= 0) {

    wetland_runoff = wetland_baseflow = 0;
    sum_runoff = sum_baseflow = 0;

    // Loop through all vegetation tiles
    for (std::vector<HRU>::iterator it = prcp->hruList.begin(); it != prcp->hruList.end(); ++it) {

      /** Solve Veg Tile only if Coverage Greater than 0% **/
      if (veg_con[it->vegIndex].Cv > 0.) {

        Cv = veg_con[it->vegIndex].Cv;
        Nbands = state->options.SNOW_BAND;
        if (veg_con[it->vegIndex].LAKE) {
          Cv *= (1 - lakefrac);
          Nbands = 1;
        }

        // For each of the prcp->snow elevation bands
        if (soil_con->AreaFract[it->bandIndex] > 0) {

          // Loop through distributed precipitation fractions
          for (int dist = 0; dist < 2; dist++) {
            cell_data_struct& cellRef = it->cell[dist];
            if (dist == 0)
              tmp_mu = prcp->mu[it->vegIndex];
            else
              tmp_mu = 1. - prcp->mu[it->vegIndex];

            if (veg_con[it->vegIndex].LAKE) {
              wetland_runoff += (cellRef.runoff * tmp_mu * Cv
                  * soil_con->AreaFract[it->bandIndex]);
              wetland_baseflow += (cellRef.baseflow * tmp_mu * Cv
                  * soil_con->AreaFract[it->bandIndex]);
              cellRef.runoff = 0;
              cellRef.baseflow = 0;
            } else {
              sum_runoff += (cellRef.runoff * tmp_mu * Cv
                  * soil_con->AreaFract[it->bandIndex]);
              sum_baseflow += (cellRef.baseflow * tmp_mu * Cv
                  * soil_con->AreaFract[it->bandIndex]);
              cellRef.runoff *= (1 - lake_con->rpercent);
              cellRef.baseflow *= (1 - lake_con->rpercent);
            }
          }
        }
      }
    }

    /** Run lake model **/
    prcp->lake_var.runoff_in = (sum_runoff * lake_con->rpercent + wetland_runoff)
        * soil_con->cell_area * 0.001; // m3
    prcp->lake_var.baseflow_in = (sum_baseflow * lake_con->rpercent
        + wetland_baseflow) * soil_con->cell_area * 0.001; // m3
    prcp->lake_var.channel_in = atmos->channel_in[state->NR] * soil_con->cell_area * 0.001; // m3
    prcp->lake_var.prec = atmos->prec[state->NR] * prcp->lake_var.sarea * 0.001; // m3
    rainonly = calc_rainonly(atmos->air_temp[state->NR], atmos->prec[state->NR], state->global_param.MAX_SNOW_TEMP, state->global_param.MIN_RAIN_TEMP, 1);
    if ((int) rainonly == ERROR) {
      return (ERROR);
    }

    /**********************************************************************
       Solve the energy budget for the lake.
     **********************************************************************/

    oldsnow = prcp->lake_var.snow.swq;
    snowprec = gauge_correction[SNOW] * (atmos->prec[state->NR] - rainonly);
    rainprec = gauge_correction[SNOW] * rainonly;
    atmos->out_prec += (snowprec + rainprec) * lake_con->Cl[0] * lakefrac;

    ErrorFlag = solve_lake(snowprec, rainprec, atmos->air_temp[state->NR],
        atmos->wind[state->NR], atmos->vp[state->NR] / 1000., atmos->shortwave[state->NR],
        atmos->longwave[state->NR], atmos->vpd[state->NR] / 1000.,
        atmos->pressure[state->NR] / 1000., atmos->density[state->NR], &(prcp->lake_var), *lake_con,
        *soil_con, state->global_param.dt, time_step_record, state->global_param.wind_h, dmy[time_step_record], fraci, state);
    if (ErrorFlag == ERROR)
      return (ERROR);

    /**********************************************************************
       Solve the water budget for the lake.
     **********************************************************************/

    ErrorFlag = water_balance(&(prcp->lake_var), *lake_con,
        state->global_param.dt, prcp, time_step_record,
        prcp->getHRUElement(lake_con->lake_idx, 0), lakefrac, *soil_con,
        *veg_con, SubsidenceUpdate, total_meltwater, state);
    if (ErrorFlag == ERROR)
      return (ERROR);

#if LINK_DEBUG
    if ( state->debug.PRT_LAKE ) {
      if ( time_step_record == 0 ) {
        // print file header
        fprintf(state->debug.fg_lake,"Date,Rec,AeroResist,BasinflowIn,BaseflowOut");
        for (int i = 0; i < MAX_LAKE_NODES; i++ )
        fprintf(state->debug.fg_lake, ",Density%i", i);
        fprintf(state->debug.fg_lake,",Evap,IceArea,IceHeight,LakeDepth,RunoffIn,RunoffOut,SurfaceArea,SnowDepth,SnowMelt");
        for (int i = 0; i < MAX_LAKE_NODES; i++ )
        fprintf(state->debug.fg_lake, ",Area%i", i);
        fprintf(state->debug.fg_lake,",SWE");
        for (int i = 0; i < MAX_LAKE_NODES; i++ )
        fprintf(state->debug.fg_lake, ",Temp%i", i);
        fprintf(state->debug.fg_lake,",IceTemp,Volume,Nodes");
        fprintf(state->debug.fg_lake,",AlbedoLake,AlbedoOver,AlbedoUnder,AtmosError,AtmosLatent,AtmosLatentSub,AtmosSensible,LongOverIn,LongUnderIn,LongUnderOut,NetLongAtmos,NetLongOver,NetLongUnder,NetShortAtmos,NetShortGrnd,NetShortOver,NetShortUnder");
        fprintf(state->debug.fg_lake,",ShortOverIn,ShortUnderIn,advection,deltaCC,deltaH,error,fusion,grnd_flux,latent,latent_sub,longwave,melt_energy,out_long_surface,refreeze_energy,sensible,shortwave");
        fprintf(state->debug.fg_lake,",Qnet,albedo,coldcontent,coverage,density,depth,mass_error,max_swq,melt,pack_temp,pack_water,store_coverage,store_swq,surf_temp,surf_water,swq,swq_slope,vapor_flux,last_snow,store_snow\n");
      }

      // print lake variables
      fprintf(state->debug.fg_lake, "%i/%i/%i %i:00:00,%i", dmy[time_step_record].month, dmy[time_step_record].day, dmy[time_step_record].year, dmy[time_step_record].hour, time_step_record);
      fprintf(state->debug.fg_lake, ",%i", prcp->lake_var.activenod);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.ldepth);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.sarea);
      for (int i = 0; i < MAX_LAKE_NODES; i++ )
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.surface[i]);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.volume);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.baseflow_in);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.baseflow_out);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.channel_in);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.evapw);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.prec);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.recharge);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.runoff_in);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.runoff_out);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.snowmlt);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.vapor_flux);
      for (int i = 0; i < MAX_LAKE_NODES; i++ )
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.temp[i]);
      for (int i = 0; i < MAX_LAKE_NODES; i++ )
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.density[i]);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.areai);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.hice);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.sdepth);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.swe);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.tempi);
      fprintf(state->debug.fg_lake, ",%f", prcp->lake_var.aero_resist);

      energy_bal_struct& energyRef = prcp->getHRUElement(lake_con->lake_idx, 0)->energy;
      // print lake prcp->energy variables
      fprintf(state->debug.fg_lake, ",%f", energyRef.AlbedoLake);
      fprintf(state->debug.fg_lake, ",%f", energyRef.AlbedoOver);
      fprintf(state->debug.fg_lake, ",%f", energyRef.AlbedoUnder);
      fprintf(state->debug.fg_lake, ",%f", energyRef.AtmosError);
      fprintf(state->debug.fg_lake, ",%f", energyRef.AtmosLatent);
      fprintf(state->debug.fg_lake, ",%f", energyRef.AtmosLatentSub);
      fprintf(state->debug.fg_lake, ",%f", energyRef.AtmosSensible);
      fprintf(state->debug.fg_lake, ",%f", energyRef.LongOverIn);
      fprintf(state->debug.fg_lake, ",%f", energyRef.LongUnderIn);
      fprintf(state->debug.fg_lake, ",%f", energyRef.LongUnderOut);
      fprintf(state->debug.fg_lake, ",%f", energyRef.NetLongAtmos);
      fprintf(state->debug.fg_lake, ",%f", energyRef.NetLongOver);
      fprintf(state->debug.fg_lake, ",%f", energyRef.NetLongUnder);
      fprintf(state->debug.fg_lake, ",%f", energyRef.NetShortAtmos);
      fprintf(state->debug.fg_lake, ",%f", energyRef.NetShortGrnd);
      fprintf(state->debug.fg_lake, ",%f", energyRef.NetShortOver);
      fprintf(state->debug.fg_lake, ",%f", energyRef.NetShortUnder);
      fprintf(state->debug.fg_lake, ",%f", energyRef.ShortOverIn);
      fprintf(state->debug.fg_lake, ",%f", energyRef.ShortUnderIn);
      fprintf(state->debug.fg_lake, ",%f", energyRef.advection);
      fprintf(state->debug.fg_lake, ",%f", energyRef.deltaCC);
      fprintf(state->debug.fg_lake, ",%f", energyRef.deltaH);
      fprintf(state->debug.fg_lake, ",%f", energyRef.error);
      fprintf(state->debug.fg_lake, ",%f", energyRef.fusion);
      fprintf(state->debug.fg_lake, ",%f", energyRef.grnd_flux);
      fprintf(state->debug.fg_lake, ",%f", energyRef.latent);
      fprintf(state->debug.fg_lake, ",%f", energyRef.latent_sub);
      fprintf(state->debug.fg_lake, ",%f", energyRef.longwave);
      fprintf(state->debug.fg_lake, ",%f", energyRef.melt_energy);
      fprintf(state->debug.fg_lake, ",%f", energyRef.out_long_surface);
      fprintf(state->debug.fg_lake, ",%f", energyRef.refreeze_energy);
      fprintf(state->debug.fg_lake, ",%f", energyRef.sensible);
      fprintf(state->debug.fg_lake, ",%f", energyRef.shortwave);

      snow_data_struct& snowRef = prcp->getHRUElement(lake_con->lake_idx, 0)->snow;
      // print lake prcp->snow variables
      fprintf(state->debug.fg_lake, ",%f", snowRef.Qnet);
      fprintf(state->debug.fg_lake, ",%f", snowRef.albedo);
      fprintf(state->debug.fg_lake, ",%f", snowRef.coldcontent);
      fprintf(state->debug.fg_lake, ",%f", snowRef.coverage);
      fprintf(state->debug.fg_lake, ",%f", snowRef.density);
      fprintf(state->debug.fg_lake, ",%f", snowRef.depth);
      fprintf(state->debug.fg_lake, ",%f", snowRef.mass_error);
      fprintf(state->debug.fg_lake, ",%f", snowRef.max_swq);
      fprintf(state->debug.fg_lake, ",%f", snowRef.melt);
      fprintf(state->debug.fg_lake, ",%f", snowRef.pack_temp);
      fprintf(state->debug.fg_lake, ",%f", snowRef.pack_water);
      fprintf(state->debug.fg_lake, ",%f", snowRef.store_coverage);
      fprintf(state->debug.fg_lake, ",%f", snowRef.store_swq);
      fprintf(state->debug.fg_lake, ",%f", snowRef.surf_temp);
      fprintf(state->debug.fg_lake, ",%f", snowRef.surf_water);
      fprintf(state->debug.fg_lake, ",%f", snowRef.swq);
      fprintf(state->debug.fg_lake, ",%f", snowRef.swq_slope);
      fprintf(state->debug.fg_lake, ",%f", snowRef.vapor_flux);
      fprintf(state->debug.fg_lake, ",%i", snowRef.last_snow);
      fprintf(state->debug.fg_lake, ",%i\n", snowRef.store_snow);

    }
#endif // LINK_DEBUG
  } // end if (state->options.LAKES && lake_con->lake_idx >= 0)

  return (0);
}

