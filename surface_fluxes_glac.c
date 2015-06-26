#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include <math.h>

int surface_fluxes_glac(
       double               BareAlbedo,
       double               height,
       double               ice0,
       double               moist0,
       int                  SubsidenceUpdate,
       double              *evap_prior_dry,
       double              *evap_prior_wet,
       HRU&                 hru,
       double              *Melt,
       double              *latent_heat_Le,
       VegConditions       *aero_resist,
       VegConditions       &displacement,
       double              *gauge_correction,
       double              *out_prec,
       double              *out_rain,
       double              *out_snow,
       VegConditions       &ref_height,
       VegConditions       &roughness,
       double              *snow_inflow,
       VegConditions       &wind_speed,
       int                  Nbands,
       int                  Ndist,
       int                  Nlayers,
       int                  rec,
       int                  veg_class,
       atmos_data_struct   *atmos,
       const dmy_struct    *dmy,
       const soil_con_struct *soil_con,
       float              lag_one,
       float              sigma_slope,
       float              fetch,
       const ProgramState *state)
/**********************************************************************

  This is a modified version of surface_fluxes.c specific for glaciers.
  This routine computes all surface fluxes, and solves the snow accumulation
  and ablation algorithm. Solutions are for the current snow band and
  vegetation type (these are defined in full_energy before the routine is called).

**********************************************************************/
{
  int                    ErrorFlag;
  int                    N_steps;
  VegConditions::VegSurfType UnderStory;
  int                    hidx;     // index of initial element of atmos array
  int                    step_inc; // number of atmos array elements to skip per surface fluxes step
  int                    endhidx;  // index of final element of atmos array
  int                    step_dt;  // time length of surface fluxes step
  double                 LongUnderIn; // incoming LW to ground surface
  double                 LongUnderOut; // outgoing LW from ground surface
  double                 NetLongSnow; // net LW over snowpack or glacier
  double                 NetShortSnow; // net SW over snowpqack or glacier
  double                 NetShortGrnd; // net SW over snow-free surface
  double                 OldTSurf; // previous snow surface temperature
  double                 ShortUnderIn; // incoming SW to understory
  double                 Tair; // air temperature
  double                 Tcanopy; // canopy air temperature
  double                 Tgrnd; // soil surface temperature
  double                 Tsurf; // ground surface temperature
  double                 VPDcanopy; // vapor pressure deficit in canopy/atmosphere
  double                 coverage; // mid-step snow cover fraction
  double                 delta_coverage; // change in snow cover fraction
  double                 ppt[2]; // precipitation/melt reaching soil surface or glacier surface
  double                 rainfall[2]; // rainfall
  double                 snowfall[2]; // snowfall
  double                 snow_flux; // heat flux through snowpack
  double                 rainOnly;

  // Step-specific quantities
  double                 step_Wdew[2];
  double                 step_melt;
  double                 step_melt_energy;  /* energy used to reduce snow coverage */
  double                 step_out_prec;
  double                 step_out_rain;
  double                 step_out_snow;
  double                 step_ppt[2];
  double                 step_prec[2];
  AeroResistUsed        *step_aero_resist;
  double step_melt_glac;

  // Quantities that need to be summed or averaged over multiple snow steps
  // energy structure
  double                 store_AlbedoOver = 0;
  double                 store_AlbedoUnder = 0;
  double                 store_AtmosLatent = 0;
  double                 store_AtmosLatentSub = 0;
  double                 store_AtmosSensible = 0;
  double                 store_LongOverIn = 0;
  double                 store_LongUnderIn = 0;
  double                 store_LongUnderOut = 0;
  double                 store_NetLongAtmos = 0;
  double                 store_NetLongOver = 0;
  double                 store_NetLongUnder = 0;
  double                 store_NetShortAtmos = 0;
  double                 store_NetShortGrnd = 0;
  double                 store_NetShortOver = 0;
  double                 store_NetShortUnder = 0;
  double                 store_ShortOverIn = 0;
  double                 store_ShortUnderIn = 0;
  double                 store_advected_sensible = 0;
  double                 store_advection = 0;
  double                 store_canopy_advection = 0;
  double                 store_canopy_latent = 0;
  double                 store_canopy_latent_sub = 0;
  double                 store_canopy_sensible = 0;
  double                 store_canopy_refreeze = 0;
  double                 store_deltaCC = 0;
  double                 store_deltaH = 0;
  double                 store_fusion = 0;
  double                 store_grnd_flux = 0;
  double                 store_latent = 0;
  double                 store_latent_sub = 0;
  double                 store_melt_energy = 0;
  double                 store_refreeze_energy = 0;
  double                 store_sensible = 0;
  double                 store_snow_flux = 0;
  double                 store_deltaCC_glac = 0;
  double                 store_glacier_flux = 0;
  // glacier structure
  double                 store_melt_glac = 0;
  double                 store_vapor_flux_glac = 0;
  double                 store_accum_glac = 0;
  // snow structure
  double                 store_canopy_vapor_flux = 0;
  double                 store_melt = 0;
  double                 store_vapor_flux = 0;
  double                 store_blowing_flux = 0;
  double                 store_surface_flux = 0;
  // veg_var structure
  double                 store_canopyevap[2];
  double                 store_throughfall[2];
  // cell structure
  double                 store_layerevap[2][MAX_LAYERS];
  double                 store_ppt[2];
  AeroResistUsed         store_aero_cond_used;
  double                 store_pot_evap[N_PET_TYPES];

  // Structures holding values for current snow step
  energy_bal_struct      step_energy; // energy fluxes at snowpack surface and glacier surface
  veg_var_struct         snow_veg_var[2]; // veg fluxes/storages in presence of snow
  veg_var_struct         soil_veg_var[2]; // veg fluxes/storages in soil energy balance
  snow_data_struct       step_snow;
  layer_data_struct      step_layer[2][MAX_LAYERS];
  glac_data_struct       step_glacier;

  // Structures holding values for current iteration
  VegConditions          temp_aero_resist;
  AeroResistUsed         temp_aero_resist_used;
  double                 stability_factor[2];
  double                 step_pot_evap[N_PET_TYPES];

  step_aero_resist = new AeroResistUsed[N_PET_TYPES];

  /***********************************************************************
   Set temporary variables - preserves original values until iterations
   are completed
   ***********************************************************************/

  coverage = hru.snow.coverage;
  step_energy = hru.energy;
  snow_veg_var[WET] = hru.veg_var[WET];
  snow_veg_var[DRY] = hru.veg_var[DRY];
  soil_veg_var[WET] = hru.veg_var[WET];
  soil_veg_var[DRY] = hru.veg_var[DRY];
  step_snow = hru.snow;
  step_glacier = hru.glacier;
  for (int lidx = 0; lidx < Nlayers; lidx++) {
    step_layer[WET][lidx] = hru.cell[WET].layer[lidx];
    step_layer[DRY][lidx] = hru.cell[DRY].layer[lidx];
  }
  for (int lidx = 0; lidx < Nlayers; lidx++) {
    step_layer[WET][lidx].evap = 0;
    step_layer[DRY][lidx].evap = 0;
  }
  soil_veg_var[WET].canopyevap = 0;
  soil_veg_var[DRY].canopyevap = 0;
  snow_veg_var[WET].canopyevap = 0;
  snow_veg_var[DRY].canopyevap = 0;
  soil_veg_var[WET].throughfall = 0;
  soil_veg_var[DRY].throughfall = 0;
  snow_veg_var[WET].throughfall = 0;
  snow_veg_var[DRY].throughfall = 0;

  /********************************
   Set-up sub-time step controls
   (May eventually want to set this up so that it is also true
   if frozen soils are present)
   ********************************/

  // Always use sub-timestep for snow + glaciers.
    hidx = 0;
    step_inc = 1;
    endhidx = hidx + state->NF;
    step_dt = state->options.SNOW_STEP;

  /*******************************************
   Initialize sub-model time step variables
   *******************************************/


  // veg_var and cell structures
  for (int dist = 0; dist < Ndist; dist++) {
    store_throughfall[dist] = 0.;
    store_canopyevap[dist] = 0.;
    for (int lidx = 0; lidx < state->options.Nlayer; lidx++) {
      store_layerevap[dist][lidx] = 0.;
    }
  }
  step_Wdew[WET] = hru.veg_var[WET].Wdew;
  step_Wdew[DRY] = hru.veg_var[WET].Wdew;
  // misc
  store_ppt[WET] = 0;
  store_ppt[DRY] = 0;
  step_prec[DRY] = 0;
  store_aero_cond_used.surface = 0;
  store_aero_cond_used.overstory = 0;
  (*snow_inflow) = 0;
  for (int p = 0; p < N_PET_TYPES; p++) {
    store_pot_evap[p] = 0;
  }
  N_steps = 0;

  /*************************
   Compute surface fluxes
   *************************/

  do { /* MPN TODO This is stupid.  Makes it into a for loop because it is a fixed number of iterations! */

    /** Solve energy balance for all sub-model time steps **/

    /* set air temperature and precipitation for this snow band */
    Tair = atmos->air_temp[hidx] + soil_con->Tfactor[hru.bandIndex];
    step_prec[WET] = atmos->prec[hidx] / hru.mu * soil_con->Pfactor[hru.bandIndex];  // precipitation in mm

    rainOnly = calc_rainonly(Tair, step_prec[WET], soil_con->MAX_SNOW_TEMP,
        soil_con->MIN_RAIN_TEMP, hru.mu);
    snowfall[WET] = gauge_correction[SNOW] * (step_prec[WET] - rainOnly);
    rainfall[WET] = gauge_correction[RAIN] * rainOnly;
    snowfall[DRY] = 0.;
    rainfall[DRY] = 0.;

    if(snowfall[WET] < 1e-5) snowfall[WET] = 0;

    step_out_prec = snowfall[WET] + rainfall[WET];
    step_out_rain = rainfall[WET];
    step_out_snow = snowfall[WET];

    // initialize ground surface temperature: set to glacier temperature
    Tgrnd = GLAC_TEMP;

    // initialize canopy terms
    Tcanopy = 0.;
    VPDcanopy = 0.;

    // Compute mass flux of blowing snow
    if (state->options.BLOWING && step_snow.swq > 0.) {
      double Ls = (677. - 0.07 * step_snow.surf_temp) * JOULESPCAL * GRAMSPKG;
      step_snow.blowing_flux = CalcBlowingSnow((double) step_dt, Tair,
          step_snow.last_snow, step_snow.surf_water, wind_speed.snowCovered, Ls,
          atmos->density[hidx], atmos->pressure[hidx], atmos->vp[hidx],
          roughness.snowCovered, ref_height.snowCovered, step_snow.depth, lag_one, sigma_slope,
          step_snow.surf_temp, hru.isArtificialBareSoil, fetch, displacement.canopyIfOverstory, roughness.canopyIfOverstory,
          &step_snow.transport);
      if ((int) step_snow.blowing_flux == ERROR) {
        return (ERROR);
      }
      step_snow.blowing_flux *= step_dt * SECPHOUR / RHO_W; /* m/time step */
    } else
      step_snow.blowing_flux = 0.0;

    temp_aero_resist = aero_resist[N_PET_TYPES];
    temp_aero_resist_used.surface = hru.cell[WET].aero_resist.surface;
    temp_aero_resist_used.overstory = hru.cell[WET].aero_resist.overstory;
    step_snow.canopy_vapor_flux = 0;
    step_snow.vapor_flux = 0;
    step_snow.surface_flux = 0;

    LongUnderOut = step_energy.LongUnderOut;

    if (step_snow.swq > 0. || snowfall[WET] > 0.) {
      /** Solve snow accumulation and ablation on the glacier surface **/

      step_melt = solve_snow_glac(BareAlbedo, Tgrnd, Tair, hru.mu,
       &hru.energy.AlbedoUnder, latent_heat_Le, &LongUnderIn, &NetLongSnow,
       &NetShortSnow, &ShortUnderIn, &OldTSurf, temp_aero_resist,
       temp_aero_resist_used, &coverage, &delta_coverage,
       displacement, &step_melt_energy, step_ppt, rainfall, ref_height,
       roughness, snow_inflow, snowfall, wind_speed, step_dt, rec, hidx,
       UnderStory, dmy, *atmos, &(step_energy), &(step_snow), soil_con,
       &step_glacier, state);

      if (step_melt == ERROR)
        return (ERROR);

      step_melt_glac = 0.;
      step_glacier.vapor_flux = 0.;
      step_energy.glacier_flux = 0.;
      step_energy.deltaCC_glac = 0.;
      step_energy.snow_flux = -step_energy.grnd_flux;
      step_energy.LongUnderOut = LongUnderIn - NetLongSnow;

    } else {

      step_melt_glac = solve_glacier(BareAlbedo, Tgrnd, Tair,
    	  &hru.energy.AlbedoUnder, latent_heat_Le, &LongUnderIn,
    	  &NetLongSnow, &NetShortSnow, &ShortUnderIn, &OldTSurf,
    	  temp_aero_resist, temp_aero_resist_used, displacement,
    	  &step_melt_energy, step_ppt, rainfall, ref_height, roughness,
          wind_speed, step_dt, rec, hidx, UnderStory, atmos,
          &step_energy, &step_glacier, soil_con, state);

      if (step_melt_glac == ERROR) {
        return (ERROR);
      }

      step_melt = 0.;
      hru.snow.snow = FALSE;
      hru.snow.store_swq = 0.;
      hru.snow.store_coverage = 1;
      hru.snow.last_snow = INVALID_INT;
      hru.snow.albedo = soil_con->NEW_SNOW_ALB;
      step_energy.deltaCC = 0.;
      step_energy.refreeze_energy = 0.;
      step_energy.snow_flux = 0.;
      step_energy.advected_sensible = 0.;
      step_energy.glacier_flux = -step_energy.grnd_flux;
      step_energy.LongUnderOut = LongUnderIn - NetLongSnow;
      step_glacier.accumulation = 0.;

    }

    /**************************************
     Compute Potential Evap
     **************************************/
    // First, determine the stability correction used in the iteration
    if (temp_aero_resist_used.surface == HUGE_RESIST)
      stability_factor[0] = HUGE_RESIST;
    else
      stability_factor[0] = temp_aero_resist_used.surface
          / aero_resist[N_PET_TYPES][UnderStory];
    if (temp_aero_resist_used.overstory == temp_aero_resist_used.surface)
      stability_factor[1] = stability_factor[0];
    else {
      if (temp_aero_resist_used.overstory == HUGE_RESIST)
        stability_factor[1] = HUGE_RESIST;
      else
        stability_factor[1] = temp_aero_resist_used.overstory
            / aero_resist[N_PET_TYPES].canopyIfOverstory;
    }

    // Next, loop over pot_evap types and apply the correction to the relevant aerodynamic resistance
    for (int p = 0; p < N_PET_TYPES; p++) {
      if (stability_factor[0] == HUGE_RESIST)
        step_aero_resist[p].surface = HUGE_RESIST;
      else
        step_aero_resist[p].surface = aero_resist[p][UnderStory]
            * stability_factor[0];
      if (stability_factor[1] == HUGE_RESIST)
        step_aero_resist[p].overstory = HUGE_RESIST;
      else
        step_aero_resist[p].overstory = aero_resist[p].canopyIfOverstory * stability_factor[1];
    }

    // Finally, compute pot_evap
    compute_pot_evap(veg_class, dmy, rec, state->global_param.dt,
        atmos->shortwave[hidx], step_energy.NetLongAtmos, Tair, VPDcanopy,
        soil_con->elevation, step_aero_resist, step_pot_evap, state);

    /**************************************
     Store sub-model time step variables
     **************************************/

    for (int dist = 0; dist < Ndist; dist++) {
      store_ppt[dist] += step_ppt[dist];
    }
    if (temp_aero_resist_used.surface > 0)
      store_aero_cond_used.surface += 1 / temp_aero_resist_used.surface;
    else
      store_aero_cond_used.surface += HUGE_RESIST;
    if (temp_aero_resist_used.overstory > 0)
      store_aero_cond_used.overstory += 1 / temp_aero_resist_used.overstory;
    else
      store_aero_cond_used.overstory += HUGE_RESIST;

    store_melt += step_melt;
    store_vapor_flux += step_snow.vapor_flux;
    store_surface_flux += step_snow.surface_flux;
    store_blowing_flux += step_snow.blowing_flux;

    out_prec[0] += step_out_prec * hru.mu;
    out_rain[0] += step_out_rain * hru.mu;
    out_snow[0] += step_out_snow * hru.mu;

    store_AlbedoUnder += step_energy.AlbedoUnder;
    store_AtmosLatent += step_energy.AtmosLatent;
    store_AtmosLatentSub += step_energy.AtmosLatentSub;
    store_AtmosSensible += step_energy.AtmosSensible;
    store_LongUnderIn += LongUnderIn;
    store_LongUnderOut += step_energy.LongUnderOut;
    store_NetLongAtmos += NetLongSnow;
    store_NetLongUnder += NetLongSnow;
    store_NetShortAtmos += NetShortSnow;
    store_NetShortUnder += NetShortSnow;
    store_ShortUnderIn += ShortUnderIn;
    store_latent += step_energy.latent;
    store_latent_sub += step_energy.latent_sub;
    store_melt_energy += step_melt_energy;
    store_sensible += step_energy.sensible;
    store_grnd_flux += step_energy.grnd_flux;
    //*****************************************************

    // glacier
    store_melt_glac += step_melt_glac;
    store_vapor_flux_glac += step_glacier.vapor_flux;
    store_accum_glac += step_glacier.accumulation;
    store_glacier_flux += step_energy.glacier_flux;
    store_deltaCC_glac += step_energy.deltaCC_glac;

    store_advected_sensible += step_energy.advected_sensible
        * (step_snow.coverage + delta_coverage);
    store_advection += step_energy.advection
        * (step_snow.coverage + delta_coverage);
    store_deltaCC += step_energy.deltaCC
        * (step_snow.coverage + delta_coverage);
    store_snow_flux += step_energy.snow_flux
        * (step_snow.coverage + delta_coverage);
    store_refreeze_energy += step_energy.refreeze_energy
        * (step_snow.coverage + delta_coverage);

    for (int p = 0; p < N_PET_TYPES; p++) {
      store_pot_evap[p] += step_pot_evap[p];
    }

    /* increment time step */
    N_steps++;
    hidx += step_inc;

  } while (hidx < endhidx);

  /************************************************
   Store glacier variables for sub-model time steps
   ************************************************/
  hru.glacier = step_glacier;
  hru.glacier.melt = store_melt_glac;
  hru.glacier.vapor_flux = store_vapor_flux_glac;
  hru.glacier.accumulation = store_accum_glac;

  /************************************************
   Store snow variables for sub-model time steps
   ************************************************/

  hru.snow = step_snow;
  hru.snow.vapor_flux = store_vapor_flux;
  hru.snow.blowing_flux = store_blowing_flux;
  hru.snow.surface_flux = store_surface_flux;
  hru.snow.canopy_vapor_flux = store_canopy_vapor_flux;
  (*Melt) = store_melt + store_melt_glac;
  hru.snow.melt = store_melt;
  for (int dist = 0; dist < 2; dist++) {
    ppt[dist] = store_ppt[dist];
  }

  hru.glacier.mass_balance = out_prec[WET]/1000. - ppt[WET] - hru.snow.vapor_flux - hru.glacier.vapor_flux; /* convert out_prec[WET] to m */
  hru.glacier.ice_mass_balance = hru.glacier.accumulation - hru.glacier.melt - hru.glacier.vapor_flux;

  /******************************************************
   Store energy flux averages for sub-model time steps
   ******************************************************/

  hru.energy = step_energy;
  hru.energy.AlbedoOver = store_AlbedoOver / (double) N_steps;
  hru.energy.AlbedoUnder = store_AlbedoUnder / (double) N_steps;
  hru.energy.AtmosLatent = store_AtmosLatent / (double) N_steps;
  hru.energy.AtmosLatentSub = store_AtmosLatentSub / (double) N_steps;
  hru.energy.AtmosSensible = store_AtmosSensible / (double) N_steps;
  hru.energy.LongOverIn = store_LongOverIn / (double) N_steps;
  hru.energy.LongUnderIn = store_LongUnderIn / (double) N_steps;
  hru.energy.LongUnderOut = store_LongUnderOut / (double) N_steps;
  hru.energy.NetLongAtmos = store_NetLongAtmos / (double) N_steps;
  hru.energy.NetLongOver = store_NetLongOver / (double) N_steps;
  hru.energy.NetLongUnder = store_NetLongUnder / (double) N_steps;
  hru.energy.NetShortAtmos = store_NetShortAtmos / (double) N_steps;
  hru.energy.NetShortGrnd = store_NetShortGrnd / (double) N_steps;
  hru.energy.NetShortOver = store_NetShortOver / (double) N_steps;
  hru.energy.NetShortUnder = store_NetShortUnder / (double) N_steps;
  hru.energy.ShortOverIn = store_ShortOverIn / (double) N_steps;
  hru.energy.ShortUnderIn = store_ShortUnderIn / (double) N_steps;
  hru.energy.advected_sensible = store_advected_sensible / (double) N_steps;
  hru.energy.canopy_advection = store_canopy_advection / (double) N_steps;
  hru.energy.canopy_latent = store_canopy_latent / (double) N_steps;
  hru.energy.canopy_latent_sub = store_canopy_latent_sub / (double) N_steps;
  hru.energy.canopy_refreeze = store_canopy_refreeze / (double) N_steps;
  hru.energy.canopy_sensible = store_canopy_sensible / (double) N_steps;
  hru.energy.deltaH = store_deltaH / (double) N_steps;
  hru.energy.fusion = store_fusion / (double) N_steps;
  hru.energy.grnd_flux = store_grnd_flux / (double) N_steps;
  hru.energy.latent = store_latent / (double) N_steps;
  hru.energy.latent_sub = store_latent_sub / (double) N_steps;
  hru.energy.melt_energy = store_melt_energy / (double) N_steps;
  hru.energy.sensible = store_sensible / (double) N_steps;
  hru.energy.glacier_flux = store_glacier_flux / (double) N_steps;
  hru.energy.deltaCC_glac = store_deltaCC_glac / (double) N_steps;
  hru.energy.advection = store_advection / (double) N_steps;
  hru.energy.deltaCC = store_deltaCC / (double) N_steps;
  hru.energy.refreeze_energy = store_refreeze_energy / (double) N_steps;
  hru.energy.snow_flux = store_snow_flux / (double) N_steps;
  hru.energy.Tfoliage = step_energy.Tfoliage;
  hru.energy.Tfoliage_fbflag = step_energy.Tfoliage_fbflag;
  hru.energy.Tfoliage_fbcount = step_energy.Tfoliage_fbcount;
  hru.energy.Tcanopy = Tcanopy;

  /**********************************************************
   Store vegetation variable sums for sub-model time steps
   **********************************************************/

  hru.veg_var[WET].throughfall = store_throughfall[WET];
  hru.veg_var[DRY].throughfall = store_throughfall[DRY];
  hru.veg_var[WET].canopyevap = store_canopyevap[WET];
  hru.veg_var[DRY].canopyevap = store_canopyevap[DRY];
  if (hru.snow.snow) {
    hru.veg_var[WET].Wdew = snow_veg_var[WET].Wdew;
    hru.veg_var[DRY].Wdew = snow_veg_var[DRY].Wdew;
  } else {
    hru.veg_var[WET].Wdew = soil_veg_var[WET].Wdew;
    hru.veg_var[DRY].Wdew = soil_veg_var[DRY].Wdew;
  }

  /**********************************************************
   Store soil layer variables for sub-model time steps
   **********************************************************/

  for (int lidx = 0; lidx < Nlayers; lidx++) {
    hru.cell[WET].layer[lidx] = step_layer[WET][lidx];
    hru.cell[DRY].layer[lidx] = step_layer[DRY][lidx];
    hru.cell[WET].layer[lidx].evap = store_layerevap[WET][lidx];
    hru.cell[DRY].layer[lidx].evap = store_layerevap[DRY][lidx];
#if EXCESS_ICE
    evap_prior[WET][lidx] = store_layerevap[WET][lidx];
    evap_prior[DRY][lidx] = store_layerevap[DRY][lidx];
#endif
  }
  if (store_aero_cond_used.surface > 0 && store_aero_cond_used.surface < HUGE_RESIST)
    hru.cell[WET].aero_resist.surface = 1 / (store_aero_cond_used.surface / (double) N_steps);
  else if (store_aero_cond_used.surface >= HUGE_RESIST)
    hru.cell[WET].aero_resist.surface = 0;
  else
    hru.cell[WET].aero_resist.surface = HUGE_RESIST;
  if (store_aero_cond_used.overstory > 0 && store_aero_cond_used.overstory < HUGE_RESIST)
    hru.cell[WET].aero_resist.overstory = 1 / (store_aero_cond_used.overstory / (double) N_steps);
  else if (store_aero_cond_used.overstory >= HUGE_RESIST)
    hru.cell[WET].aero_resist.overstory = 0;
  else
    hru.cell[WET].aero_resist.overstory = HUGE_RESIST;
  for (int p = 0; p < N_PET_TYPES; p++)
    hru.cell[WET].pot_evap[p] = store_pot_evap[p] / (double) N_steps;

  delete [] step_aero_resist;

  /********************************************************
   Compute Runoff, Baseflow, and Soil Moisture Transport
   ********************************************************/

  /* calculate glacier outflow and HRU runoff; all glacier outflow is assumed to be surface runoff */
  hru.glacier.inflow = ppt[WET] + ppt[DRY];
  ppt[WET] = 0.;
  ppt[DRY] = 0.;
  //  hru.cell[WET].runoff = 0.;
  //  hru.cell[DRY].runoff = 0.;
  //  hru.cell[WET].baseflow = 0.;
  //  hru.cell[DRY].baseflow = 0.;

  double wt = hru.glacier.water_storage + hru.glacier.inflow;
  double kt = soil_con->GLAC_KMIN + soil_con->GLAC_DK * exp(-soil_con->GLAC_A * hru.snow.swq * 1000.); /* multiply by 1000 to convert swe to millimetres */
  double qt = kt * wt;
  wt -= qt;   /* always positive as qt always less than wt */
  hru.glacier.water_storage = wt;
  hru.glacier.outflow = qt;
  hru.glacier.outflow_coef = kt;

#if EXCESS_ICE
  if(SubsidenceUpdate != 2) {
#endif
  hru.cell[WET].inflow = 0.;
  hru.cell[DRY].inflow = 0.;

  ErrorFlag = runoff(&hru.cell[WET], &hru.cell[DRY], &hru.energy, soil_con, ppt,
      SubsidenceUpdate, hru.mu, hru.bandIndex, rec, state);

  hru.cell[WET].runoff += (qt * 1000.); /* convert to mm; keeps units consistent with surface_fluxes.c */

  return (ErrorFlag);
#if EXCESS_ICE
} else {
	hru.cell[WET].runoff = qt * 1000.;
}
#endif

  return (0);
}

