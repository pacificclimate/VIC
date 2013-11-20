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
       double               current_prcp_mu,
       double               surf_atten,
       double              *Melt,
       double              *latent_heat_Le,
       double             **aero_resist,
       double              *displacement,
       double              *gauge_correction,
       double              *out_prec,
       double              *out_rain,
       double              *out_snow,
       double              *ref_height,
       double              *roughness,
       double              *snow_inflow,
       double              *wind,
       const float         *root,
       int                  Nbands,
       int                  Ndist,
       int                  Nlayers,
       int                  Nveg,
       int                  band,
       int                  dp,
       int                  iveg,
       int                  rec,
       int                  veg_class,
       atmos_data_struct   *atmos,
       const dmy_struct    *dmy,
       energy_bal_struct   *energy,
       hru_data_struct    *cell_dry,
       hru_data_struct    *cell_wet,
       snow_data_struct    *snow,
       const soil_con_struct *soil_con,
       veg_var_struct      *veg_var_dry,
       veg_var_struct      *veg_var_wet,
       glac_data_struct   *glacier,
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
  int                    UnderStory;
  int                    dist;
  int                    hidx;     // index of initial element of atmos array
  int                    step_inc; // number of atmos array elements to skip per surface fluxes step
  int                    endhidx;  // index of final element of atmos array
  int                    step_dt;  // time length of surface fluxes step
  int                    lidx;
  int                    p,q;
  double                 Ls;
  double                 LongUnderIn; // inmoing LW to ground surface
  double                 LongUnderOut; // outgoing LW from ground surface
  double                 NetLongSnow; // net LW over snowpack
  double                 NetShortSnow; // net SW over understory
  double                 NetShortGrnd; // net SW over snow-free surface
  double                 OldTSurf; // previous snow surface temperature
  double                 ShortUnderIn; // incoming SW to understory
  double                 Tair; // air temperature
  double                 Tcanopy; // canopy air temperature
  double                 Tgrnd; // soil surface temperature
  double                 Tsurf; // ground surface temperature
  double                 VPDcanopy; // vapor pressure deficit in canopy/atmos
  double                 coverage; // mid-step snow cover fraction
  double                 delta_coverage; // change in snow cover fraction
  double                 delta_snow_heat; // change in snowpack heat storage
  double                 ppt[2]; // precipitation/melt reaching soil surface or glacier surface
  double                 rainfall[2]; // rainfall
  double                 snowfall[2]; // snowfall
  double                 snow_flux; // heat flux through snowpack

  // Step-specific quantities
  double                 step_Wdew[2];
  double                 step_melt;
  double                 step_melt_energy;  /* energy used to reduce snow coverage */
  double                 step_out_prec;
  double                 step_out_rain;
  double                 step_out_snow;
  double                 step_ppt[2];
  double                 step_prec[2];
  double               **step_aero_resist;
  double step_melt_glac;

  // Quantities that need to be summed or averaged over multiple snow steps
  // energy structure
  double                 store_AlbedoOver;
  double                 store_AlbedoUnder;
  double                 store_AtmosLatent;
  double                 store_AtmosLatentSub;
  double                 store_AtmosSensible;
  double                 store_LongOverIn;
  double                 store_LongUnderIn;
  double                 store_LongUnderOut;
  double                 store_NetLongAtmos;
  double                 store_NetLongOver;
  double                 store_NetLongUnder;
  double                 store_NetShortAtmos;
  double                 store_NetShortGrnd;
  double                 store_NetShortOver;
  double                 store_NetShortUnder;
  double                 store_ShortOverIn;
  double                 store_ShortUnderIn;
  double                 store_advected_sensible;
  double                 store_advection;
  double                 store_canopy_advection;
  double                 store_canopy_latent;
  double                 store_canopy_latent_sub;
  double                 store_canopy_sensible;
  double                 store_canopy_refreeze;
  double                 store_deltaCC;
  double                 store_deltaH;
  double                 store_fusion;
  double                 store_grnd_flux;
  double                 store_latent;
  double                 store_latent_sub;
  double                 store_melt_energy;
  double                 store_refreeze_energy;
  double                 store_sensible;
  double                 store_snow_flux;
  double                 store_deltaCC_glac;
  double                 store_glacier_flux;
  // glacier structure
  double                 store_melt_glac;
  double                 store_vapor_flux_glac;
  double                 store_mass_balance;
  double                 store_ice_mass_balance;
  double                 store_accum_glac;
  double                 store_inflow_glac;
  // snow structure
  double                 store_canopy_vapor_flux;
  double                 store_melt;
  double                 store_vapor_flux;
  double                 store_blowing_flux;
  double                 store_surface_flux;
  // veg_var structure
  double                 store_canopyevap[2];
  double                 store_throughfall[2];
  // cell structure
  double                 store_layerevap[2][MAX_LAYERS];
  double                 store_ppt[2];
  double                 store_aero_cond_used[2];
  double                 store_pot_evap[N_PET_TYPES];

  // Structures holding values for current snow step
  energy_bal_struct      step_energy; // energy fluxes at snowpack surface and glacier surface
  veg_var_struct         snow_veg_var[2]; // veg fluxes/storages in presence of snow
  veg_var_struct         soil_veg_var[2]; // veg fluxes/storages in soil energy balance
  snow_data_struct       step_snow;
  layer_data_struct      step_layer[2][MAX_LAYERS];
  glac_data_struct       step_glacier;

  // Structures holding values for current iteration
  double                 temp_aero_resist[3];
  double                 temp_aero_resist_used[2];
  double                 stability_factor[2];
  double                 step_pot_evap[N_PET_TYPES];

  step_aero_resist = (double**) calloc(N_PET_TYPES, sizeof(double*));
  for (p = 0; p < N_PET_TYPES; p++) {
    step_aero_resist[p] = (double*) calloc(2, sizeof(double));
  }

  /***********************************************************************
   Set temporary variables - preserves original values until iterations
   are completed
   ***********************************************************************/

  energy->advection = 0;
  energy->deltaCC = 0;
  if (snow->swq > 0) {
    snow_flux = energy->snow_flux;
  } else
    snow_flux = 0;
  energy->refreeze_energy = 0;
  coverage = snow->coverage;
  step_energy = (*energy);
  snow_veg_var[WET] = (*veg_var_wet);
  snow_veg_var[DRY] = (*veg_var_dry);
  soil_veg_var[WET] = (*veg_var_wet);
  soil_veg_var[DRY] = (*veg_var_dry);
  step_snow = (*snow);
  step_glacier = (*glacier);
  for (lidx = 0; lidx < Nlayers; lidx++) {
    step_layer[WET][lidx] = cell_wet->layer[lidx];
    step_layer[DRY][lidx] = cell_dry->layer[lidx];
  }
  for (lidx = 0; lidx < Nlayers; lidx++) {
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

  // energy structure
  store_AlbedoOver = 0;
  store_AlbedoUnder = 0;
  store_AtmosLatent = 0;
  store_AtmosLatentSub = 0;
  store_AtmosSensible = 0;
  store_LongOverIn = 0;
  store_LongUnderIn = 0;
  store_LongUnderOut = 0;
  store_NetLongAtmos = 0;
  store_NetLongOver = 0;
  store_NetLongUnder = 0;
  store_NetShortAtmos = 0;
  store_NetShortGrnd = 0;
  store_NetShortOver = 0;
  store_NetShortUnder = 0;
  store_ShortOverIn = 0;
  store_ShortUnderIn = 0;
  store_advected_sensible = 0;
  store_advection = 0;
  store_canopy_advection = 0;
  store_canopy_latent = 0;
  store_canopy_latent_sub = 0;
  store_canopy_sensible = 0;
  store_canopy_refreeze = 0;
  store_deltaCC = 0;
  store_deltaH = 0;
  store_fusion = 0;
  store_grnd_flux = 0;
  store_latent = 0;
  store_latent_sub = 0;
  store_melt_energy = 0;
  store_refreeze_energy = 0;
  store_sensible = 0;
  store_snow_flux = 0;
  store_deltaCC_glac = 0;
  store_glacier_flux = 0;
  // glacier structure
  store_melt_glac = 0;
  store_vapor_flux_glac = 0;
  store_mass_balance = 0;
  store_ice_mass_balance = 0;
  store_accum_glac = 0;
  store_inflow_glac = 0;
  // snow structure
  store_canopy_vapor_flux = 0;
  store_melt = 0;
  store_vapor_flux = 0;
  store_surface_flux = 0;
  store_blowing_flux = 0;
  // veg_var and cell structures
  for (dist = 0; dist < Ndist; dist++) {
    store_throughfall[dist] = 0.;
    store_canopyevap[dist] = 0.;
    for (lidx = 0; lidx < state->options.Nlayer; lidx++) {
      store_layerevap[dist][lidx] = 0.;
    }
  }
  step_Wdew[WET] = veg_var_wet->Wdew;
  step_Wdew[DRY] = veg_var_wet->Wdew;
  // misc
  store_ppt[WET] = 0;
  store_ppt[DRY] = 0;
  step_prec[DRY] = 0;
  store_aero_cond_used[0] = 0;
  store_aero_cond_used[1] = 0;
  (*snow_inflow) = 0;
  for (p = 0; p < N_PET_TYPES; p++)
    store_pot_evap[p] = 0;
  N_steps = 0;

  /*************************
   Compute surface fluxes
   *************************/

  do { /* MPN TODO This is stupid.  Makes it into a for loop because it is a fixed number of iterations! */

    /** Solve energy balance for all sub-model time steps **/

    /* set air temperature and precipitation for this snow band */
    Tair = atmos->air_temp[hidx] + soil_con->Tfactor[band];
    step_prec[WET] = atmos->prec[hidx] / current_prcp_mu * soil_con->Pfactor[band];

    // initialize ground surface temperaure
    Tgrnd = GLAC_TEMP;

    // initialize canopy terms
    Tcanopy = 0.;
    VPDcanopy = atmos->vpd[hidx];

    // Compute mass flux of blowing snow
    if (state->options.BLOWING && step_snow.swq > 0.) {
      Ls = (677. - 0.07 * step_snow.surf_temp) * JOULESPCAL * GRAMSPKG;
      step_snow.blowing_flux = CalcBlowingSnow((double) step_dt, Tair,
          step_snow.last_snow, step_snow.surf_water, wind[2], Ls,
          atmos->density[hidx], atmos->pressure[hidx], atmos->vp[hidx],
          roughness[2], ref_height[2], step_snow.depth, lag_one, sigma_slope,
          step_snow.surf_temp, iveg, Nveg, fetch, displacement[1], roughness[1],
          &step_snow.transport);
      if ((int) step_snow.blowing_flux == ERROR) {
        return (ERROR);
      }
      step_snow.blowing_flux *= step_dt * SECPHOUR / RHO_W; /* m/time step */
    } else
      step_snow.blowing_flux = 0.0;



        /** Iterate for understory solution - itererates to find snow flux **/

        // Initialize structures for new iteration

        for (q = 0; q < 3; q++) {
          temp_aero_resist[q] = aero_resist[N_PET_TYPES][q];
        }
        temp_aero_resist_used[0] = cell_wet->aero_resist[0];
        temp_aero_resist_used[1] = cell_wet->aero_resist[1];
        step_snow.canopy_vapor_flux = 0;
        step_snow.vapor_flux = 0;
        step_snow.surface_flux = 0;
        /* iter_snow.blowing_flux has already been reset to step_snow.blowing_flux */
        LongUnderOut = step_energy.LongUnderOut;

        if (step_snow.swq > 0) {
          /** Solve snow accumulation, ablation and interception **/


          step_melt = 0; // TODO: step_melt_glac = ... ?
          // TODO: check these arguments (especially -snow_flux?)
          /*step_melt = solve_snow(LongUnderOut,
              soil_con->MIN_RAIN_TEMP, soil_con->MAX_SNOW_TEMP, Tgrnd, Tair, dp, current_prcp_mu,
              step_prec[WET], (-snow_flux), state->global_param.wind_h, &energy->AlbedoUnder,
              latent_heat_Le, &LongUnderIn, &NetLongSnow, &NetShortGrnd,
              &NetShortSnow, &ShortUnderIn, &OldTSurf, temp_aero_resist,
              temp_aero_resist_used, &coverage, &delta_coverage, &delta_snow_heat,
              displacement, gauge_correction, &step_melt_energy, &step_out_prec,
              &step_out_rain, &step_out_snow, step_ppt, rainfall, ref_height,
              roughness, snow_inflow, snowfall, &surf_atten, wind, root,
              state->options.Nnode, Nveg, iveg, band, step_dt, rec, hidx,
              veg_class, &UnderStory, dmy, *atmos, &(step_energy),
              &(step_snow), soil_con, state);
          */
          if (step_melt == ERROR)
            return (ERROR);
        } else {

          step_melt_glac = 0;//TODO: step_melt_glac = solve_glac(...);


          if (step_melt_glac == ERROR) {
            return (ERROR);
          }
        }


    /**************************************
     Compute Potential Evap
     **************************************/
    // First, determine the stability correction used in the iteration
    if (temp_aero_resist_used[0] == HUGE_RESIST)
      stability_factor[0] = HUGE_RESIST;
    else
      stability_factor[0] = temp_aero_resist_used[0]
          / aero_resist[N_PET_TYPES][UnderStory];
    if (temp_aero_resist_used[1] == temp_aero_resist_used[0])
      stability_factor[1] = stability_factor[0];
    else {
      if (temp_aero_resist_used[1] == HUGE_RESIST)
        stability_factor[1] = HUGE_RESIST;
      else
        stability_factor[1] = temp_aero_resist_used[1]
            / aero_resist[N_PET_TYPES][1];
    }

    // Next, loop over pot_evap types and apply the correction to the relevant aerodynamic resistance
    for (p = 0; p < N_PET_TYPES; p++) {
      if (stability_factor[0] == HUGE_RESIST)
        step_aero_resist[p][0] = HUGE_RESIST;
      else
        step_aero_resist[p][0] = aero_resist[p][UnderStory]
            * stability_factor[0];
      if (stability_factor[1] == HUGE_RESIST)
        step_aero_resist[p][1] = HUGE_RESIST;
      else
        step_aero_resist[p][1] = aero_resist[p][1] * stability_factor[1];
    }

    // Finally, compute pot_evap
    compute_pot_evap(veg_class, dmy, rec, state->global_param.dt, atmos->shortwave[hidx],
        step_energy.NetLongAtmos, Tair, VPDcanopy, soil_con->elevation,
        step_aero_resist, step_pot_evap, state);

    /**************************************
     Store sub-model time step variables
     **************************************/

    for (dist = 0; dist < Ndist; dist++) {
      store_ppt[dist] += step_ppt[dist];
    }
    if (temp_aero_resist_used[0] > 0)
      store_aero_cond_used[0] += 1 / temp_aero_resist_used[0];
    else
      store_aero_cond_used[0] += HUGE_RESIST;
    if (temp_aero_resist_used[1] > 0)
      store_aero_cond_used[1] += 1 / temp_aero_resist_used[1];
    else
      store_aero_cond_used[1] += HUGE_RESIST;

    store_melt += step_melt;
    store_vapor_flux += step_snow.vapor_flux;
    store_surface_flux += step_snow.surface_flux;
    store_blowing_flux += step_snow.blowing_flux;

    out_prec[0] += step_out_prec * current_prcp_mu;
    out_rain[0] += step_out_rain * current_prcp_mu;
    out_snow[0] += step_out_snow * current_prcp_mu;

    store_AlbedoUnder += step_energy.AlbedoUnder;
    store_AtmosLatent += step_energy.AtmosLatent;
    store_AtmosLatentSub += step_energy.AtmosLatentSub;
    store_AtmosSensible += step_energy.AtmosSensible;
    store_LongUnderIn += LongUnderIn;
    store_LongUnderOut += step_energy.LongUnderOut;
    store_NetLongAtmos += step_energy.NetLongAtmos;
    store_NetLongUnder += step_energy.NetLongUnder;
    store_NetShortAtmos += step_energy.NetShortAtmos;
    store_NetShortUnder += step_energy.NetShortUnder;
    store_ShortUnderIn += step_energy.ShortUnderIn;
    store_latent += step_energy.latent;
    store_latent_sub += step_energy.latent_sub;
    store_melt_energy += step_melt_energy;
    store_sensible += step_energy.sensible;
    // glacier
    store_melt_glac += step_melt_glac;
    store_vapor_flux_glac += step_glacier.vapor_flux;
    store_mass_balance += step_glacier.mass_balance;
    store_ice_mass_balance += step_glacier.ice_mass_balance;
    store_accum_glac += step_glacier.accumulation;
    store_inflow_glac += step_glacier.inflow;
    store_glacier_flux += step_energy.glacier_flux;
    store_deltaCC_glac += step_energy.deltaCC_glac;

    if (step_snow.snow) {
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
    }
    for (p = 0; p < N_PET_TYPES; p++)
      store_pot_evap[p] += step_pot_evap[p];

    /* increment time step */
    N_steps++;
    hidx += step_inc;

  } while (hidx < endhidx);

  /************************************************
   Store glacier variables for sub-model time steps
   ************************************************/
  (*glacier) = step_glacier;
  glacier->melt = store_melt_glac;
  glacier->vapor_flux = store_vapor_flux_glac;
  glacier->mass_balance = store_mass_balance;
  glacier->ice_mass_balance = store_ice_mass_balance;
  glacier->accumulation = store_accum_glac;
  glacier->inflow = store_inflow_glac;

  /************************************************
   Store snow variables for sub-model time steps
   ************************************************/

  (*snow) = step_snow;
  snow->vapor_flux = store_vapor_flux;
  snow->blowing_flux = store_blowing_flux;
  snow->surface_flux = store_surface_flux;
  snow->canopy_vapor_flux = store_canopy_vapor_flux;
  (*Melt) = store_melt + store_melt_glac;
  snow->melt = store_melt;
  for (dist = 0; dist < 2; dist++)
    ppt[dist] = store_ppt[dist];

  /******************************************************
   Store energy flux averages for sub-model time steps
   ******************************************************/

  (*energy) = step_energy;
  energy->AlbedoOver = store_AlbedoOver / (double) N_steps;
  energy->AlbedoUnder = store_AlbedoUnder / (double) N_steps;
  energy->AtmosLatent = store_AtmosLatent / (double) N_steps;
  energy->AtmosLatentSub = store_AtmosLatentSub / (double) N_steps;
  energy->AtmosSensible = store_AtmosSensible / (double) N_steps;
  energy->LongOverIn = store_LongOverIn / (double) N_steps;
  energy->LongUnderIn = store_LongUnderIn / (double) N_steps;
  energy->LongUnderOut = store_LongUnderOut / (double) N_steps;
  energy->NetLongAtmos = store_NetLongAtmos / (double) N_steps;
  energy->NetLongOver = store_NetLongOver / (double) N_steps;
  energy->NetLongUnder = store_NetLongUnder / (double) N_steps;
  energy->NetShortAtmos = store_NetShortAtmos / (double) N_steps;
  energy->NetShortGrnd = store_NetShortGrnd / (double) N_steps;
  energy->NetShortOver = store_NetShortOver / (double) N_steps;
  energy->NetShortUnder = store_NetShortUnder / (double) N_steps;
  energy->ShortOverIn = store_ShortOverIn / (double) N_steps;
  energy->ShortUnderIn = store_ShortUnderIn / (double) N_steps;
  energy->advected_sensible = store_advected_sensible / (double) N_steps;
  energy->canopy_advection = store_canopy_advection / (double) N_steps;
  energy->canopy_latent = store_canopy_latent / (double) N_steps;
  energy->canopy_latent_sub = store_canopy_latent_sub / (double) N_steps;
  energy->canopy_refreeze = store_canopy_refreeze / (double) N_steps;
  energy->canopy_sensible = store_canopy_sensible / (double) N_steps;
  energy->deltaH = store_deltaH / (double) N_steps;
  energy->fusion = store_fusion / (double) N_steps;
  energy->grnd_flux = store_grnd_flux / (double) N_steps;
  energy->latent = store_latent / (double) N_steps;
  energy->latent_sub = store_latent_sub / (double) N_steps;
  energy->melt_energy = store_melt_energy / (double) N_steps;
  energy->sensible = store_sensible / (double) N_steps;
  energy->glacier_flux = store_glacier_flux / (double) N_steps;
  energy->deltaCC_glac = store_deltaCC_glac / (double) N_steps;
  if (snow->snow) {
    energy->advection = store_advection / (double) N_steps;
    energy->deltaCC = store_deltaCC / (double) N_steps;
    energy->refreeze_energy = store_refreeze_energy / (double) N_steps;
    energy->snow_flux = store_snow_flux / (double) N_steps;
  }
  energy->Tfoliage = step_energy.Tfoliage;
  energy->Tfoliage_fbflag = step_energy.Tfoliage_fbflag;
  energy->Tfoliage_fbcount = step_energy.Tfoliage_fbcount;
  energy->Tcanopy = Tcanopy;

  /**********************************************************
   Store vegetation variable sums for sub-model time steps
   **********************************************************/

    veg_var_wet->throughfall = store_throughfall[WET];
    veg_var_dry->throughfall = store_throughfall[DRY];
    veg_var_wet->canopyevap = store_canopyevap[WET];
    veg_var_dry->canopyevap = store_canopyevap[DRY];
    if (snow->snow) {
      veg_var_wet->Wdew = snow_veg_var[WET].Wdew;
      veg_var_dry->Wdew = snow_veg_var[DRY].Wdew;
    } else {
      veg_var_wet->Wdew = soil_veg_var[WET].Wdew;
      veg_var_dry->Wdew = soil_veg_var[DRY].Wdew;
    }

  /**********************************************************
   Store soil layer variables for sub-model time steps
   **********************************************************/

  for (lidx = 0; lidx < Nlayers; lidx++) {
    cell_wet->layer[lidx] = step_layer[WET][lidx];
    cell_dry->layer[lidx] = step_layer[DRY][lidx];
    cell_wet->layer[lidx].evap = store_layerevap[WET][lidx];
    cell_dry->layer[lidx].evap = store_layerevap[DRY][lidx];
#if EXCESS_ICE
    evap_prior_wet[lidx] = store_layerevap[WET][lidx];
    evap_prior_dry[lidx] = store_layerevap[DRY][lidx];
#endif
  }
  if (store_aero_cond_used[0] > 0 && store_aero_cond_used[0] < HUGE_RESIST)
    cell_wet->aero_resist[0] = 1 / (store_aero_cond_used[0] / (double) N_steps);
  else if (store_aero_cond_used[0] >= HUGE_RESIST)
    cell_wet->aero_resist[0] = 0;
  else
    cell_wet->aero_resist[0] = HUGE_RESIST;
  if (store_aero_cond_used[1] > 0 && store_aero_cond_used[1] < HUGE_RESIST)
    cell_wet->aero_resist[1] = 1 / (store_aero_cond_used[1] / (double) N_steps);
  else if (store_aero_cond_used[1] >= HUGE_RESIST)
    cell_wet->aero_resist[1] = 0;
  else
    cell_wet->aero_resist[1] = HUGE_RESIST;
  for (p = 0; p < N_PET_TYPES; p++)
    cell_wet->pot_evap[p] = store_pot_evap[p] / (double) N_steps;

  for (p = 0; p < N_PET_TYPES; p++) {
    free((char *) step_aero_resist[p]);
  }
  free((char *) step_aero_resist);

  /********************************************************
   Compute Runoff, Baseflow, and Soil Moisture Transport
   ********************************************************/

#if EXCESS_ICE
  if(SubsidenceUpdate != 2) {
#endif
  cell_wet->inflow = 0.;
  cell_dry->inflow = 0.;
  glacier->water_storage += glacier->inflow;

  // TODO: implement and call runoff_glac with the correct parameters.
  //ErrorFlag = runoff_glac(cell_wet, cell_dry, energy, soil_con, ppt,
  //   SubsidenceUpdate, current_prcp_mu, band, rec, iveg, state);

  return (ErrorFlag);
#if EXCESS_ICE
}
#endif

  return (0);
}

