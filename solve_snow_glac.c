#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

double solve_snow_glac(
      double               BareAlbedo,
      double               Tgrnd, // soil surface temperature
      double               air_temp, // air temperature
      double               precipitation_mu,
      double              *AlbedoUnder,
      double              *latent_heat_Le,
      double              *LongUnderIn, // surface incoming LW
      double              *NetLongSnow, // net LW at snow surface
      double              *NetShortSnow, // net SW at snow surface
      double              *ShortUnderIn, // surface incoming SW
      double              *Torg_snow,
      VegConditions       &aero_resist,
      AeroResistUsed      &aero_resist_used,
      double              *coverage, // best guess snow coverage
      double              *delta_coverage, // cover fract change
      VegConditions       &displacement,
      double              *melt_energy,
      double              *ppt,
      double              *rainfall,
      VegConditions       &ref_height,
      VegConditions       &roughness,
      double              *snow_inflow,
      double              *snowfall,
      VegConditions       &wind_speed,
      int                  dt,
      int                  rec,
      int                  hidx,
      VegConditions::VegSurfType &UnderStory,
      const dmy_struct    *dmy,
      const atmos_data_struct &atmos,
      energy_bal_struct   *energy,
      snow_data_struct    *snow,
      const soil_con_struct *soil_con,
      glac_data_struct    *glacier,
      const ProgramState*  state) {
/*********************************************************************

  This routine was written to handle the various calls and data
  handling needed to solve the various components of the new VIC
  snow code for both the full_energy and water_balance models.

  This function is specific for glaciers and based off of the
  solve_snow function.

  Returns snow, veg_var, and energy variables for each elevation
  band.  Variable ppt[] is defined for elevation bands with snow.

*********************************************************************/

  char                ErrStr[MAXSTRING];
  char                FIRST_SOLN[1];
  int                 ErrorFlag;
  float               tempstep;
  double              ShortOverIn;
  double              melt;
  double              old_coverage;
  double              old_depth;
  double              old_swq;
  double              rainonly;
  double              tmp_Wdew[2];
  double              tmp_grnd_flux;
  int                 hour;
  int                 day_in_year;

  hour        = dmy[rec].hour;
  day_in_year = dmy[rec].day_in_year;


  /* initialize moisture variables */
  melt     = 0.;
  ppt[WET] = 0.;
  ppt[DRY] = 0.;

  /* initialize storage for energy consumed in changing snowpack cover fraction */
  (*melt_energy)     = 0.;

  /** Compute latent heats **/
  (*latent_heat_Le) = (2.501e6 - 0.002361e6 * air_temp);

  /** verify that distributed precipitation fraction equals 1 if
      snow is present or falling **/
    if ( precipitation_mu != 1 && state->options.FULL_ENERGY ) {
      fprintf(stderr,"ERROR: Snow model cannot be used if mu (%f) is not equal to 1.\n\tsolve_snow.c: record = %i",
        precipitation_mu, rec);
      return( ERROR );
    }
    else if ( precipitation_mu != 1 ) {
      fprintf(stderr,"WARNING: Snow is falling, but mu not equal to 1 (%f)\n",
        precipitation_mu);
      fprintf(stderr,"\trec = %i, hour = %i\n",rec,hour);
    }

  /* initialize understory radiation inputs */
  (*ShortUnderIn) = atmos.shortwave[hidx];
  (*LongUnderIn)  = atmos.longwave[hidx];

  snow->snow = TRUE; // snow is present during time step

  old_coverage = snow->coverage; // store previous coverage fraction

  energy->NetLongOver = 0;
  energy->LongOverIn = 0;

  (*snow_inflow) = rainfall[WET] + snowfall[WET];

  old_swq = snow->swq; /* store swq for density calculations */
  UnderStory = VegConditions::SNOW_COVERED_CASE; /* ground snow is present of accumulating during time step */

#if SPATIAL_SNOW
  /* make snowpack uniform at mean depth */
  if ( snowfall[WET] > 0 ) snow->coverage = 1;
  if (snow->coverage > 0 && snowfall[WET] == 0) {
    if ( snow->coverage < 1) {
      /* rain falls evenly over grid cell */
      ppt[WET] = rainfall[WET] * (1.0 - snow->coverage);
      rainfall[WET] *= snow->coverage;
    }
  }
#endif

  /** compute understory albedo and net shortwave radiation **/
  if (snow->swq > 0. && snowfall[WET] == 0.) {
    // age snow albedo if no new snowfall
    // ignore effects of snow dropping from canopy; only consider fresh snow from sky
    snow->last_snow++;
    snow->albedo = snow_albedo(snowfall[WET], snow->swq, snow->depth,
        snow->albedo, snow->coldcontent, (double) dt, snow->last_snow,
        snow->MELTING, soil_con, state);
    (*AlbedoUnder) = (*coverage * snow->albedo + (1. - *coverage) * BareAlbedo);
  } else {
    // set snow albedo to new snow albedo
    snow->last_snow = 0;
    snow->albedo = soil_con->NEW_SNOW_ALB;
    (*AlbedoUnder) = snow->albedo;
  }
  (*NetShortSnow) = (1.0 - *AlbedoUnder) * (*ShortUnderIn);

  /** Call snow pack accumulation and ablation algorithm **/
  ErrorFlag = snow_melt_glac((*latent_heat_Le), (*NetShortSnow),
      Tgrnd, roughness, aero_resist[UnderStory], aero_resist_used,
      air_temp, *coverage, (double) dt * SECPHOUR, atmos.density[hidx],
      displacement[UnderStory], *LongUnderIn,
      atmos.pressure[hidx], rainfall[WET], snowfall[WET], atmos.vp[hidx],
      atmos.vpd[hidx], wind_speed[UnderStory], ref_height[UnderStory], NetLongSnow,
      Torg_snow, &melt, &energy->error, &energy->advected_sensible,
      &energy->advection, &energy->deltaCC, &energy->grnd_flux, &energy->latent,
      &energy->latent_sub, &energy->refreeze_energy, &energy->sensible, rec,
      snow, soil_con, glacier, state);

  if (ErrorFlag == ERROR) return (ERROR);

  // store melt water
  ppt[WET] += melt;

  // store snow albedo
  energy->AlbedoUnder = *AlbedoUnder;

  /** Compute Snow Parameters **/
  if (snow->swq > 0.) {

    /** Calculate Snow Density **/
    if (IS_VALID(snow->surf_temp) && snow->surf_temp <= 0)
      // snowpack present, compress and age density
      snow->density = snow_density(snow, snowfall[WET], old_swq, Tgrnd,
          air_temp, (double) dt, state);
    else
    // no snowpack present, start with new snow density
    if (snow->last_snow == 0)
      snow->density = new_snow_density(air_temp, state);

    /** Calculate Snow Depth (H.B.H. 7.2.1) **/
    old_depth = snow->depth;
    snow->depth = 1000. * snow->swq / snow->density;

    /** Record if snowpack is melting this time step **/
    if (snow->coldcontent >= 0 && ((soil_con->lat >= 0 && (day_in_year > 60 // ~ March 1
    && day_in_year < 273)) // ~ October 1
    || (soil_con->lat < 0 && (day_in_year < 60 // ~ March 1
    || day_in_year > 273)) // ~ October 1
    ))
      snow->MELTING = TRUE;
    else if (snow->MELTING && snowfall[WET] > TraceSnow)
      snow->MELTING = FALSE;

    /** Check for Thin Snowpack which only Partially Covers Grid Cell
     exists only if not snowing and snowpack has started to melt **/
#if SPATIAL_SNOW
  snow->coverage = calc_snow_coverage(&snow->store_snow,
              soil_con->depth_full_snow_cover,
              old_coverage, snow->swq,
              old_swq, snow->depth, old_depth,
              melt*0.001 + snow->vapor_flux,
              &snow->max_swq, snowfall,
              &snow->store_swq,
              &snow->swq_slope,
              &snow->store_coverage);

#else

    if (snow->swq > 0)
      snow->coverage = 1.;
    else
      snow->coverage = 0.;
#endif

  } else {
    snow->coverage = 0.;
  }

  *delta_coverage = old_coverage - snow->coverage;

  if (*delta_coverage != 0) {

    /* returns mixed surface albedo if snow cover fraction has
     decreased (old_coverage is cover fraction for previous
     time step, snow->coverage is cover fraction for current
     time step. */
    if (old_coverage > snow->coverage) {
      /* melt has occured */
      *coverage = (old_coverage);
      (*AlbedoUnder) = (*coverage - snow->coverage) / (1. - snow->coverage)
          * snow->albedo;
      (*AlbedoUnder) += (1. - *coverage) / (1. - snow->coverage) * BareAlbedo;

      /* compute snowpack energy used in reducing coverage area */
      (*melt_energy) = (*delta_coverage)
          * (energy->advection - energy->deltaCC + energy->latent
              + energy->latent_sub + energy->sensible + energy->refreeze_energy
              + energy->advected_sensible);
    } else if (old_coverage < snow->coverage) {
#if VERBOSE
      if (snow->coverage != 1.)
        fprintf(stderr,
            "WARNING: snow cover fraction has increased, but it is not equal to 1 (%f).\n",
            snow->coverage);
#endif // VERBOSE
      *coverage = snow->coverage;
      *delta_coverage = 0;
    } else {
      *coverage = snow->coverage;
      *delta_coverage = 0.;
    }
  } else if (old_coverage == 0 && snow->coverage == 0) {
    // snow falls and melts all in one time step
    *delta_coverage = 1.;
    *coverage = 0.;
    (*melt_energy) = (energy->advection - energy->deltaCC + energy->latent
        + energy->latent_sub + energy->sensible + energy->refreeze_energy
        + energy->advected_sensible);
  }

  /** Compute energy balance components for snowpack */

  (*NetLongSnow) *= (snow->coverage);
  (*NetShortSnow) *= (snow->coverage);
  energy->latent *= (snow->coverage + *delta_coverage);
  energy->latent_sub *= (snow->coverage + *delta_coverage);
  energy->sensible *= (snow->coverage + *delta_coverage);

  if (snow->swq == 0) {

    /** Reset Snow Pack Variables after Complete Melt **/

    /*** NOTE *coverage should not be zero the time step the
     snowpack melts - FIX THIS ***/

    snow->density = 0.;
    snow->depth = 0.;
    snow->surf_water = 0;
    snow->pack_water = 0;
    snow->surf_temp = 0;
    snow->pack_temp = 0;
    snow->coverage = 0;
    snow->swq_slope = 0;
    snow->store_snow = TRUE;
    snow->MELTING = FALSE;

  }

  snowfall[WET] = 0; /* all falling snow has been added to the pack */
  rainfall[WET] = 0; /* all rain has been added to the pack */

  energy->melt_energy *= -1.;

  return (melt);
}
