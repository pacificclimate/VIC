#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

double solve_glacier(double               LongUnderOut, // LW from understory
      /* double               MIN_RAIN_TEMP, */
      /* double               MAX_SNOW_TEMP, */
      /* double               Tcanopy, // canopy air temperature */
      double               Tgrnd, // glacier slab temperature
      double               air_temp, // air temperature
      double               mu,
      double               prec,
      double               snow_grnd_flux,
      double               wind_h,
      double              *AlbedoUnder,
      double              *Evap,
      double              *Le,
      double              *LongUnderIn, // surface incoming LW
      double              *NetLongSnow, // net LW at glacier surface
      double              *NetShortGrnd, // net SW reaching ground
      double              *NetShortSnow, // net SW at glaciersurface
      double              *ShortUnderIn, // surfave incoming SW
      double              *Torg_snow,
      double              *aero_resist,
      double              *aero_resist_used,
      /* double              *coverage, // best guess snow coverage */
      /* double              *delta_coverage, // cover fract change */
      /* double              *delta_snow_heat, // change in pack heat */
      double              *displacement,
      double              *gauge_correction,
      double              *melt_energy,
      double              *out_prec,
      double              *out_rain,
      double              *out_snow,
      double              *ppt,
      double              *rainfall,
      double              *ref_height,
      double              *roughness,
      /* double              *snow_inflow, */
      double              *snowfall,
      /* double              *surf_atten, */
      double              *wind,
      /* float               *root, */
      /* int                  INCLUDE_SNOW, */
      int                  Nveg,
      int                  iveg,
      int                  band,
      int                  dt,
      int                  rec,
      int                  hidx,
      int                 *UnderStory,
      dmy_struct          *dmy,
      atmos_data_struct   *atmos,
      energy_bal_struct   *energy,
      /* snow_data_struct    *snow, */
      glac_data_struct   *glacier) {
/*********************************************************************
  solve_snow.c                Keith Cherkauer       July 2, 1998

  This routine was written to handle the various calls and data
  handling needed to solve the various components of the new VIC
  glacier code for both the full_energy and water_balance models.

  Returns snow, veg_var, and energy variables for each elevation
  band.  Variable ppt[] is defined for elevation bands with snow.

*********************************************************************/

  extern option_struct   options;
  extern veg_lib_struct *veg_lib;

  int                 ErrorFlag;
  double              ShortOverIn;
  double              melt;
  /* double              old_coverage; */
  /* double              old_depth; */
  /* double              old_swq; */
  /* double              rainonly; */
  /* double              tmp_Wdew[2]; */
  double              tmp_grnd_flux;
  double              store_snowfall;
  double              tmp_ref_height;
  /* int                 month; */
  /* int                 hour; */
  double              density;
  double              longwave;
  double              pressure;
  double              shortwave;
  double              vp;
  double              vpd;

  /* month       = dmy[rec].month; */
  /* hour        = dmy[rec].hour; */

  density   = atmos->density[hidx];
  longwave  = atmos->longwave[hidx];
  pressure  = atmos->pressure[hidx];
  shortwave = atmos->shortwave[hidx];
  vp        = atmos->vp[hidx];
  vpd       = atmos->vpd[hidx];

  /* initialize moisture variables */
  melt     = 0.;
  ppt[WET] = 0.;
  ppt[DRY] = 0.;

  /* initialize storage for energy consumed in changing snowpack
     cover fraction */
  (*melt_energy)     = 0.;

  /** Compute latent heats **/
  (*Le) = (2.501e6 - 0.002361e6 * air_temp);

  /* initialize glacier surface radiation inputs */
  (*ShortUnderIn) = shortwave;
  (*LongUnderIn)  = longwave;

   energy->NetLongOver = 0;
   energy->LongOverIn  = 0;

   /* (*NetShortGrnd) = 0.; */

   /* (*snow_inflow) += rainfall[WET] + snowfall[WET]; */

   (*UnderStory) = 2;         /* ground snow is present or accumulating during time step */

   /** compute net shortwave radiation **/
   (*AlbedoUnder) = GLAC_ALBEDO;
   (*NetShortSnow) = (1.0 - *AlbedoUnder) * (*ShortUnderIn);

   /** Call snow pack accumulation and ablation algorithm **/
   ErrorFlag = glacier_melt((*Le), (*NetShortSnow), Tgrnd,
    roughness, aero_resist[*UnderStory], aero_resist_used,
    air_temp, (double)dt * SECPHOUR, density,
    displacement[*UnderStory], snow_grnd_flux,
    *LongUnderIn, pressure, rainfall[WET], vp, vpd,
    wind[*UnderStory], ref_height[*UnderStory],
    NetLongSnow, Torg_snow, &melt, &energy->error,
    &energy->advected_sensible, &energy->advection,
    &energy->deltaCC, &tmp_grnd_flux, &energy->latent,
    &energy->latent_sub,
    &energy->sensible,
    rec, iveg, band, glacier);
   if ( ErrorFlag == ERROR ) return ( ERROR );

   // store melt water and rainfall
   ppt[WET] += (melt + rainfall[WET]);

   // store glacier albedo
   energy->AlbedoUnder = *AlbedoUnder;

   rainfall[WET] = 0; /* all rain has been added to the glacier */

  energy->melt_energy = 0.;

  return(melt);

}



