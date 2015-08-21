#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

double solve_glacier(
      double BareAlbedo,
      double               Tgrnd,               // glacier slab temperature
      double               air_temp,            // air temperature
      double              *AlbedoUnder,
      double              *Le,
      double              *LongUnderIn,         // surface incoming LW
      double              *NetLongGlac,         // net LW at glacier surface
      double              *NetShortGlac,        // net SW at glacier surface
      double              *ShortUnderIn,        // surface incoming SW
      double              *Torg_snow,
      VegConditions       &aero_resist,
      AeroResistUsed      &aero_resist_used,
      VegConditions       &displacement,
      double              *melt_energy,
      double              *ppt,
      double              *rainfall,
      VegConditions       &ref_height,
      VegConditions       &roughness,
      VegConditions       &wind_speed,
      int                  dt,
      int                  rec,
      int                  hidx,
      VegConditions::VegSurfType &UnderStory,
      atmos_data_struct   *atmos,
      energy_bal_struct   *energy,
      glac_data_struct   *glacier,
      const soil_con_struct* soil,
      const ProgramState *state) {
/*********************************************************************

  This routine was written to handle the various calls and data
  handling needed to solve the various components of the new VIC
  glacier code for both the full_energy and water_balance models.

  Returns snow, veg_var, and energy variables for each elevation
  band.  Variable ppt[] is defined for elevation bands with snow.

*********************************************************************/

  int                 ErrorFlag;
  double              ShortOverIn;
  double              melt;
  double              tmp_grnd_flux;
  double              store_snowfall;
  double              tmp_ref_height;
  double              density;
  double              longwave;
  double              pressure;
  double              shortwave;
  double              vp;
  double              vpd;

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

  /* initialize storage for energy consumed in changing snowpack cover fraction */
  (*melt_energy)     = 0.;

  /** Compute latent heats **/
  (*Le) = (2.501e6 - 0.002361e6 * air_temp);

  /* initialize glacier surface radiation inputs */
  (*ShortUnderIn) = shortwave;
  (*LongUnderIn)  = longwave;

  /** compute net shortwave radiation **/
  (*AlbedoUnder) = BareAlbedo;
  (*NetShortGlac) = (1.0 - *AlbedoUnder) * (*ShortUnderIn);

   UnderStory = VegConditions::GLACIER_SURFACE_CASE;         /* glacier is present */

   /** Call glacier ablation algorithm **/
   ErrorFlag = glacier_melt((*Le), (*NetShortGlac), Tgrnd,
    roughness, aero_resist[UnderStory], aero_resist_used,
    air_temp, (double)dt * SECPHOUR, density,
    displacement[UnderStory],
    *LongUnderIn, pressure, rainfall[WET], vp, vpd,
    wind_speed[UnderStory], ref_height[UnderStory],
    NetLongGlac, Torg_snow, &melt, &energy->error,
    &energy->advection,
    &energy->deltaCC_glac, &energy->grnd_flux, &energy->latent,
    &energy->latent_sub,
    &energy->sensible,
    rec, glacier, soil, state);
   if ( ErrorFlag == ERROR ) return ( ERROR );

   // store melt water and rainfall
   ppt[WET] = (melt + rainfall[WET]/1000.); /* convert rainfall to m */

   // store glacier albedo
   energy->AlbedoUnder = *AlbedoUnder;

   rainfall[WET] = 0; /* all rain has been added to the glacier */

  return(melt);

}



