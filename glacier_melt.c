#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "GlacierEnergyBalance.h"
#include "vicNl.h"
// Forward declaration for the error function which just prints out all the variables.
double ErrorPrintGlacierEnergyBalance(double TSurf, int rec,
    double Dt, double Ra, AeroResistUsed& Ra_used, double Displacement, double Z,
    VegConditions& roughness, double AirDens, double EactAir, double LongSnowIn, double Lv,
    double Press, double Rain, double NetShortUnder, double Vpd, double Wind,
    double OldTSurf, double IceDepth, double IceWE, double Tair, double TGrnd,
    double *AdvectedEnergy,
    double *DeltaColdContent, double *GroundFlux, double *LatentHeat,
    double *LatentHeatSub, double *NetLongUnder, double *SensibleHeat,
    double *vapor_flux, char *ErrorString);

/*****************************************************************************
  Function name: glacier_melt()

  Purpose      : Calculate glacier accumulation and melt using an energy balance
                 approach for a two layer snow model

  Required     :
    double delta_t               - Model timestep (secs)
    double z2           - Reference height (m)
    double displacement          - Displacement height (m)
    double aero_resist           - Aerodynamic resistance (uncorrected for
                                   stability) (s/m)
    double *aero_resist_used     - Aerodynamic resistance (corrected for
                                   stability) (s/m)
    double atmos->density        - Density of air (kg/m3)
    double atmos->vp             - Actual vapor pressure of air (Pa)
    double Le           - Latent heat of vaporization (J/kg3)
    double atmos->net_short      - Net exchange of shortwave radiation (W/m2)
    double atmos->longwave       - Incoming long wave radiation (W/m2)
    double atmos->pressure       - Air pressure (Pa)
    double RainFall              - Amount of rain (m)
    double Snowfall              - Amount of snow (m)
    double atmos->air_temp       - Air temperature (C)
    double atmos->vpd            - Vapor pressure deficit (Pa)
    double wind                  - Wind speed (m/s)
    double snow->pack_water      - Liquid water content of snow pack
    double snow->surf_water  - Liquid water content of surface layer
    double snow->swq             - Snow water equivalent at current pixel (m)
    double snow->vapor_flux;     - Mass flux of water vapor to or from the
                                   intercepted snow (m/time step)
    double snow->pack_temp       - Temperature of snow pack (C)
    double snow->surf_temp       - Temperature of snow pack surface layer (C)
    double snow->melt_energy     - Energy used for melting and heating of
                                   snow pack (W/m2)

  Modifies     :
    double *melt                 - Amount of snowpack outflow (initially is m, but converted to mm for output)
    double snow->pack_water      - Liquid water content of snow pack
    double snow->surf_water  - Liquid water content of surface layer
    double snow->swq             - Snow water equivalent at current pixel (m)
    double snow->vapor_flux;     - Mass flux of water vapor to or from the
                                   intercepted snow (m/time step)
    double snow->pack_temp       - Temperature of snow pack (C)
    double snow->surf_temp       - Temperature of snow pack surface layer (C)
    double snow->melt_energy     - Energy used for melting and heating of
                                   snow pack (W/m2)*/

int glacier_melt(double Le,
    double NetShort,  // net SW at absorbed by glacier
    double Tgrnd,
    VegConditions& roughness,  // roughness
    double aero_resist,  // aerodynamic resistance
    AeroResistUsed& aero_resist_used,  // stability-corrected aerodynamic resistance
    double air_temp,  // air temperature
    double delta_t,  // time step in secs
    double density,  // atmospheric density
    double displacement,  // surface displacement
    double LongIn,  // incoming longwave radiation
    double pressure, double rainfall,
    double vp,
    double vpd,
    double wind,
    double z2,
    double *NetLong,
    double *OldTSurf,
    double *melt,
    double *save_Qnet,
    double *save_advection,
    double *save_deltaCC_glac,
    double *save_grnd_flux,
    double *save_latent,
    double *save_latent_sub,
    double *save_sensible,
    int rec,
    glac_data_struct *glacier,
    const soil_con_struct* soil,
    const ProgramState *state)
{
  double error;
  double MassBalanceError;       /* Mass balance error (m) */
  double Qnet;                   /* Net energy exchange at the surface (W/m2) */
  double GlacMelt;               /* Amount of ice melt during time interval (m water equivalent) */
  double GlacCC;                 /* Cold content of glacier surface layer (J) */
  double RainFall;
  double advection;
  double deltaCC_glac;
  double latent_heat;
  double latent_heat_sub;
  double sensible_heat;
  double melt_energy = 0.;
  double grnd_flux;

  char ErrorString[MAXSTRING];

  /* SnowFall = snowfall / 1000.;*/ /* convert to m */
 RainFall = rainfall / 1000.; /* convert to m */

  (*OldTSurf) = glacier->surf_temp;

  /* Calculate the surface energy balance for surf_temp = 0.0 */
  GlacierEnergyBalance glacierEnergy(delta_t, aero_resist, aero_resist_used,
         displacement, z2, roughness,
         density, vp, LongIn, Le, pressure,
         RainFall, NetShort, vpd,
         wind, (*OldTSurf),
         soil->GLAC_SURF_THICK,
         soil->GLAC_SURF_WE, air_temp, Tgrnd,
         &advection,
         &deltaCC_glac,
         &grnd_flux, &latent_heat,
         &latent_heat_sub, NetLong,
         &sensible_heat,
         &glacier->vapor_flux);

  Qnet = glacierEnergy.calculate((double)0.0);

  /* If Qnet == 0.0, then set the surface temperature to 0.0 */
  /* if (std::abs(Qnet) < 2e-7) { */
  if (Qnet == 0.0) {
    glacier->surf_temp = 0.;
    melt_energy = NetShort + (*NetLong) + sensible_heat
        + latent_heat + latent_heat_sub + advection - deltaCC_glac;
    GlacMelt = melt_energy / (Lf * RHO_W) * delta_t;
    GlacCC = 0.;
  }

  /* Else, GlacierEnergyBalance(T=0.0) <= 0.0 */
  else {
    /* Calculate surface layer temperature using "Brent method" */
    GlacierEnergyBalance glacierIterative(delta_t, aero_resist, aero_resist_used,
        displacement, z2, roughness,
        density, vp, LongIn, Le, pressure,
        RainFall, NetShort, vpd,
        wind, (*OldTSurf),
        soil->GLAC_SURF_THICK,
        soil->GLAC_SURF_WE, air_temp, Tgrnd,
        &advection,
        &deltaCC_glac,
        &grnd_flux, &latent_heat,
        &latent_heat_sub, NetLong,
        &sensible_heat,
        &glacier->vapor_flux);

    glacier->surf_temp = glacierIterative.root_brent(
        (double) (glacier->surf_temp - SNOW_DT),
        (double) (glacier->surf_temp + SNOW_DT), ErrorString);

    if (glacierIterative.resultIsError(glacier->surf_temp)) {
      if (state->options.TFALLBACK) {
        glacier->surf_temp = *OldTSurf;
        glacier->surf_temp_fbflag = 1;
        glacier->surf_temp_fbcount++;
      } else {
        error = ErrorPrintGlacierEnergyBalance(glacier->surf_temp, rec,
            delta_t, aero_resist, aero_resist_used, displacement, z2, roughness,
            density, vp, LongIn, Le, pressure, RainFall, NetShort, vpd, wind,
            (*OldTSurf), soil->GLAC_SURF_THICK, soil->GLAC_SURF_WE, air_temp, Tgrnd,
            &advection, &deltaCC_glac, &grnd_flux,
            &latent_heat, &latent_heat_sub, NetLong,
            &sensible_heat, &glacier->vapor_flux, ErrorString);
        return (ERROR);
      }
    }

    if (glacierIterative.resultIsError(glacier->surf_temp) == false) {  // Result is valid
      GlacierEnergyBalance glacierEnergy(delta_t, aero_resist, aero_resist_used,
          displacement, z2, roughness,
          density, vp, LongIn, Le, pressure,
          RainFall, NetShort, vpd,
          wind, (*OldTSurf),
          soil->GLAC_SURF_THICK,
          soil->GLAC_SURF_WE, air_temp, Tgrnd,
          &advection,
          &deltaCC_glac,
          &grnd_flux, &latent_heat,
          &latent_heat_sub, NetLong,
          &sensible_heat,
          &glacier->vapor_flux);

      Qnet = glacierEnergy.calculate(glacier->surf_temp);

      /* since we iterated, the surface layer is below freezing and no snowmelt */
      GlacMelt = 0.0;
      GlacCC = CH_ICE * glacier->surf_temp * soil->GLAC_SURF_WE;

    }
  }

  melt[0] = GlacMelt;

  /* Mass balance test */
  /* MassBalanceError = (InitialSwq - snow->swq) + (RainFall + SnowFall)
   - melt[0] + snow->vapor_flux; */

  /*  printf("%d %d %g\n", y, x, MassBalanceError);*/

  /*melt[0] *= 1000.;*/ /* converts back to mm */
  glacier->cold_content = GlacCC;
  glacier->vapor_flux *= -1.;
  *save_advection = advection;
  *save_deltaCC_glac = deltaCC_glac;
  *save_grnd_flux = grnd_flux;
  *save_latent = latent_heat;
  *save_latent_sub = latent_heat_sub;
  *save_sensible = sensible_heat;
  *save_Qnet = Qnet;

  return (0);
}

double ErrorPrintGlacierEnergyBalance(double TSurf,
    /* General Model Parameters */
    int rec,
    double Dt,                      /* Model time step (sec) */
    double Ra,                      /* Aerodynamic resistance (s/m) */
    AeroResistUsed& Ra_used,        /* Aerodynamic resistance (s/m) after stability correction */
    /* Vegetation Parameters */
    double Displacement,            /* Displacement height (m) */
    double Z,                       /* Reference height (m) */
    VegConditions& roughness,       /* surface roughness height (m) */
    /* Atmospheric Forcing Variables */
    double AirDens,                 /* Density of air (kg/m3) */
    double EactAir,                 /* Actual vapor pressure of air (Pa) */
    double LongSnowIn,              /* Incoming longwave radiation (W/m2) */
    double Lv,                      /* Latent heat of vaporization (J/kg3) */
    double Press,                   /* Air pressure (Pa) */
    double Rain,                    /* Rain fall (m/timestep) */
    double NetShortUnder,           /* Net incident shortwave radiation (W/m2) */
    double Vpd,                     /* Vapor pressure deficit (Pa) */
    double Wind,                    /* Wind speed (m/s) */
    /* Snowpack Variables */
    double OldTSurf,                /* Surface temperature during previous time step */
    double IceDepth,                /* Depth of glacier surface layer (m) */
    double IceWE,                   /* Liquid water in the glacier surface layer (m) */
    /* Energy Balance Components */
    double Tair,                    /* Canopy air / Air temperature (C) */
    double TGrnd,                   /* Ground surface temperature (C) */
    double *AdvectedEnergy,         /* Energy advected by precipitation (W/m2) */
    double *DeltaColdContent,       /* Change in cold content of surface layer (W/m2) */
    double *GroundFlux,             /* Ground Heat Flux (W/m2) */
    double *LatentHeat,             /* Latent heat exchange at surface (W/m2) */
    double *LatentHeatSub,          /* Latent heat of sublimation exchange at surface (W/m2) */
    double *NetLongUnder,           /* Net longwave radiation at snowpack surface (W/m^2) */
    double *SensibleHeat,           /* Sensible heat exchange at surface (W/m2) */
    double *vapor_flux,             /* Mass flux of water vapor to or from the intercepted snow (m/timestep) */
    char *ErrorString)
{

  /* print variables */
  fprintf(stderr, "%s", ErrorString);
  fprintf(stderr, "ERROR: glacier_melt failed to converge to a solution in root_brent.  Variable values will be dumped to the screen, check for invalid values.\n");

  /* general model terms */
  fprintf(stderr, "TSurf = %f\n", TSurf);
  fprintf(stderr, "rec = %i\n", rec);
  fprintf(stderr, "Dt = %f\n",Dt);

  /* land surface parameters */
  fprintf(stderr,"Ra = %f\n",Ra);
  fprintf(stderr, "Ra_used = %f\n", Ra_used.surface);
  fprintf(stderr,"Displacement = %f\n",Displacement);
  fprintf(stderr,"Z = %f\n",Z);
  fprintf(stderr,"Z0 = %f\n",roughness.snowFree);

  /* meteorological terms */
  fprintf(stderr,"AirDens = %f\n",AirDens);
  fprintf(stderr,"EactAir = %f\n",EactAir);
  fprintf(stderr,"LongIn = %f\n",LongSnowIn);
  fprintf(stderr,"Lv = %f\n",Lv);
  fprintf(stderr,"Press = %f\n",Press);
  fprintf(stderr,"Rain = %f\n",Rain);
  fprintf(stderr,"ShortRad = %f\n",NetShortUnder);
  fprintf(stderr,"Vpd = %f\n",Vpd);
  fprintf(stderr,"Wind = %f\n",Wind);

  /* glacier terms */
  fprintf(stderr,"OldTSurf = %f\n",OldTSurf);
  fprintf(stderr, "IceDepth = %f\n", IceDepth);
  fprintf(stderr, "IceWE = %f\n", IceWE);
  fprintf(stderr,"Tair = %f\n",Tair);
  fprintf(stderr,"TGrnd = %f\n",TGrnd);
  fprintf(stderr,"AdvectedEnergy = %f\n",AdvectedEnergy[0]);
  fprintf(stderr,"DeltaColdContent = %f\n",DeltaColdContent[0]);
  fprintf(stderr,"GroundFlux = %f\n",GroundFlux[0]);
  fprintf(stderr,"LatentHeat = %f\n",LatentHeat[0]);
  fprintf(stderr,"LatentHeatSub = %f\n",LatentHeatSub[0]);
  fprintf(stderr,"NetLong = %f\n",NetLongUnder[0]);
  fprintf(stderr,"SensibleHeat = %f\n",SensibleHeat[0]);
  fprintf(stderr,"VaporMassFlux = %f\n",vapor_flux[0]);

  fprintf(stderr,"Finished dumping glacier_melt variables.\nTry increasing SNOW_DT to get model to complete cell.\nThen check output for instabilities.\n");

  return(ERROR);

}
