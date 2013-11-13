
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

/*****************************************************************************
  Function name: GlacierEnergyBalance()

  Solves the glacier energy balance. See SnowPackEnergyBalance for comparison.

*****************************************************************************/
double GlacierEnergyBalance(double TSurf, va_list ap) {

  extern option_struct options;

  const char *Routine = "GlacierEnergyBalance";

  /* Define Variable Argument List */

  /* General Model Parameters */
  double Dt;                      /* Model time step (sec) */
  double Ra;                      /* Aerodynamic resistance (s/m) */
  double *Ra_used;                /* Aerodynamic resistance (s/m) after stability correction */

  /* Vegetation Parameters */
  double Displacement;            /* Displacement height (m) */
  double Z;                       /* Reference height (m) */
  double *Z0;                      /* surface roughness height (m) */

  /* Atmospheric Forcing Variables */
  double AirDens;                 /* Density of air (kg/m3) */
  double EactAir;                 /* Actual vapor pressure of air (Pa) */
  double LongSnowIn;               /* Incoming longwave radiation (W/m2) */
  double Lv;                      /* Latent heat of vaporization (J/kg3) */
  double Press;                   /* Air pressure (Pa) */
  double Rain;                    /* Rain fall (m/timestep) */
  double NetShortUnder;           /* Net incident shortwave radiation
             (W/m2) */
  double Vpd;       /* Vapor pressure deficit (Pa) */
  double Wind;                    /* Wind speed (m/s) */

  /* Snowpack Variables */
  double OldTSurf;                /* Surface temperature during previous time
             step */
  double IceDepth;               /* Depth of glacier surface layer (m) */
  double IceWE;      /* Liquid water in the glacier surface layer (m) */

  /* Energy Balance Components */
  double Tair;                    /* Canopy air / Air temperature (C) */
  double TGrnd;                   /* Ground surface temperature (C) */

  double *AdvectedEnergy;         /* Energy advected by precipitation (W/m2) */
  double *AdvectedSensibleHeat;   /* Sensible heat advected from snow-free
             area into snow covered area (W/m^2) */
  double *DeltaColdContent;       /* Change in cold content of surface
             layer (W/m2) */
  double *GroundFlux;     /* Ground Heat Flux (W/m2) */
  double *LatentHeat;     /* Latent heat exchange at surface (W/m2) */
  double *LatentHeatSub;  /* Latent heat of sublimation exchange at
             surface (W/m2) */
  double *NetLongUnder;           /* Net longwave radiation at snowpack
             surface (W/m^2) */
  double *SensibleHeat;     /* Sensible heat exchange at surface
             (W/m2) */
  double *vapor_flux;             /* Mass flux of water vapor to or from the
             intercepted snow (m/timestep) */
  double *blowing_flux;           /* Mass flux of water vapor from blowing snow. (m/timestep) */
  double *surface_flux;           /* Mass flux of water vapor from pack snow. (m/timestep) */

  /* Internal Routine Variables */

  double Density;                 /* Density of water/ice at TMean (kg/m3) */
  /* double LongRadOut; */      /* long wave radiation emitted by surface
             (W/m2) */
  double NetRad;      /* Net radiation exchange at surface
             (W/m2) */
  double RestTerm;      /* Rest term in surface energy balance
             (W/m2) */
  double TMean;                   /* Average temperature for time step (C) */
  double Tmp;
  double VaporMassFlux;           /* Mass flux of water vapor to or from the
             intercepted snow (kg/m2s) */
  double BlowingMassFlux;         /* Mass flux of water vapor from blowing snow. (kg/m2s) */
  double SurfaceMassFlux;         /* Mass flux of water vapor from pack snow. (kg/m2s) */

  /* Assign the elements of the array to the appropriate variables.  The list
     is traversed as if the elements are doubles, because:

     In the variable-length part of variable-length argument lists, the old
     ``default argument promotions'' apply: arguments of type double are
     always promoted (widened) to type double, and types char and short int
     are promoted to int. Therefore, it is never correct to invoke
     va_arg(argp, double); instead you should always use va_arg(argp,
     double).

     (quoted from the comp.lang.c FAQ list)
     */

  /* General Model Parameters */
  Dt           = (double) va_arg(ap, double);
  Ra           = (double) va_arg(ap, double);
  Ra_used      = (double *) va_arg(ap, double *);

  /* Vegetation Parameters */
  Displacement = (double) va_arg(ap, double);
  Z            = (double) va_arg(ap, double);
  Z0           = (double *) va_arg(ap, double *);

  /* Atmospheric Forcing Variables */
  AirDens       = (double) va_arg(ap, double);
  EactAir       = (double) va_arg(ap, double);
  LongSnowIn    = (double) va_arg(ap, double);
  Lv            = (double) va_arg(ap, double);
  Press         = (double) va_arg(ap, double);
  Rain          = (double) va_arg(ap, double);
  NetShortUnder = (double) va_arg(ap, double);
  Vpd           = (double) va_arg(ap, double);
  Wind          = (double) va_arg(ap, double);

  /* Glacier Variables */
  OldTSurf           = (double) va_arg(ap, double);
  IceDepth          = (double) va_arg(ap, double);
  IceWE     = (double) va_arg(ap, double);

  /* Energy Balance Components */
  Tair = (double) va_arg(ap, double);
  TGrnd   = (double) va_arg(ap, double);

  AdvectedEnergy        = (double *) va_arg(ap, double *);
  AdvectedSensibleHeat  = (double *)va_arg(ap, double *);
  DeltaColdContent      = (double *) va_arg(ap, double *);
  GroundFlux            = (double *) va_arg(ap, double *);
  LatentHeat            = (double *) va_arg(ap, double *);
  LatentHeatSub         = (double *) va_arg(ap, double *);
  NetLongUnder          = (double *) va_arg(ap, double *);
  SensibleHeat          = (double *) va_arg(ap, double *);
  vapor_flux            = (double *) va_arg(ap, double *);
  blowing_flux          = (double *) va_arg(ap, double *);
  surface_flux          = (double *) va_arg(ap, double *);

  /* Calculate active temp for energy balance as average of old and new  */

  TMean = TSurf;
  Density = RHO_W;

  /* Correct aerodynamic conductance for stable conditions
     Note: If air temp >> glacier temp then aero_cond -> 0 (i.e. very stable)
     velocity (vel_2m) is expected to be in m/sec */

  /* Apply the stability correction to the aerodynamic resistance
     NOTE: In the old code 2m was passed instead of Z-Displacement.  I (bart)
     think that it is more correct to calculate ALL fluxes at the same
     reference level */


  if (Wind > 0.0)
    Ra_used[0] = Ra / StabilityCorrection(Z, 0.f, TMean, Tair, Wind, Z0[2]);
  else
    Ra_used[0] = HUGE_RESIST;

  /* Calculate longwave exchange and net radiation */

  Tmp = TMean + KELVIN;
  (*NetLongUnder) = LongSnowIn - STEFAN_B * Tmp * Tmp * Tmp * Tmp;
  NetRad = NetShortUnder + (*NetLongUnder);

  /* Calculate the sensible heat flux */

  *SensibleHeat = AirDens * Cp * (Tair - TMean) / Ra_used[0];

  (*AdvectedSensibleHeat) = 0;

  /* Convert sublimation terms from m/timestep to kg/m2s */
  VaporMassFlux = *vapor_flux * Density / Dt;
  BlowingMassFlux = *blowing_flux * Density / Dt; /* blowing_flux = 0. for a bare glacier */
  SurfaceMassFlux = *surface_flux * Density / Dt;

  /* Calculate the mass flux of ice to or from the surface layer */

  /* Calculate the saturated vapor pressure,
     (Equation 3.32, Bras 1990) */

  latent_heat_from_snow(AirDens, EactAir, Lv, Press, Ra_used[0], TMean, Vpd,
      LatentHeat, LatentHeatSub, &VaporMassFlux, &BlowingMassFlux,
      &SurfaceMassFlux);

  /* Convert sublimation terms from kg/m2s to m/timestep */
  *vapor_flux = VaporMassFlux * Dt / Density;
  *blowing_flux = BlowingMassFlux * Dt / Density;
  *surface_flux = SurfaceMassFlux * Dt / Density;

  /* Calculate advected heat flux from rain
     Equation 7.3.12 from H.B.H. for rain falling on melting snowpack */

  if ( TMean == 0 )
    *AdvectedEnergy = (CH_WATER * (Tair) * Rain) / (Dt);
  else
    *AdvectedEnergy = 0.;

  /* Calculate change in cold content */
  *DeltaColdContent = CH_ICE * IceWE * (TSurf - OldTSurf) / (Dt);

  /* Calculate Ground Heat Flux */
  /* Estimate of ice thermal conductivity (at atmospheric pressure) adapted from Slack (1980), Table 1; assumes
       linear relationship between Tmean and K below -75C */
  *GroundFlux = (K_ICE + Tmean*(-0.0142)) * (TGrnd - TMean) / IceDepth / (Dt);

  /* Calculate energy balance error at the snowpack surface */
  Fbal = NetRad + *SensibleHeat + *LatentHeat + *LatentHeatSub + *AdvectedEnergy + *AdvectedSensibleHeat;
  RestTerm = Fbal - *DeltaColdContent + *GroundFlux;

  if (TSurf == 0.0 && Fbal >= 0.) RestTerm = 0.;

  return RestTerm;
}
