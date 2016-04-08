
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "GlacierEnergyBalance.h"
#include "vicNl.h"

/*****************************************************************************
  Function name: GlacierEnergyBalance()

  Solves the glacier energy balance. See SnowPackEnergyBalance for comparison.

*****************************************************************************/
double GlacierEnergyBalance::calculate(double TSurf) {

  const char *Routine = "GlacierEnergyBalance";

  /* Internal Routine Variables */
  double Density;                 /* Density of water/ice at TSurf (kg/m3) */
  double NetRad;                  /* Net radiation exchange at surface (W/m2) */
  double RestTerm;                /* Rest term in surface energy balance (W/m2) */
  double TMean;                   /* Average ice surface layer temperature for time step (C) */
  double OldTMean;
  double Tmp;
  double VaporMassFlux;           /* Mass flux of water vapor to or from the intercepted snow (kg/m2s) */
  double Fbal;                    /* Energy balance at glacier surface */
  double temp_IceDepth;			  /* Local variable to hold scaled value for IceDepth */


  TMean = (TSurf + TGrnd)/2;
  OldTMean = (OldTSurf + TGrnd)/2;
  Density = RHO_W;
  temp_IceDepth = IceDepth/1000.;	/* convert IceDepth to mm */

  /* Correct aerodynamic conductance for stable conditions
     Note: If air temp >> glacier temp then aero_cond -> 0 (i.e. very stable)
     velocity (vel_2m) is expected to be in m/sec */

  /* Apply the stability correction to the aerodynamic resistance
     NOTE: In the old code 2m was passed instead of Z-Displacement.  I (bart)
     think that it is more correct to calculate ALL fluxes at the same
     reference level */
  if (Wind > 0.0)
    Ra_used.surface = Ra / StabilityCorrection(Z, 0.f, TSurf, Tair, Wind, roughness.snowCovered);
  else
    Ra_used.surface = HUGE_RESIST;

  /* Calculate longwave exchange and net radiation */

  Tmp = TSurf + KELVIN;
  (*NetLongUnder) = LongSnowIn - STEFAN_B * Tmp * Tmp * Tmp * Tmp;
  NetRad = NetShortUnder + (*NetLongUnder);

  /* Calculate the sensible heat flux */
  *SensibleHeat = AirDens * Cp * (Tair - TSurf) / Ra_used.surface;

  /* Convert sublimation terms from m/timestep to kg/m2s */
  VaporMassFlux = *vapor_flux * Density / Dt;

  /* Calculate the mass flux of ice to or from the surface layer */

  /* Calculate the saturated vapor pressure,
     (Equation 3.32, Bras 1990) */
   latent_heat_from_glacier(AirDens, EactAir, Lv, Press, Ra_used.surface, TSurf, Vpd,
      LatentHeat, LatentHeatSub, &VaporMassFlux);

  /* Convert sublimation terms from kg/m2s to m/timestep */
  *vapor_flux = VaporMassFlux * Dt / Density;

  /* Calculate advected heat flux from rain
     Equation 7.3.12 from H.B.H. for rain falling on melting snowpack */
  if ( TSurf == 0 )
    *AdvectedEnergy = (CH_WATER * (Tair) * Rain) / (Dt);
  else
    *AdvectedEnergy = 0.;

  /* Calculate change in cold content */
  *DeltaColdContent = CH_ICE * temp_IceDepth * (TMean - OldTMean) / (Dt);

  /* Calculate Ground Heat Flux */
  /* Estimate of ice thermal conductivity (at atmospheric pressure) adapted from Slack (1980), Table 1; assumes
       linear relationship between TSurf and K above -75C */
  *GroundFlux = (GLAC_K_ICE + TSurf*(-0.0142)) * (TGrnd - TSurf) / temp_IceDepth;

  /* Calculate energy balance error at the glacier surface */
  Fbal = NetRad + *SensibleHeat + *LatentHeat + *LatentHeatSub + *AdvectedEnergy;
  RestTerm = Fbal - *DeltaColdContent + *GroundFlux;

  /* Melting occurs when surface at melting point and surface energy flux is positive */
  if (TSurf == 0.0 && RestTerm >= 0.) RestTerm = 0.;

  return RestTerm;
}
