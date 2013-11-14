
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

  double Density;                 /* Density of water/ice at TMean (kg/m3) */
  double NetRad;                  /* Net radiation exchange at surface (W/m2) */
  double RestTerm;                /* Rest term in surface energy balance (W/m2) */
  double TMean;                   /* Average temperature for time step (C) */
  double Tmp;
  double VaporMassFlux;           /* Mass flux of water vapor to or from the intercepted snow (kg/m2s) */
  double BlowingMassFlux;         /* Mass flux of water vapor from blowing snow. (kg/m2s) */
  double SurfaceMassFlux;         /* Mass flux of water vapor from pack snow. (kg/m2s) */


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
  *GroundFlux = (GLAC_K_ICE + TMean*(-0.0142)) * (TGrnd - TMean) / IceDepth / (Dt);

  /* Calculate energy balance error at the snowpack surface */
  double Fbal = NetRad + *SensibleHeat + *LatentHeat + *LatentHeatSub + *AdvectedEnergy + *AdvectedSensibleHeat;
  RestTerm = Fbal - *DeltaColdContent + *GroundFlux;

  if (TSurf == 0.0 && Fbal >= 0.) RestTerm = 0.;

  return RestTerm;
}
