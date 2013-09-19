/*
 * SUMMARY:      IceEnergyBalance.c - Calculate lake ice energy balance
 * USAGE:        Part of the lake algorithm
 *
 * AUTHOR:       Laura Bowling
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              nijssen@u.washington.edu
 * ORIG-DATE:     March 16, 2001
 * LAST-MOD: Mon Apr 21 15:51:12 2003 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Calculate ice energy balance
 * DESCRIP-END.
 * FUNCTIONS:    IceEnergyBalance()
 * COMMENTS:     
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include "IceEnergyBalance.h"

static char vcid[] = "$Id$";

/*****************************************************************************
  Function name: IceEnergyBalance()

  Purpose      : Calculate the surface energy balance for the snow pack

  Required     :
    double TSurf           - new estimate of effective surface temperature

  Returns      :
    double RestTerm        - Rest term in the energy balance

  Modifies     : 
    double *RefreezeEnergy - Refreeze energy (W/m2) 
    double *VaporMassFlux  - Mass flux of water vapor to or from the
                            intercepted snow (m/s)

  Comments     :
    Reference:  Bras, R. A., Hydrology, an introduction to hydrologic
                science, Addisson Wesley, Inc., Reading, etc., 1990.

  Modifications:
  16-Jul-04 Renamed VaporMassFlux to vapor_flux to denote fact that its
	    units are m/timestep rather than kg/m2s.  Created new variable
	    VaporMassFlux with units of kg/m2s.  After VaporMassFlux is
	    computed, vapor_flux is derived from it by unit conversion.
	    vapor_flux is the variable that is passed in/out of this
	    function.							TJB
  24-Aug-04 Added logic to handle blowing_flux and vapor_flux.		TJB
  28-Sep-04 Added Ra_used to store the aerodynamic resistance used in flux
	    calculations.						TJB
  04-Oct-04 Merged with Laura Bowling's updated lake model code.  Now
	    blowing snow sublimation is calculated for lakes.		TJB
  2006-Sep-23 Replaced redundant STEFAN constant with STEFAN_B.  TJB
  2006-Nov-07 Removed LAKE_MODEL option. TJB

*****************************************************************************/
double IceEnergyBalance::calculate(double TSurf)
{

  const char *Routine = "IceEnergyBalance";

  double Density;                /* Density of water/ice at TMean (kg/m3) */
  double EsSnow;                 /* saturated vapor pressure in the snow pack
                                   (Pa)  */  
  
  double Ls;                     /* Latent heat of sublimation (J/kg) */
  double NetRad;			/* Net radiation exchange at surface (W/m2) */
  double RestTerm;		/* Rest term in surface energy balance
				   (W/m2) */
  double  TMean;                /* Average temperature for time step (C) */
  double qnull;
  double VaporMassFlux;          /* Total mass flux of water vapor to or from
                                    snow (kg/m2s) */
  double BlowingMassFlux;        /* Mass flux of water vapor to or from
                                    blowing snow (kg/m2s) */
  double SurfaceMassFlux;        /* Mass flux of water vapor to or from
                                    snow pack (kg/m2s) */
  
  /* Calculate active temp for energy balance as average of old and new  */
  
/*   TMean = 0.5 * (OldTSurf + TSurf); */
  TMean = TSurf;
  Density = RHO_W;
  
  /* Correct aerodynamic conductance for stable conditions
     Note: If air temp >> snow temp then aero_cond -> 0 (i.e. very stable)
     velocity (vel_2m) is expected to be in m/sec */                        
  
  /* Apply the stability correction to the aerodynamic resistance 
     NOTE: In the old code 2m was passed instead of Z-Displacement.  I (bart)
     think that it is more correct to calculate ALL fluxes at the same
     reference level */


  if (Wind > 0.0)
    *Ra_used = Ra / StabilityCorrection(Z, 0.f, TMean, Tair, Wind, Z0);
    /* *Ra_used = Ra / StabilityCorrection(2.f, 0.f, TMean, Tair, Wind, Z0);*/
  else
    *Ra_used = HUGE_RESIST;

  /* Calculate longwave exchange and net radiation */
 
  *LongRadOut = LongRadIn - STEFAN_B * (TMean+273.15) * (TMean+273.15) 
    * (TMean+273.15) * (TMean+273.15);
  NetRad = ShortRad + *LongRadOut;
  
  /* Calculate the sensible heat flux */

  *SensibleHeat = AirDens * CP_PM * (Tair - TMean) / *Ra_used;

  /* Calculate the mass flux of ice to or from the surface layer */
 
  /* Calculate the saturated vapor pressure in the snow pack, 
     (Equation 3.32, Bras 1990) */

  EsSnow = svp(TMean) /* * 1000. */;

/*   EsSnow = 610.78 * exp((double)((17.269 * TMean) / (237.3 + TMean))); */

/*   if (TMean < 0.0) */
/*     EsSnow *= 1.0 + .00972 * TMean + .000042 * pow((double)TMean,(double)2.0); */
  
  /* Calculate sublimation terms and latent heat flux */

  /* blowing_flux was calculated outside of the root_brent iteration */
  BlowingMassFlux = *blowing_flux * Density / (Dt * SECPHOUR);

  latent_heat_from_snow(AirDens, EactAir, Lv, Press, Ra, TMean, Vpd,
                        LatentHeat, LatentHeatSub, &VaporMassFlux, &BlowingMassFlux,
                        &SurfaceMassFlux);

  /* Convert sublimation terms from kg/m2s to m/timestep */
  *vapor_flux = VaporMassFlux * Dt * SECPHOUR / Density;
  *surface_flux = SurfaceMassFlux * Dt * SECPHOUR / Density;

  /* Calculate advected heat flux from rain 
     WORK IN PROGRESS:  Should the following read (Tair - Tsurf) ?? */
  
  // Temporary fix for lake model.
  *AdvectedEnergy = (CH_WATER * Tair * Rain) / (Dt*SECPHOUR);
  //*AdvectedEnergy = 0.0;
  
  /* Calculate change in cold content */
  /* No change in cold content in lake model */

  /* Changes for lake model start here. Actually, equals qo-Io (P&H eq. 7)*/
    /* because Io (SWnet) is included in NetRad below. */
  qnull = (1/AvgCond)*(Tfreeze - TMean + SWconducted);
  *qf   = qnull;

  /* Changes for lake model end here. */
  
  /* Calculate net energy exchange at the snow surface */
  
  RestTerm = ( NetRad + *SensibleHeat + *LatentHeat + *LatentHeatSub + *AdvectedEnergy 
	       + qnull );

  *RefreezeEnergy = (SurfaceLiquidWater * Lf * Density)/(Dt * SECPHOUR);

  /* Melting, or partially refreeze surface water. */
  if (TSurf == 0.0 && RestTerm > -(*RefreezeEnergy)) {
    *RefreezeEnergy = -RestTerm;  /* available energy input over cold content
                                    used to melt, i.e. Qrf is negative value
                                    (energy out of pack)*/ 
    RestTerm = 0.0;
  }
  else {                         /* Pack is getting colder. */
    RestTerm += *RefreezeEnergy; /* add this positive value to the pack */
  }
  
  return RestTerm;
}
