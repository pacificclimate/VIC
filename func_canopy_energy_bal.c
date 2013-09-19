#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include "canopy_energy_bal.h"

static char vcid[] = "$Id$";

double CanopyEnergyBal::calculate(double Tfoliage)
/*********************************************************************
  func_canopy_energy_bal    Keith Cherkauer         January 27, 2001

  This routine iterates to determine the temperature of the canopy,
  and solve the resulting fluxes between the canopy and the atmosphere 
  and the canopy and the ground.

  Modifications:
  2004-Sep-28 Added Ra_used to store the aerodynamic resistance used in
	      flux calculations.					TJB
  2009-Jan-16 Modified aero_resist_used and Ra_used to become arrays of
	      two elements (surface and overstory); added
	      options.AERO_RESIST_CANSNOW.				TJB
  2009-May-17 Added AR_406_LS to options.AERO_RESIST_CANSNOW.		TJB
  2009-Sep-14 Replaced 0.622 with EPS in equation for vapor flux.	TJB

 ********************************************************************/
{
  /* Internal Variables */
  double  EsSnow;
  double  LongOut;
  double  Ls;
  double  RestTerm;
  double  Tmp;
  double  prec[2];

  /* Calculate the net radiation at the canopy surface, using the canopy 
     temperature.  The outgoing longwave is subtracted twice, because the 
     canopy radiates in two directions */

  Tmp = Tfoliage + KELVIN;
  *LongOverOut = STEFAN_B * (Tmp * Tmp * Tmp * Tmp);
  *NetRadiation = NetShortOver + LongOverIn + LongUnderOut 
    - 2 * (*LongOverOut);

  *NetLongOver  = LongOverIn - (*LongOverOut);

  if ( IntSnow > 0 ) {

    Ra_used[0] = Ra[0];
    Ra_used[1] = Ra[1];

    /** Added multiplication by 10 to incorporate change in canopy resistance due
	to smoothing by intercepted snow **/
    if (state->options.AERO_RESIST_CANSNOW == AR_COMBO || state->options.AERO_RESIST_CANSNOW == AR_406 || state->options.AERO_RESIST_CANSNOW == AR_406_LS || state->options.AERO_RESIST_CANSNOW == AR_406_FULL)
      Ra_used[1] *= 10.;

    /** Calculate the vapor mass flux between intercepted snow in 
	the canopy and the surrounding air mass **/
  
    EsSnow = svp(Tfoliage); 
    
    /* Apply stability correction to aerodynamic resistance */
    if (state->options.AERO_RESIST_CANSNOW == AR_COMBO || state->options.AERO_RESIST_CANSNOW == AR_410) {
      if (Wind[1] > 0.0) {
        Ra_used[1] /= StabilityCorrection(ref_height[1], displacement[1], Tfoliage, 
				          Tcanopy, Wind[1], roughness[1]);
      }
      else
        Ra_used[1] = HUGE_RESIST;
    }

    *VaporMassFlux = AirDens * ( EPS / Press ) * (EactAir - EsSnow) 
      / Ra_used[1] / RHO_W; 

    if (Vpd == 0.0 && *VaporMassFlux < 0.0)
      *VaporMassFlux = 0.0;
  
    /* Calculate the latent heat flux */
    
    Ls = (677. - 0.07 * Tfoliage) * JOULESPCAL * GRAMSPKG;
    *LatentHeatSub = Ls * *VaporMassFlux * RHO_W;
    *LatentHeat = 0;
    *Evap = 0;
    veg_var_wet->throughfall = 0;

    if (state->options.AERO_RESIST_CANSNOW == AR_406)
      Ra_used[1] /= 10;
  }
  else {

    if (state->options.AERO_RESIST_CANSNOW == AR_406_FULL || state->options.AERO_RESIST_CANSNOW == AR_410 || state->options.AERO_RESIST_CANSNOW == AR_COMBO) {
      Ra_used[0] = Ra[0];
      Ra_used[1] = Ra[1];
    }
    else {
      Ra_used[0] = Ra[0];
      Ra_used[1] = Ra[0];
    }

    Wdew[WET] = IntRain * 1000.;
    prec[WET] = *Rainfall * 1000;
    prec[DRY] = 0;
    *Evap = canopy_evap(layer_wet, layer_dry, veg_var_wet, veg_var_dry, FALSE, 
			veg_class, month, precipitation_mu, Wdew, delta_t, *NetRadiation, 
			Vpd, NetShortOver, Tcanopy, Ra_used[1], displacement[1], 
			roughness[1], ref_height[1], elevation, prec, 
			depth, Wcr, Wpwp, frost_fract, root, state);
    Wdew[WET] /= 1000.;

    *LatentHeat = latent_heat_Le * *Evap * RHO_W;
    *LatentHeatSub = 0;

  }

  /* Calculate the sensible heat flux */

  *SensibleHeat = AirDens * Cp * (Tcanopy - Tfoliage) / Ra_used[1];

  /* Calculate the advected energy */

  *AdvectedEnergy = (4186.8 * Tcanopy * Rainfall[0]) / (delta_t);

  /* Calculate the amount of energy available for refreezing */
  
  RestTerm = *SensibleHeat + *LatentHeat + *LatentHeatSub + *NetRadiation 
    + *AdvectedEnergy;
  
  if ( IntSnow > 0 ) {
    /* Intercepted snow present, check if excess energy can be used to 
       melt or refreeze it */

    *RefreezeEnergy = (IntRain * Lf * RHO_W) / (delta_t);

    if (Tfoliage == 0.0 && RestTerm > -(*RefreezeEnergy)) {
      *RefreezeEnergy = -RestTerm;  /* available energy input over cold content
				       used to melt, i.e. Qrf is negative value
				       (energy out of pack)*/ 
      RestTerm = 0.0;
    }
    else {
      RestTerm += *RefreezeEnergy; /* add this positive value to the pack */
    }

  }
  else *RefreezeEnergy = 0;

  return (RestTerm);

}
