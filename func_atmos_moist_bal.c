#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vicNl.h"

static char vcid[] = "$Id$";

double func_atmos_moist_bal(double VPcanopy, double InLatentHeat, double Lv,
    double Ra, double atmos_density, double gamma, double atmospheric_vp // atmospheric vapor pressure
    ) {
/**********************************************************************
  func_atmos_moist_bal.c      Keith Cherkauer        March 2, 2001

  This routine solves the atmospheric exchange moisture balance.

**********************************************************************/
  double LatentHeat;
 
  // internal routine variables
  double  Error;

  // compute sensible heat flux between canopy and atmosphere
  LatentHeat = Lv * atmos_density * Cp * (atmospheric_vp - VPcanopy) / ( gamma * Ra );

  // compute energy balance error
  Error = InLatentHeat - LatentHeat;

  return ( Error );

}
