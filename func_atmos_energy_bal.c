#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vicNl.h"
#include "atmos_energy_bal.h"

static char vcid[] = "$Id$";

double AtmosEnergyBal::calculate(double Tcanopy) {
/**********************************************************************
  func_atmos_energy_bal.c      Keith Cherkauer        February 6, 2001

  This routine solves the atmospheric exchange energy balance.

**********************************************************************/

 
  // internal routine variables
  double  Error;

  // compute sensible heat flux between canopy and atmosphere
  (*SensibleHeat) = atmos_density * Cp * (Tair - Tcanopy) / Ra;

  // compute energy balance error
  //Error = NetRadiation + LatentHeat + (*SensibleHeat);
  Error = InSensible - (*SensibleHeat);

  return ( Error );

}
