#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

static char vcid[] = "$Id: latent_heat_from_snow.c,v 5.5 2004/08/27 22:15:27 vicadmin Exp $";

void latent_heat_from_glacier(double  AirDens,
         double  EactAir,
         double  Lv,
         double  Press,
         double  Ra,
         double  TMean,
         double  Vpd,
         double *LatentHeat,
         double *LatentHeatSublimation,
         double *VaporMassFlux) {
/**********************************************************************
  latent_heat_from_snow.c       Laura Bowling

  Split out of the snowpack energy balance, this subroutine computes
  the latent heat from the snowpack.


***********************************************************************/

  double EsSnow;
  double Ls;

  EsSnow = svp(TMean);

  // SurfaceMassFlux and BlowingMassFlux in kg/m2s

  *VaporMassFlux = AirDens * ( EPS / Press ) * ( EactAir - EsSnow ) / Ra;

  if ( Vpd == 0.0 && *VaporMassFlux < 0.0 )
    *VaporMassFlux = 0.0;


  if ( TMean >= 0.0 ) {
    /* Melt conditions: use latent heat of vaporization */
    *LatentHeat = Lv * (*VaporMassFlux);
    *LatentHeatSublimation = 0;
  }
  else {
    /* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990 */
    Ls = (677. - 0.07 * TMean) * JOULESPCAL * GRAMSPKG;
    *LatentHeatSublimation = Ls * (*VaporMassFlux);
    *LatentHeat = 0;
  }

}
