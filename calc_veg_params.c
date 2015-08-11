#include <stdio.h>
#include "vicNl.h"
#include <math.h>

static char vcid[] = "$Id$";

double calc_veg_displacement(double height, double L) {
/**********************************************************************
  calc_veg_displacement		Markus Schnorbus		August 9, 2015

  This subroutine estimates the displacement height of vegetation
  as a function of vegetation height and canopy density (via leaf area
  index), as suggested by Choudhury and Monteith (1988).
**********************************************************************/

  double value;
  double X;

  X = COEF_DRAG*L;
  value = 1.1 * height * log(1 + pow(X,0.25));

  return (value);

}

double calc_veg_height(double displacement, double L) {
/**********************************************************************
  calc_veg_height		Markus Schnorbus		August 9, 2015

  This subroutine backs the vegetation height out of the given
  displacement by simply reversing calc_veg_displacement().
**********************************************************************/

  double value;
  double X;

  X = COEF_DRAG*L;
  value = displacement / (1.1 * log(1 + pow(X,0.25)));

  return (value);

}

double calc_veg_roughness(double height, double displacement, double Z0_SOIL, double L) {
/**********************************************************************
  calc_veg_roughness		Markus Schnorbus		August 9, 2015

  This subroutine estimates the roughness height of vegetation
  as a function of vegetation height and canopy density (via leaf area
  index), as suggested by Choudhury and Monteith (1988).
**********************************************************************/

  double value;
  double X;

  X = COEF_DRAG*L;
  if (X <= 0.2)
	  value = Z0_SOIL + 0.3 * height * pow(X,0.5);
  else
	  value = 0.3*height*(1-displacement/height);

  return (value);

}
