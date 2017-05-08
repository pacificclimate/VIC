#include <stdlib.h>
#include <stdio.h>
#include "vicNl.h"

static char vcid[] = "$Id$";

double get_dist(double lat1, double long1, double lat2, double long2)
/*******************************************************************************
  Function: double get_dist(double lat1, double long1, double lat2, double long2)
  Returns : distance between two locations

  Modifications:
  2007-Nov-06 Moved to separate file from read_lakeparam.c.		TJB
  2017-May-03 Updated to use geocentric earth radius and
              haversine great-circle distance formula           MAS
********************************************************************************/
{
  double theta1;
  double phi1;
  double theta2;
  double phi2;
  double dtor;
  double term1;
  double term2;
  double term3;
  double RAD;
  double distance;

  dtor = 2.0*PI/360.0;
  theta1 = dtor*long1;
  phi1 = dtor*lat1;
  theta2 = dtor*long2;
  phi2 = dtor*lat2;
  
  //Geocentric Earth radius (https://en.wikipedia.org/wiki/Earth_radius)
  term1 = pow(pow(RADIUS_E, 2.0)*cos(phi1), 2.0);
  term2 = pow(pow(RADIUS_P, 2.0)*sin(phi1), 2.0);
  term3 = (pow(RADIUS_E*cos(phi1),2.0) + pow(RADIUS_P*sin(phi), 2.0));
  RAD = sqrt((term1 + term2)/term3);
  //Distance using haversine formula(https://en.wikipedia.org/wiki/Haversine_formula)
  distance = 2.0*RAD*(asin(sqrt(pow(sin((phi2-phi1)/2), 2.0) + cos(phi1)*cos(phi2)*pow(sin((theta2-theta1)/2), 2.0))));
  
  return distance;
}  

