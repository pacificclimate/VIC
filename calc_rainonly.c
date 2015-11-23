#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vicNl.h"

#define MIN_PREC     1.e-5      /* smallest amount of precipitation that
				   	   	   	   	   is allowed to fall as snow or rain in
				   	   	   	   	   a mixed precipitation event */

static char vcid[] = "$Id$";

double calc_rainonly(double air_temp,
		     double prec,
		     double MAX_SNOW_TEMP,
		     double MIN_RAIN_TEMP,
		     double precipitation_mu,
		     const ProgramState* state) {
/**********************************************************************
  calc_rainonly.c	Keith Cherkauer		March 7, 1998

  Determines from the air temperature what fraction of incoming
  precipitation is frozen and unfrozen (snow and rain).

  Modifications:
  09-22-98 Modified to filter out very small fractions of snow
           or rain in mixed precipitation.  Minimum value MIN_PREC
	   is adjusted to account for the size of mu (minimum
	   is based of fractional precipitation with mu=1, since
	   snow cannot be solved for when mu<1).                  KAC
  10-May-04 Changed test
		if ( MAX_SNOW_TEMP < MIN_RAIN_TEMP )
	    to
		if ( MAX_SNOW_TEMP <= MIN_RAIN_TEMP )
            to avoid possibility of dividing by zero.		TJB
  10-May-04 Changed test
		else if(air_temp > MAX_SNOW_TEMP)
	    to
		else if(air_temp >= MAX_SNOW_TEMP)
	    to fix situation in which, if air_temp = MAX_SNOW_TEMP,
	    rainfall (rainonly) was set to 0 and snowfall was set
	    to 100% of precip, causing function to fail.	TJB
  2007-Apr-04 Modified to handle grid cell errors by returning to the
           main subroutine, rather than ending the simulation.   GCT/KAC
  2015-Nov-19 Modified to handle two precipitation partitioning algorithms
  	    algorithms based on value of TEMP_TH_TYPE parameter. Choices are:
  	    VIC_412 - original VIC v4.1.2 linear partitioning between MAX_SNOW_TEMP
  	    		  and MIN_RAIN_TEMP;
  	    KIENZLE - S-shaped curve with two parameters, the temperature at which
  	    		  50% of precipitation occurs as snow (TT) and the temperature
  	    		  range over which mixed precipitation occurs (TR). This
  	    		  algorithm follows Kienzle (2008).  MAS
**********************************************************************/

  double rainonly = 0;

  if(state->options.TEMP_TH_TYPE == VIC_412){

    if ( MAX_SNOW_TEMP <= MIN_RAIN_TEMP ) {
      fprintf( stderr, "ERROR: For method VIC_412, MAX_SNOW_TEMP must be greater than MIN_RAIN_TEMP");
      return (ERROR);
    }
    if(air_temp < MAX_SNOW_TEMP && air_temp > MIN_RAIN_TEMP) {
      rainonly = (air_temp - MIN_RAIN_TEMP)
          / (MAX_SNOW_TEMP - MIN_RAIN_TEMP) * prec;
    }
    else if(air_temp >= MAX_SNOW_TEMP) {
      rainonly = prec;
    }
  }

  else if(state->options.TEMP_TH_TYPE == KIENZLE){

	if(MIN_RAIN_TEMP <= 0.){
	  fprintf( stderr, "ERROR: For method KIENZLE, temperature range (MIN_RAIN_TEMP) must be greater than 0");
	  return(ERROR);
	}

	double rfrac = 0.;
	double TT = MAX_SNOW_TEMP;
	double TR = MIN_RAIN_TEMP;
	double D = 1.4 * TR;
	double E1 = 5. *   pow((air_temp-TT)/D, 3.0);
	double E2 = 6.76 * pow((air_temp-TT)/D, 2.0);
	double E3 = 3.19 * (air_temp-TT)/D;

	if(air_temp <= TT){
		rfrac = E1 + E2 + E3 + 0.5;
	} else {
		rfrac = E1 - E2 + E3 + 0.5;
	}
	/* Ensure rain fraction constrained between 1 and 0 */
	if(rfrac < 0.) rfrac = 0.;
	if(rfrac > 1.) rfrac = 1.;

	rainonly = rfrac * prec;
  }

  if(rainonly < MIN_PREC) rainonly = 0.;
  if((prec-rainonly) < MIN_PREC) rainonly = prec;
  if(precipitation_mu < 1) rainonly = prec;

  return(rainonly);

}

#undef MIN_PREC
