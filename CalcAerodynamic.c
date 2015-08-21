/*
 * SUMMARY:      CalcAerodynamic.c - Calculate the aerodynamic resistances
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu, pstorck@u.washington.edu
 * ORIG-DATE:    Thu Mar 27 18:00:10 1997
 * LAST-MOD: Thu Mar  8 13:24:10 2001 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Calculate the aerodynamic resistances
 * DESCRIP-END.
 * FUNCTIONS:    CalcAerodynamic()
 * COMMENTS:     Modified for use with the vicNl model 3-12-98
 *		 by Keith Cherkauer
 */

#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include <math.h>

static char vcid[] = "$Id$";

  
/*****************************************************************************
  Function name: CalcAerodynamic()

  Purpose      : Calculate the aerodynamic resistance for each vegetation 
                 layer, and the wind 2m above the layer boundary.  In case of 
                 an overstory, also calculate the wind in the overstory.
                 The values are normalized based on a reference height wind 
                 speed, Uref, of 1 m/s.  To get wind speeds and aerodynamic 
                 resistances for other values of Uref, you need to multiply 
                 the here calculated wind speeds by Uref and divide the 
                 here calculated aerodynamic resistances by Uref
                 
  Required     :
    int NVegLayers - Number of vegetation layers
    char OverStory - flag for presence of overstory.  Only used if NVegLayers 
                     is equal to 1
    double Zref[0]     - Reference height for windspeed
    double n        - Attenuation coefficient for wind in the overstory
    double Height  - Height of the vegetation layers (top layer first)
    double Trunk    - Multiplier for Height[0] that indicates the top of the
                     trunk space

    VegConditions &aero_resist    contains aerodynamic resistance values for each vegetation condition
    VegConditions &wind_speed     contains wind speeds for each vegetation condition
    VegConditions &displacement   contains displacement height values for each vegetation condition
    VegConditions &ref_height     contains reference height values for each vegetation condition
    VegConditions &roughness      contains roughness length values for each vegetation condition

  Returns      : int

  Modifies     :
    VegConditions &aero_resist
    VegConditions &wind_speed
    VegConditions &displacement
    VegConditions &ref_height
    VegConditions &roughness
   
  Comments     :
*****************************************************************************/
int  CalcAerodynamic(char    OverStory,     /* overstory flag */
                     double  Height,        /* vegetation height */
                     double  Trunk,         /* trunk ratio parameter */
                     double  Z0_SNOW,       /* snow roughness */
                     double  Z0_SOIL,       /* soil roughness */
                     double  n,             /* wind attenuation parameter */
                     VegConditions& aero_resist,   /* aerodynamic resistances */
                     VegConditions& wind_speed,    /* adjusted wind speed */
                     VegConditions& displacement,  /* vegetation displacement */
                     VegConditions& ref_height,    /* vegetation reference height */
                     VegConditions& roughness)     /* vegetation roughness */
{
  /******************************************************************
  Modifications:
  2007-Apr-04 Modified to catch and return error flags from surface_fluxes
              subroutine.                                      GCT/KAC
  2009-Jun-09 Modified to use extension of veg_lib structure to contain
	      bare soil information.				TJB
  *******************************************************************/


  double d_Lower;
  double d_Upper;
  double K2;
  double Uh;
  double Ut;
  double Uw;
  double Z0_Lower;
  double Z0_Upper;
  double Zt;
  double Zw;
  double tmp_wind;

  tmp_wind = wind_speed.snowFree;

  K2 = von_K * von_K;
  
  /* No OverStory, thus maximum one soil layer */
  
  if ( OverStory == FALSE ) {
    
    /* vegetation cover */
    Z0_Lower = roughness.snowFree;
    d_Lower  = displacement.snowFree;
    
    /* No snow */
    wind_speed.snowFree  = log((2. + Z0_Lower)/Z0_Lower)/log((ref_height.snowFree - d_Lower)/Z0_Lower);
    /****** DHSVM ******/
    aero_resist.snowFree = log((2. + Z0_Lower)/Z0_Lower) * log((ref_height.snowFree - d_Lower)/Z0_Lower)
            /K2;
    /***** Old VIC *****
    aero_resist.snowFree = log((2. + (1.0/0.63 - 1.0) * d_Lower) / Z0_Lower)
          * log((2. + (1.0/0.63 - 1.0) * d_Lower) / (0.1*Z0_Lower)) / K2;
    ******************/

    /* Set aerodynamic resistance terms for canopy  - not used */
    ref_height.canopyIfOverstory   = ref_height.snowFree;
    roughness.canopyIfOverstory    = roughness.snowFree;
    displacement.canopyIfOverstory = displacement.snowFree;
    wind_speed.canopyIfOverstory = wind_speed.snowFree;
    aero_resist.canopyIfOverstory = aero_resist.snowFree;

    /* Set aerodynamic resistance terms for snow */
    ref_height.snowCovered = ref_height.snowFree;
    roughness.snowCovered = Z0_SNOW;
    displacement.snowCovered = 0.;
    wind_speed.snowCovered = log((2. + Z0_SNOW)/Z0_SNOW)/log(ref_height.snowCovered/Z0_SNOW);
    aero_resist.snowCovered = log((2. + Z0_SNOW)/Z0_SNOW) * log(ref_height.snowCovered/Z0_SNOW)/K2;
    ref_height.snowCovered = 2. + Z0_SNOW;

    /* Set aerodynamic resistance terms for Glacier */
    ref_height.glacierSurface = ref_height.snowFree;
    roughness.glacierSurface = Z0_Lower;
    displacement.glacierSurface = 0.;
    wind_speed.glacierSurface =  log((2. + Z0_Lower)/Z0_Lower)/log(ref_height.glacierSurface / Z0_Lower);
    aero_resist.glacierSurface = log((2. + Z0_Lower)/Z0_Lower)*log(ref_height.glacierSurface / Z0_Lower) / K2;
    ref_height.glacierSurface = 2. + Z0_Lower;

  }
  
  /* Overstory present, one or two vegetation layers possible */
  else {
    Z0_Upper = roughness.snowFree;
    d_Upper  = displacement.snowFree;
    
    Z0_Lower = Z0_SOIL;
    d_Lower  = 0;
    
    Zw = 1.5 * Height - 0.5 * d_Upper;
    Zt = Trunk * Height;

    if (Zt < (Z0_Lower+d_Lower)) {
      fprintf(stderr,"ERROR: CalcAerodynamic - Trunk space height below \"center\" of lower boundary");
      return( ERROR );
    }

    /* Resistance for overstory */
    aero_resist.canopyIfOverstory = log((ref_height.snowFree-d_Upper)/Z0_Upper)/K2
          * (Height/(n*(Zw-d_Upper)) 
          * (exp(n*(1-(d_Upper+Z0_Upper)/Height))-1)
          + (Zw-Height)/(Zw-d_Upper)
          + log((ref_height.snowFree-d_Upper)/(Zw-d_Upper)));
    
    /* Wind at different levels in the profile */
    Uw = log((Zw-d_Upper)/Z0_Upper) / log((ref_height.snowFree-d_Upper)/Z0_Upper);
    Uh = Uw - (1-(Height-d_Upper)/(Zw-d_Upper))
       / log((ref_height.snowFree-d_Upper)/Z0_Upper);
    wind_speed.canopyIfOverstory = Uh * exp(n * ((Z0_Upper+d_Upper)/Height - 1.));
    Ut = Uh * exp(n * (Zt/Height - 1.));
    
    /* resistance at the lower boundary */
    wind_speed.snowFree  = Ut*log((2. + Z0_Lower)/Z0_Lower)/log(Zt/Z0_Lower);
    aero_resist.snowFree = log((2. + Z0_Lower)/Z0_Lower) * log(Zt/Z0_Lower) / (K2*Ut);
    
    /***** Old VIC *****
    wind_speed.snowFree  = log((2. + Z0_Upper)/Z0_Upper)/log((ref_height.snowFree - d_Upper)/Z0_Upper);
    aero_resist.snowFree = log((2. + (1.0/0.63 - 1.0) * d_Upper) / Z0_Upper)
          * log((2. + (1.0/0.63 - 1.0) * d_Upper) / (0.1*Z0_Upper)) / K2;
    ******************/



    /* Snow */
    /* case 1: the wind profile to a height of 2m above the lower boundary is 
       entirely logarithmic */
    if (Zt > (2. + Z0_SNOW)) {
      wind_speed.snowCovered = Ut*log((2.+Z0_SNOW)/Z0_SNOW)/log(Zt/Z0_SNOW);
      aero_resist.snowCovered = log((2.+Z0_SNOW)/Z0_SNOW) * log(Zt/Z0_SNOW)/(K2*Ut);
    }
    
    /* case 2: the wind profile to a height of 2m above the lower boundary 
       is part logarithmic and part exponential, but the top of the overstory 
       is more than 2 m above the lower boundary */
    else if (Height > (2. + Z0_SNOW)) {
      wind_speed.snowCovered = Uh * exp(n * ((2. + Z0_SNOW)/Height - 1.));
      aero_resist.snowCovered = log(Zt/Z0_SNOW) * log(Zt/Z0_SNOW)/
        (K2*Ut) +
        Height * log((ref_height.snowFree-d_Upper)/Z0_Upper) / (n*K2*(Zw-d_Upper)) *
        (exp(n*(1-Zt/Height)) - exp(n*(1-(Z0_SNOW+2.)/Height)));
    }
    
    /* case 3: the top of the overstory is less than 2 m above the lower 
       boundary.  The wind profile above the lower boundary is part 
       logarithmic and part exponential, but only extends to the top of the 
       overstory */
    else {
      wind_speed.snowCovered = Uh;
      aero_resist.snowCovered = log(Zt/Z0_SNOW) * log(Zt/Z0_SNOW)/
        (K2*Ut) +
        Height * log((ref_height.snowFree-d_Upper)/Z0_Upper) / (n*K2*(Zw-d_Upper)) *
        (exp(n*(1-Zt/Height)) - 1);
      fprintf(stderr, "WARNING:  Top of overstory is less than 2 meters above the lower boundary\n");
    }

    /** Set aerodynamic resistance terms for canopy */
    /* not currently used */
    ref_height.canopyIfOverstory   = ref_height.snowFree;
    roughness.canopyIfOverstory    = roughness.snowFree;
    displacement.canopyIfOverstory = displacement.snowFree;
    ref_height.snowFree   = 2. + Z0_Lower;
    roughness.snowFree    = Z0_Lower;
    displacement.snowFree = d_Lower;

    /** Set aerodynamic resistance terms for snow */
    ref_height.snowCovered   = 2. + Z0_SNOW;
    roughness.snowCovered    = Z0_SNOW;
    displacement.snowCovered = 0.;

    /** Set aerodynamic resistance terms for glacier - not used */
    ref_height.glacierSurface = 2. + Z0_Lower;
    roughness.glacierSurface = Z0_Lower;
    displacement.glacierSurface = 0.;
    //aero_resist.glacierSurface = HUGE_RESIST;
    //wind_speed.glacierSurface remains INVALID

  }

  if ( tmp_wind > 0. ) {
    wind_speed.snowFree *= tmp_wind;
    aero_resist.snowFree /= tmp_wind;
    if(IS_VALID(wind_speed.canopyIfOverstory)) {
      wind_speed.canopyIfOverstory *= tmp_wind;
      aero_resist.canopyIfOverstory /= tmp_wind;
    }
    if(IS_VALID(wind_speed.snowCovered)) {
      wind_speed.snowCovered *= tmp_wind;
      aero_resist.snowCovered /= tmp_wind;
    }
    if(IS_VALID(wind_speed.glacierSurface)) {
      wind_speed.glacierSurface *= tmp_wind;
      aero_resist.glacierSurface /= tmp_wind;
    }
  }
  else {
    wind_speed.snowFree *= tmp_wind;
    aero_resist.snowFree = HUGE_RESIST;
    if(IS_VALID(wind_speed.canopyIfOverstory))
      wind_speed.canopyIfOverstory *= tmp_wind;
    aero_resist.canopyIfOverstory = HUGE_RESIST;
    if(IS_VALID(wind_speed.snowCovered))
      wind_speed.snowCovered *= tmp_wind;
    aero_resist.snowCovered = HUGE_RESIST;
    if(IS_VALID(wind_speed.glacierSurface))
      wind_speed.glacierSurface *= tmp_wind;
    aero_resist.glacierSurface = HUGE_RESIST;
  }
  return (0);

}
