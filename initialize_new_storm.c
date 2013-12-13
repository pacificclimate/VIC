#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include <math.h>

static char vcid[] = "$Id$";

int  initialize_new_storm(HRU& hru,
			  int                 rec,
			  double              new_mu,
			  const ProgramState *state) {
/**********************************************************************
  initialize_new_storm.c	Keith Cherkauer		June 24, 1998

  This subroutine averages all soil moisture components before the
  start of a new storm, so that the model loses memory of the previous
  storm's location (which was never really known).

  Modifications:
  07-13-98 modified to redistribute Wdew within all defined
           elevation bands                                         KAC
  6-8-2000 modified to work with spatially distributed frost       KAC
  2007-Apr-04 Modified to return to main subroutine on cell error GCT/KAC
**********************************************************************/

  double old_mu = hru.mu;

  /** Redistribute Soil Moisture **/
  for (int layer = 0; layer < state->options.Nlayer; layer++) {

    hru_data_struct& cellWet = hru.cell[WET];
    hru_data_struct& cellDry = hru.cell[DRY];

    double temp_wet = cellWet.layer[layer].moist;
    double temp_dry = cellDry.layer[layer].moist;
    unsigned char error = average_moisture_for_storm(&temp_wet, &temp_dry,
        old_mu, new_mu);
    if (error) {
      fprintf(stderr,
          "moist does not balance before new storm: %f -> %f record %i\n",
          cellWet.layer[layer].moist * new_mu
              + cellDry.layer[layer].moist * (1. - new_mu), temp_wet + temp_dry,
          rec);
      return (ERROR);
    }
    cellWet.layer[layer].moist = temp_wet;
    cellDry.layer[layer].moist = temp_dry;

#if SPATIAL_FROST
    for (int frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
      temp_wet = cellWet.layer[layer].soil_ice[frost_area];
      temp_dry = cellDry.layer[layer].soil_ice[frost_area];
#else
    temp_wet = cellWet.layer[layer].soil_ice;
    temp_dry = cellDry.layer[layer].soil_ice;
#endif
    error = average_moisture_for_storm(&temp_wet, &temp_dry, old_mu, new_mu);
    if (error) {
#if SPATIAL_FROST
      fprintf(stderr,"ice does not balance before new storm: %f -> %f record %i\n",
          cellWet.layer[layer].soil_ice[frost_area] * new_mu
          + cellDry.layer[layer].soil_ice[frost_area]
          * (1. - new_mu), temp_wet + temp_dry, rec);
#else
      fprintf(stderr,
          "ice does not balance before new storm: %f -> %f record %i\n",
          cellWet.layer[layer].soil_ice * new_mu
              + cellDry.layer[layer].soil_ice * (1. - new_mu),
          temp_wet + temp_dry, rec);
#endif
      return (ERROR);
    }
#if SPATIAL_FROST
    cellWet.layer[layer].soil_ice[frost_area] = temp_wet;
    cellDry.layer[layer].soil_ice[frost_area] = temp_dry;
  }
#else
    cellWet.layer[layer].soil_ice = temp_wet;
    cellDry.layer[layer].soil_ice = temp_dry;
#endif 
  }

  /****************************************
   Redistribute Stored Water in Vegetation
   ****************************************/
  if (hru.isArtificialBareSoil == false) {
    double temp_wet = hru.veg_var[WET].Wdew;
    double temp_dry = hru.veg_var[DRY].Wdew;
    unsigned char error = average_moisture_for_storm(&temp_wet, &temp_dry, old_mu, new_mu);
    if (error) {
      fprintf(stderr,"Wdew does not balance before new storm: %f -> %f record %i\n",
          hru.veg_var[WET].Wdew * new_mu + hru.veg_var[DRY].Wdew * (1. - new_mu), temp_wet + temp_dry, rec);
      return (ERROR);
    }
    hru.veg_var[WET].Wdew = temp_wet;
    hru.veg_var[DRY].Wdew = temp_dry;
  }
  return (0);
}

unsigned char average_moisture_for_storm(double *wet_value,
					 double *dry_value,
					 double  old_mu,
					 double  new_mu) {
/**********************************************************************
  This subroutine averages total soil moisture between the wet and dry
  fractions of the soil column
**********************************************************************/

  unsigned char error;
  double temp_wet;
  double temp_dry;
  double diff;

  temp_wet = *wet_value * old_mu;
  temp_dry = *dry_value * (1. - old_mu);
  *wet_value = temp_wet + temp_dry;
  *dry_value = temp_wet + temp_dry;
  diff = (temp_wet+temp_dry) - (*wet_value*new_mu + *dry_value*(1.-new_mu));
  if(fabs(diff) > 1.e-10) error = TRUE;
  else error = FALSE;

  return (error);
}
