#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

dist_prcp_struct make_dist_prcp(int  nveg, const int NUM_LAYERS, const int NUM_SNOW_BAND)
/**********************************************************************
	read_dist_prcp	Keith Cherkauer		May 21, 1996

  This routine creates an array of structures which will store 
  necessary information about the distribution of precipitation, moisture,
  evaporation, and dew.  Mu represents the fractional area of the grid 
  that receives precipitation (wet), while 1-mu is the corresponding 
  area that receives no precipitation.  The value of mu changes with
  the intensity of incoming precipitation, and is set in the routine
  dist_prec.

  modifications:
  11-18-02 Modified to allocate vegetation variables for the 
           wetland vegetation class.                             LCB
  01-Nov-04 Updated arglist to make_energy_bal() as part of fix for
	    QUICK_FLUX state file compatibility.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.  TJB
  2009-Jul-31 Removed extra lake/wetland tile.			TJB
**********************************************************************/
{
  dist_prcp_struct temp(NUM_SNOW_BAND);
  int              Nitems;

  Nitems = nveg + 1;

  temp.mu     = (double *)calloc(Nitems,sizeof(double));
  for (int i = 0; i < Nitems; i++ ) temp.mu[i] = 1;

  temp.hruElements.reserve(Nitems * NUM_SNOW_BAND);
  for (int i = 0; i < Nitems; i++) {
    for (int j = 0; j < NUM_SNOW_BAND; j++) {
      HRUElement e;
      e.vegIndex = i;
      e.bandIndex = j;
      e.energy.frozen = FALSE;
      temp.hruElements.push_back(e);
    }
  }

  return (temp);

}
