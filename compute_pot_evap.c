#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vicNl.h"

static char vcid[] = "$Id$";

void compute_pot_evap(int veg_class, 
		      const dmy_struct *dmy, 
		      int rec, 
		      int dt, 
		      double shortwave,
		      double net_longwave,
		      double tair, 
		      double vpd,
		      double elevation,
		      double **aero_resist,
		      double *pot_evap,
		      const ProgramState* state)
/****************************************************************************
                                                                           
  compute_pot_evap: computes potential evaporation for several different
                    reference land cover types (which are defined in
                    vicNl_def.h and global.h).

  modifications:
  2009-Aug-21 Fixed bug in assignment of albedo for natural vegetation.	TJB

****************************************************************************/
{
  int NVegLibTypes;
  int i;
  double albedo;
  double net_short;
  double net_rad;
  double rs;
  double rarc;
  double RGL;
  double lai;
  double gsm_inv;
  char ref_crop;
  double rc;
  double ra;

  /************************************************
  Estimate and store potential evap estimates using penman equation
  ************************************************/

  NVegLibTypes = state->veg_lib[0].NVegLibTypes;
  for (i=0; i<N_PET_TYPES; i++) {
    if (i < N_PET_TYPES_NON_NAT) {
      rs = state->veg_lib[NVegLibTypes+i].rmin;
      rarc = state->veg_lib[NVegLibTypes+i].rarc;
      RGL = state->veg_lib[NVegLibTypes+i].RGL;
      lai = state->veg_lib[NVegLibTypes+i].LAI[dmy[rec].month-1];
      albedo = state->veg_lib[NVegLibTypes+i].albedo[dmy[rec].month-1];
    }
    else {
      rs = state->veg_lib[veg_class].rmin;
      if (i == PET_VEGNOCR) rs = 0;
      rarc = state->veg_lib[veg_class].rarc;
      RGL = state->veg_lib[veg_class].RGL;
      lai = state->veg_lib[veg_class].LAI[dmy[rec].month-1];
      albedo = state->veg_lib[veg_class].albedo[dmy[rec].month-1];
    }
    gsm_inv = 1.0;
    ref_crop = ref_veg_ref_crop[i];
    rc = calc_rc(rs, net_short, RGL, tair, vpd, lai, gsm_inv, ref_crop);
    if (i < N_PET_TYPES_NON_NAT || !state->veg_lib[veg_class].overstory)
      ra = aero_resist[i][0];
    else
      ra = aero_resist[i][1];
    net_short = (1.0 - albedo) * shortwave;
    net_rad = net_short + net_longwave;
    pot_evap[i] = penman(tair, elevation, net_rad, vpd, ra, rc, rarc) * dt/24.0;
  }

}
