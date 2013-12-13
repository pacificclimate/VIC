#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include <math.h>

static char vcid[] = "$Id$";

void prepare_full_energy(HRU& hru,
			 int               Nnodes,
			 const soil_con_struct  *soil_con,
			 double           *moist0,
			 double           *ice0,
			 const ProgramState* state) {
/*******************************************************************
  prepare_full_energy.c      Keith Cherkauer       January 20, 2000

  This subroutine returns the soil thermal properties, moisture 
  and ice contents for the top two layers for use with the QUICK_FLUX
  ground heat flux solution.

  Modifications:
  01-20-00 split into separate file, formerly at the end of 
           full_energy.c                                      KAC
  03-12-03 modified so that ice content is set to zero unless
           the frozen soil algorithm is implemented and active
           in the current grid cell.                          KAC
  2007-Aug-09 Added features for EXCESS_ICE option.			JCA
  2008-Jan-23 Changed ice0 from a scalar to an array.  Previously,
	      when state->options.SNOW_BAND > 1, the value of ice0 computed
	      for earlier bands was always overwritten by the value
	      of ice0 computed for the final band (even if the final
	      band had 0 area).						JS via TJB
  2008-May-05 Changed moist from a scalar to an array (moist0).  Previously,
	      when state->options.SNOW_BAND > 1, the value of moist computed
	      for earlier bands was always overwritten by the value
	      of moist computed for the final band (even if the final
	      band had 0 area).						KAC via TJB
  2011-Jun-03 Added state->options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.				TJB

*******************************************************************/

  layer_data_struct *layer;

  layer = (layer_data_struct *) calloc(state->options.Nlayer,
      sizeof(layer_data_struct));

  int band = hru.bandIndex;

  if (soil_con->AreaFract[band] > 0.0) {

    /* Compute average soil moisture values for distributed precipitation */

    for (int i = 0; i < state->options.Nlayer; i++)
      layer[i] = find_average_layer(&(hru.cell[WET].layer[i]),
          &(hru.cell[DRY].layer[i]), soil_con->depth[i], hru.mu, state);

    /* Compute top soil layer moisture content (mm/mm) */

    moist0[band] = layer[0].moist / (soil_con->depth[0] * 1000.);

    /* Compute top soil layer ice content (mm/mm) */

    if (state->options.FROZEN_SOIL && soil_con->FS_ACTIVE) {
      if ((hru.energy.T[0] + hru.energy.T[1]) / 2. < 0.) {
        ice0[band] = moist0[band]
            - maximum_unfrozen_water((hru.energy.T[0] + hru.energy.T[1]) / 2.,
                soil_con->max_moist[0] / (soil_con->depth[0] * 1000.),
                soil_con->bubble[0], soil_con->expt[0]);
        if (ice0[band] < 0.)
          ice0[band] = 0.;
      } else
        ice0[band] = 0.;
    } else {
      ice0[band] = 0.;
    }

    /** Compute Soil Thermal Properties **/
    compute_soil_layer_thermal_properties(layer, soil_con,
        state->options.Nlayer);

    /** Save Thermal Conductivities for Energy Balance **/
    hru.energy.kappa[0] = layer[0].kappa;
    hru.energy.Cs[0] = layer[0].Cs;
    hru.energy.kappa[1] = layer[1].kappa;
    hru.energy.Cs[1] = layer[1].Cs;

  } else {
    ice0[band] = 0.;
  }

  free((char *)layer);

}
