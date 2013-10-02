#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

void prepare_full_energy(int               iveg,
			 int               Nveg,
			 int               Nnodes,
			 dist_prcp_struct *prcp,
			 soil_con_struct  *soil_con, //TODO: make this const since it is read only
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

  layer = (layer_data_struct *)calloc(state->options.Nlayer,
				      sizeof(layer_data_struct));
  for (std::vector<HRUElement>::iterator it = prcp->hruElements.begin(); it != prcp->hruElements.end(); ++it) {

    // Only apply this function to the specified vegetation type.
    if (it->vegIndex != iveg) continue;

    int band = it->bandIndex;

    if (soil_con->AreaFract[band] > 0.0) {

      /* Compute average soil moisture values for distributed precipitation */

      for(int i=0;i<state->options.Nlayer;i++)
        layer[i] = find_average_layer(&(it->cell[WET].layer[i]),
				      &(it->cell[DRY].layer[i]),
				      soil_con->depth[i], prcp->mu[iveg], state);
    
      /* Compute top soil layer moisture content (mm/mm) */

      moist0[band] = layer[0].moist / ( soil_con->depth[0] * 1000. );

      /* Compute top soil layer ice content (mm/mm) */

      if(state->options.FROZEN_SOIL && soil_con->FS_ACTIVE){
        if((it->energy.T[0]
	    + it->energy.T[1])/2.<0.) {
	  ice0[band] = moist0[band] 
	    - maximum_unfrozen_water((it->energy.T[0]
				      + it->energy.T[1]) / 2.,
				     soil_con->porosity[0], soil_con->effective_porosity[0],
				     soil_con->max_moist[0]
				     / (soil_con->depth[0] * 1000.),
				     soil_con->bubble[0], soil_con->expt[0]);
	  if(ice0[band]<0.) ice0[band]=0.;
        }
        else ice0[band]=0.;
      }
      else {
        ice0[band] = 0.;
      }

      /** Compute Soil Thermal Properties **/
      compute_soil_layer_thermal_properties(layer,soil_con->depth,
					    soil_con->bulk_dens_min,
					    soil_con->soil_dens_min,
					    soil_con->quartz,
					    soil_con->bulk_density,
					    soil_con->soil_density,
					    soil_con->organic,
					    soil_con->frost_fract,
					    state->options.Nlayer);
    
      /** Save Thermal Conductivities for Energy Balance **/
      it->energy.kappa[0] = layer[0].kappa;
      it->energy.Cs[0]    = layer[0].Cs;
      it->energy.kappa[1] = layer[1].kappa;
      it->energy.Cs[1]    = layer[1].Cs;

    } else {
      ice0[band] = 0.;
    }

  }

  free((char *)layer);

}
