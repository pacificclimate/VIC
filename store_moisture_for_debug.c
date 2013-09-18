#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

#if LINK_DEBUG
void store_moisture_for_debug(int                 iveg,
		              int                 Nveg,
			      double             *mu,
			      cell_data_struct ***cell,
			      veg_var_struct   ***veg_var,
			      snow_data_struct  **snow,
			      soil_con_struct    *soil_con,
			      const ProgramState* state) {
/****************************************************************
  This subroutine was written to save the current water storage
  terms for use in calculating the model water balance error
****************************************************************/

  int               Ndist;
  int               i;
  int               band;
  int               dist;
  int               Nbands;

  if(state->options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
  Nbands = state->options.SNOW_BAND;

  for(band=0;band<Nbands;band++) {
    if(soil_con->AreaFract[band]>0) {
      for(dist=0;dist<Ndist;dist++) 
	for(i=0;i<state->options.Nlayer+3;i++)
	  state->debug.store_moist[dist][band][i] = 0.;
      if(iveg<Nveg) {
	for(dist=0;dist<Ndist;dist++)
	  state->debug.store_moist[dist][band][0]
	    += veg_var[dist][iveg][band].Wdew;
	state->debug.store_moist[WET][band][0] += (snow[iveg][band].snow_canopy)
	  * 1000.;
      }
      for(dist=0;dist<Ndist;dist++)
	state->debug.store_moist[dist][band][state->options.Nlayer+2]
	  += state->debug.store_moist[dist][band][0];
      state->debug.store_moist[WET][band][1] += (snow[iveg][band].swq*1000.);
      for(dist=0;dist<Ndist;dist++)
	state->debug.store_moist[dist][band][state->options.Nlayer+2]
	  += state->debug.store_moist[dist][band][1];
      for(i=0;i<state->options.Nlayer;i++) {
	for(dist=0;dist<Ndist;dist++) {
	  state->debug.store_moist[dist][band][i+2]
	    = cell[dist][iveg][band].layer[i].moist;
	  state->debug.store_moist[dist][band][state->options.Nlayer+2]
	    += state->debug.store_moist[dist][band][i+2];
	}
      }
    }
  }
}
#endif

