#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

#if LINK_DEBUG
void store_moisture_for_debug(int                       iveg,
		                          int                       Nveg,
                              double                    *mu,
                              std::vector<HRUElement>&  hruElements,
                              const soil_con_struct           *soil_con,
                              const ProgramState*       state) {
/****************************************************************
  This subroutine was written to save the current water storage
  terms for use in calculating the model water balance error
****************************************************************/

  int               Ndist;
  int               i;
  int               dist;
  int               Nbands;

  if (state->options.DIST_PRCP)
    Ndist = 2;
  else
    Ndist = 1;
  Nbands = state->options.SNOW_BAND;

  for (std::vector<HRUElement>::iterator it = hruElements.begin(); it != hruElements.end(); ++it) {

    if (it->vegIndex != iveg) continue;
    int band = it->bandIndex;

    if (soil_con->AreaFract[band] > 0) {
      for (dist = 0; dist < Ndist; dist++)
        for (i = 0; i < state->options.Nlayer + 3; i++)
          state->debug.store_moist[dist][band][i] = 0.;
      if (iveg < Nveg) {
        for (dist = 0; dist < Ndist; dist++)
          state->debug.store_moist[dist][band][0] +=
              it->veg_var[dist].Wdew;
        state->debug.store_moist[WET][band][0] += (it->snow.snow_canopy)
            * 1000.;
      }
      for (dist = 0; dist < Ndist; dist++)
        state->debug.store_moist[dist][band][state->options.Nlayer + 2] +=
            state->debug.store_moist[dist][band][0];
      state->debug.store_moist[WET][band][1] += (it->snow.swq * 1000.);
      for (dist = 0; dist < Ndist; dist++)
        state->debug.store_moist[dist][band][state->options.Nlayer + 2] +=
            state->debug.store_moist[dist][band][1];
      for (i = 0; i < state->options.Nlayer; i++) {
        for (dist = 0; dist < Ndist; dist++) {
          state->debug.store_moist[dist][band][i + 2] =
              it->cell[dist].layer[i].moist;
          state->debug.store_moist[dist][band][state->options.Nlayer + 2] +=
              state->debug.store_moist[dist][band][i + 2];
        }
      }
    }
  }
}
#endif

