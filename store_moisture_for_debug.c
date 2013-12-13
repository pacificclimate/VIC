#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include <math.h>

static char vcid[] = "$Id$";

#if LINK_DEBUG
void store_moisture_for_debug(const HRU&         hru,
                              const soil_con_struct     *soil_con,
                              const ProgramState*       state) {
/****************************************************************
  This subroutine was written to save the current water storage
  terms for use in calculating the model water balance error
****************************************************************/

  int               Ndist;
  int               Nbands;

  if (state->options.DIST_PRCP)
    Ndist = 2;
  else
    Ndist = 1;
  Nbands = state->options.SNOW_BAND;

  int band = hru.bandIndex;

  if (soil_con->AreaFract[band] > 0) {
    for (int dist = 0; dist < Ndist; dist++)
      for (int i = 0; i < state->options.Nlayer + 3; i++)
        state->debug.store_moist[dist][band][i] = 0.;
    if (hru.isArtificialBareSoil == false) {
      for (int dist = 0; dist < Ndist; dist++)
        state->debug.store_moist[dist][band][0] += hru.veg_var[dist].Wdew;
      state->debug.store_moist[WET][band][0] += (hru.snow.snow_canopy) * 1000.;
    }
    for (int dist = 0; dist < Ndist; dist++)
      state->debug.store_moist[dist][band][state->options.Nlayer + 2] +=
          state->debug.store_moist[dist][band][0];
    state->debug.store_moist[WET][band][1] += (hru.snow.swq * 1000.);
    for (int dist = 0; dist < Ndist; dist++)
      state->debug.store_moist[dist][band][state->options.Nlayer + 2] +=
          state->debug.store_moist[dist][band][1];
    for (int i = 0; i < state->options.Nlayer; i++) {
      for (int dist = 0; dist < Ndist; dist++) {
        state->debug.store_moist[dist][band][i + 2] =
            hru.cell[dist].layer[i].moist;
        state->debug.store_moist[dist][band][state->options.Nlayer + 2] +=
            state->debug.store_moist[dist][band][i + 2];
      }
    }
  }
}
#endif

