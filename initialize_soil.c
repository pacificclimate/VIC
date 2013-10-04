#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void initialize_soil (std::vector<HRU>& elements,
                      int                dist,
                      soil_con_struct   *soil_con,
                      veg_con_struct    *veg_con,
                      int                veg_num,
                      const ProgramState *state)
/**********************************************************************
	initialize_soil		Keith Cherkauer		July 31, 1996

  This routine initializes the soil variable arrays for each new
  grid cell.

  modifications:
  11-18-02 Modified to initialize wetland soil moisture.          LCB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2007-Aug-10 Added features for EXCESS_ICE option.			JCA
  2009-Mar-16 Modified to use min_liq (minimum allowable liquid water
	      content) instead of resid_moist.  For unfrozen soil,
	      min_liq = resid_moist.					TJB
  2009-Jul-31 Replaced extra lake/wetland veg tile with reference to
	      veg_con[j].LAKE.						TJB
  2009-Dec-11 Removed min_liq and options.MIN_LIQ.			TJB
  2011-Mar-01 Now initializes more cell data structure terms, including
	      asat and zwt.						TJB
**********************************************************************/
{
  int lindex, frost_area;
  double tmp_moist[MAX_LAYERS];
  double tmp_runoff;
  
  for (std::vector<HRU>::iterator it = elements.begin();
      it != elements.end(); ++it) {

    cell_data_struct& cellRef = it->cell[dist];

    cellRef.baseflow = 0;
    cellRef.runoff = 0;
    for (lindex = 0; lindex < state->options.Nlayer; lindex++) {
      cellRef.layer[lindex].evap = 0;
      cellRef.layer[lindex].moist = soil_con->init_moist[lindex];
      if (cellRef.layer[lindex].moist > soil_con->max_moist[lindex])
        cellRef.layer[lindex].moist = soil_con->max_moist[lindex];
      tmp_moist[lindex] = cellRef.layer[lindex].moist;
#if SPATIAL_FROST
      for (frost_area=0; frost_area<FROST_SUBAREAS; frost_area++) {
        cellRef.layer[lindex].soil_ice[frost_area] = 0;
      }
#else
      cellRef.layer[lindex].soil_ice = 0;
#endif
    }
    compute_runoff_and_asat(soil_con, tmp_moist, 0, &(cellRef.asat),
        &tmp_runoff, state);
    wrap_compute_zwt(soil_con, &(cellRef), state);
  }
}
