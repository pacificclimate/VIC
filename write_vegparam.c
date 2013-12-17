#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

static char vcid[] = "$Id$";

void write_vegparam(const cell_info_struct& cell, const ProgramState* state)
/**********************************************************************
	write_vegparam		Dag Lohmann	January 1996

  This routine writes vegetation parameters to stdout, used primarily 
  for debugging, and making sure the model is reading the proper 
  parameters..

  Modifications:
  5/21/96	Routine was modified to allow for variable
		number of layers				KAC
  4-12-98  Updated for new standard vegetation parameters       KAC

**********************************************************************/
{
  int i, j, l;
  int vegclass;

  printf("Vegetation Parameters:\n");
  printf("\tnumber of HRUs = %d\n",  (int)cell.prcp.hruList.size());

  for (std::vector<HRU>::const_iterator hru = cell.prcp.hruList.begin(); hru != cell.prcp.hruList.end(); ++hru) {
    vegclass = hru->veg_con.vegIndex;
    printf("\n\tveg_class            = %d\n",  state->veg_lib[vegclass].veg_class);
    printf("\tCv                   = %f\n", hru->veg_con.Cv);
    if(state->veg_lib[vegclass].overstory)
      printf("\tOverstory            = TRUE\n");
    else 
      printf("\tOverstory            = FALSE\n");
    printf("\trarc                 = %f s/m\n", state->veg_lib[vegclass].rarc);
    printf("\trmin                 = %f s/m\n", state->veg_lib[vegclass].rmin);
    for(l=0;l<state->options.ROOT_ZONES;l++)
      printf("\tzone_depth _fract%d   = %f %f\n",l+1,
          hru->veg_con.zone_depth[l],hru->veg_con.zone_fract[l]);
    for(l=0;l<state->options.Nlayer;l++)
      printf("\troot_percent%d        = %f\n",l+1,hru->veg_con.root[l]);
    for (j = 0; j < 12; j++) 
      printf("\tLAI[%02d]             = %f\n",j,state->veg_lib[vegclass].LAI[j]);
    for (j = 0; j < 12; j++) 
      printf("\talbedo[%02d]          = %f\n",j,state->veg_lib[vegclass].albedo[j]);
    for (j = 0; j < 12; j++) 
      printf("\tdisplacement[%02d]    = %f m\n",j,
	     state->veg_lib[vegclass].displacement[j]);
    for (j = 0; j < 12; j++) 
      printf("\troughness[%02d]       = %f m\n",j,state->veg_lib[vegclass].roughness[j]);
    printf("\twind_h                  = %f s/m\n", state->veg_lib[vegclass].wind_h);
  }
}

