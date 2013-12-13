#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
 
static char vcid[] = "$Id$";

void free_vegcon(cell_info_struct& cell)
/**********************************************************************
  free_vegcon.c	            Keith Cherkauer	  September 25, 1998

  This subroutine frees all components of the veg_con structure.

**********************************************************************/
{
  for(std::vector<HRU>::iterator hru = cell.prcp.hruList.begin(); hru != cell.prcp.hruList.end(); ++hru) {
    free(hru->veg_con.zone_depth);
    free(hru->veg_con.zone_fract);
  }
}
