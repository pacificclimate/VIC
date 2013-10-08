#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void free_dist_prcp(dist_prcp_struct *prcp, 
		    int               Nveg)
/**********************************************************************
	free_dist_prcp	Keith Cherkauer		March 1998

  This routine frees all memory allocated down the distributed 
  precipitation data structure.  This include all grid cell specific
  variables (soil, vegetation, energy, snow).

  modifications:
  06-24-98 modified to account for redesign of distributed precipitation
           data structures                                          KAC
  2007-Apr-21 Replaced loop over Nveg to loop over Nitems, so that lake-
	      specific veg tiles could be freed.			TJB
  2009-Jul-31 Removed extra veg tile for lake/wetland.			TJB

**********************************************************************/
{
  int Ndist = 2;
  int Nitems = Nveg + 1;

  free(prcp->mu);

}
