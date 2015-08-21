#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

static char vcid[] = "$Id$";

void initialize_veg(std::vector<HRU>& elements, int dist)
/**********************************************************************
  initialize_veg		Dag Lohmann	 January 1996

  This routine initializes the vegetation variable array.

  Modifications:
  07-13-98 modified to initialize vegetation structure for all 
           defined elevation bands                                 KAC
  11-18-02 modified to get the maximum number of vegetation types
           passed to it.  This allows the maximum number of vegetation
           types to include the wetland vegetation fraction when the 
           lake model is active.                                  LCB

**********************************************************************/
{
  for (std::vector<HRU>::iterator it = elements.begin(); it != elements.end(); ++it) {
    it->veg_var[dist].Wdew = 0.0;
    it->veg_var[dist].throughfall = 0.0;
  }
}
