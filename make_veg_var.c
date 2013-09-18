#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

veg_var_struct **make_veg_var(int veg_type_num, const int NUM_SNOW_BAND)
/**********************************************************************
	make_veg_var	Dag Lohman		January 1996

  This routine makes an array of vegitation variables for each vegitation
  type.

  Modifications:
  07-13-98 modified to add structure definitions for all defined 
           elevation bands                                       KAC

**********************************************************************/
{
  int              i;
  veg_var_struct **temp;

  temp = (veg_var_struct **) calloc(veg_type_num, 
				    sizeof(veg_var_struct *));
  for(i=0;i<veg_type_num;i++)
    temp[i] = (veg_var_struct *) calloc(NUM_SNOW_BAND,
					sizeof(veg_var_struct));
  return temp;
}
