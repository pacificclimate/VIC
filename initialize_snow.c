#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

static char vcid[] = "$Id$";

void initialize_snow (std::vector<HRU>& elements)
/**********************************************************************
	initialize_snow		Keith Cherkauer		January 22, 1997

  This routine initializes the snow variable arrays for each new
  grid cell.

  VARIABLES INITIALIZED:
    snow[i][j].snow;	          TRUE = snow, FALSE = no snow 
    snow[i][j].last_snow;         time steps since last snowfall 
    snow[i][j].snow_canopy;       amount of snow on canopy (m) 
    snow[i][j].swq;               snow water equivalent of the entire pack (m) 
    snow[i][j].surf_water;        liquid water content of the surface 
                                  layer (m) 
    snow[i][j].pack_water;        liquid water content of the snow pack (m) 
    snow[i][j].surf_temp;         depth averaged temperature of the snow pack
                                  surface layer (C) 
    snow[i][j].pack_temp;         depth averaged temperature of the snow pack
                                  (C) 
    snow[i][j].vapor_flux;        depth of water evaporation, sublimation, or 
                                  condensation from snow pack (m) 
    snow[i][j].canopy_vapor_flux; depth of water evaporation, sublimation, or 
                                  condensation from intercepted snow (m) 
    snow[i][j].albedo;            snow surface albedo (fraction) 
    snow[i][j].coldcontent;       cold content of snow pack 
    snow[i][j].mass_error;        snow mass balance error 
    snow[i][j].density;	          snow density (kg/m^3) 
    snow[i][j].depth;	          snow depth (m) 
    snow[i][j].tmp_int_storage;   temporary canopy storage, used in 
                                  snow_canopy 
    snow[i][j].Qnet;              Net energy error in snow model 
    snow[i][j].band_elev;         median elevation of the current snow band 
    snow[i][j].prec_frac;         fraction of precipitation that falls in the
	  		          current snow band 

  modifications:
  07-09-98 modified to initialize snow variables for each defined
           snow elevation band.                                   KAC
  01-11-99 modified to read new initial snow conditions file format KAC
  04-17-00 removed call for read_initial_snow properties file, the
           file read is now incorporated into a single model state
           file.                                                  KAC
  xx-xx-01 modified to initialize spatially distributed snow variables. KAC
  11-18-02 modified to initialize blowing_snow variable.			LCB
  2006-Oct-16 Removed unused init_snow file.				TJB
  2009-Sep-19 Initializing last_snow to MISSING.			TJB
  2009-Sep-28 Added initialization of some terms that previously had
	      not been initialized.					TJB
  2010-Apr-24 Added initialization of surf_temp_fbcount and fbflag.	TJB

**********************************************************************/
{
  for (std::vector<HRU>::iterator it = elements.begin(); it != elements.end(); ++it) {
      // State variables
      it->snow.albedo            = 0.0;
      it->snow.canopy_albedo     = 0.0;
      it->snow.coldcontent       = 0.0;
      it->snow.coverage          = 0.0;
      it->snow.density           = 0.0;
      it->snow.depth             = 0.0;
      it->snow.last_snow         = INVALID_INT;
      it->snow.max_swq           = 0.0;
      it->snow.MELTING           = FALSE;
      it->snow.pack_temp         = 0.0;
      it->snow.pack_water        = 0.0;
      it->snow.snow              = FALSE;
      it->snow.snow_canopy       = 0.0;
      it->snow.store_coverage    = 0.0;
      it->snow.store_snow        = FALSE;
      it->snow.store_swq         = 0.0;
      it->snow.surf_temp         = 0.0;
      it->snow.surf_temp_fbflag  = 0;
      it->snow.surf_temp_fbcount = 0;
      it->snow.surf_water        = 0.0;
      it->snow.swq               = 0.0;
      it->snow.swq_slope         = 0.0;
      it->snow.tmp_int_storage   = 0.0;
      // Fluxes
      it->snow.blowing_flux      = 0.0;
      it->snow.canopy_vapor_flux = 0.0;
      it->snow.mass_error        = 0.0;
      it->snow.melt              = 0.0;
      it->snow.Qnet              = 0.0;
      it->snow.surface_flux      = 0.0;
      it->snow.transport         = 0.0;
      it->snow.vapor_flux        = 0.0;
  }
}
