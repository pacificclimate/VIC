#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vicNl.h"

static char vcid[] = "$Id$";

void calc_root_fractions(std::vector<HRU>& hruList,
			 soil_con_struct  *soil_con,
			 const ProgramState* state)
/**********************************************************************
  calc_root_fraction.c    Keith Cherkauer      September 24, 1998

  This routine computes the fraction of roots in each soil layer based
  on the root zone distribution defined in the vegetation parameter
  file.  Roots are assumed to be linearly distributed within each
  root zone.

**********************************************************************/
{
  char   ErrStr[MAXSTRING];
  int    Nhrus;
  int    veg;
  int    layer;
  int    zone;
  int    i;
  float  sum_depth;
  float  sum_fract;
  float  dum;
  double Zstep;
  double Zsum;
  double Lstep;
  double Lsum;
  double Zmin_fract;
  double Zmin_depth;
  double Zmax;

  Nhrus      = hruList.size();

  for(std::vector<HRU>::iterator hru = hruList.begin(); hru != hruList.end(); ++hru) {
    sum_depth  = 0;
    sum_fract  = 0;
    layer      = 0;
    Lstep      = soil_con->depth[layer];
    Lsum       = Lstep;
    Zsum       = 0;
    zone       = 0;
    
    while(zone<state->options.ROOT_ZONES) {
      Zstep = (double)hru->veg_con.zone_depth[zone];
      if((Zsum + Zstep) <= Lsum && Zsum >= Lsum - Lstep) {
	/** CASE 1: Root Zone Completely in Soil Layer **/
	sum_fract += hru->veg_con.zone_fract[zone];
      }
      else {
	/** CASE 2: Root Zone Partially in Soil Layer **/
	if(Zsum < Lsum - Lstep) {
	  /** Root zone starts in previous soil layer **/
	  Zmin_depth = Lsum - Lstep;
	  Zmin_fract = linear_interp(Zmin_depth,Zsum,Zsum+Zstep,0,
	      hru->veg_con.zone_fract[zone]);
	}
	else {
	  /** Root zone starts in current soil layer **/
	  Zmin_depth = Zsum;
	  Zmin_fract = 0.;
	}
	if(Zsum + Zstep <= Lsum) {
	  /** Root zone ends in current layer **/
	  Zmax = Zsum + Zstep;
	}
	else {
	  /** Root zone extends beyond bottom of current layer **/
	  Zmax = Lsum;
	}
        sum_fract += linear_interp(Zmax, Zsum, Zsum + Zstep, 0,
            hru->veg_con.zone_fract[zone]) - Zmin_fract;
      }

      /** Update Root Zone and Soil Layer **/
      if(Zsum + Zstep < Lsum) {
	Zsum += Zstep;
	zone ++;
      }
      else if(Zsum + Zstep == Lsum) {
	Zsum += Zstep;
	zone ++;
	if(layer<state->options.Nlayer) {
	  hru->veg_con.root[layer] = sum_fract;
	  sum_fract = 0.;
	}
	layer++;
	if(layer<state->options.Nlayer) {
	  Lstep  = soil_con->depth[layer];
	  Lsum  += Lstep;
	}
	else if(layer==state->options.Nlayer) {
	  Lstep  = Zsum + Zstep - Lsum;
	  if(zone<state->options.ROOT_ZONES-1) {
	    for(i=zone+1;i<state->options.ROOT_ZONES;i++) {
	      Lstep += hru->veg_con.zone_depth[i];
	    }
	  }
	  Lsum  += Lstep;
	}
      }
      else if(Zsum + Zstep > Lsum) {
	if(layer<state->options.Nlayer) {
	  hru->veg_con.root[layer] = sum_fract;
	  sum_fract = 0.;
	}
	layer++;
	if(layer<state->options.Nlayer) {
	  Lstep  = soil_con->depth[layer];
	  Lsum  += Lstep;
	}
	else if(layer==state->options.Nlayer) {
	  Lstep  = Zsum + Zstep - Lsum;
	  if(zone<state->options.ROOT_ZONES-1) {
	    for(i=zone+1;i<state->options.ROOT_ZONES;i++) {
	      Lstep += hru->veg_con.zone_depth[i];
	    }
	  }
	  Lsum  += Lstep;
	}
      }
	
    }

    if(sum_fract > 0 && layer >= state->options.Nlayer) {
      hru->veg_con.root[state->options.Nlayer-1] += sum_fract;
    }
    else if(sum_fract > 0) {
      hru->veg_con.root[layer] += sum_fract;
    }

    dum=0.;
    for (layer=0;layer<state->options.Nlayer;layer++) {
      if(hru->veg_con.root[layer] < 1.e-4) hru->veg_con.root[layer] = 0.;
      dum += hru->veg_con.root[layer];
    }
    if(dum == 0.0){
      sprintf(ErrStr,"Root fractions sum equals zero: %f , Vege Class: %d\n",
	      dum, hru->veg_con.veg_class);
      nrerror(ErrStr);
    }
    for (layer=0;layer<state->options.Nlayer;layer++) {
      hru->veg_con.root[layer] /= dum;
    }

  }

}

