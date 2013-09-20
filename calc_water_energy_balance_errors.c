#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double calc_water_balance_error(int    rec,
				double inflow,
				double outflow,
				double storage,
				int    Nrecs,
				CellBalanceErrors* errors) {
  /***************************************************************
  calc_water_balance_error  Keith Cherkauer        April 1998

  This subroutine computes the overall model water balance, and 
  warns the model user if large errors are found.

  Modifications:
  2007-Aug-22 Added error as return value.  JCA
***************************************************************/

  double error;

  if(rec<0) {
    errors->water_last_storage = storage;
    errors->water_cum_error    = 0.;
    errors->water_max_error    = 0.;
    return(0.0);
  }
  else {
    error = inflow - outflow - (storage - errors->water_last_storage);
    errors->water_cum_error += error;
    if(fabs(error)>fabs(errors->water_max_error) && fabs(error)>1e-5) {
      errors->water_max_error = error;
      fprintf(stderr,"Maximum Moist Error:\t%i\t%.5f\t%.5f\n",
	      rec,error,errors->water_cum_error);
    }
    if(rec==Nrecs-1) {
      fprintf(stderr,"Total Cumulative Water Error for Grid Cell = %.4f\n",
	      errors->water_cum_error);
    }
    errors->water_last_storage = storage;

    return(error);
  }

}

void calc_energy_balance_error(int    rec,
			       double net_rad,
			       double latent,
			       double sensible,
			       double grnd_flux,
			       double snow_fluxes,
			       int    Nrecs,
			       CellBalanceErrors* errors) {
/***************************************************************
  calc_energy_balance_error   Keith Cherkauer     April 1998

  This subroutine computes the overall model energy balance, and
  reports the maximum time step error above a thresehold to the
  user.  The total cumulative error for the grid cell is also 
  computed and reported at the end of the model run.

***************************************************************/

  double error;

  if(rec<0) {
    errors->energy_cum_error = 0;
    errors->energy_max_error = 0;
  }
  else {
    error = net_rad - latent - sensible - grnd_flux + snow_fluxes;
    errors->energy_cum_error += error;
    if(fabs(error)>fabs(errors->energy_max_error) && fabs(error)>0.001) {
      errors->energy_max_error = error;
      if ( rec > 0 ) 
	fprintf(stderr,"Maximum Energy Error:\t%i\t%.4f\t%.4f\n",
		rec,error,errors->energy_cum_error/(double)rec);
      else 
	fprintf(stderr,"Maximum Energy Error:\t%i\t%.4f\t%.4f\n",
		rec,error,errors->energy_cum_error);
    }
    if(rec==Nrecs-1) {
      fprintf(stderr,"Total Cumulative Energy Error for Grid Cell = %.4f\n",
	      errors->energy_cum_error/(double)rec);
    }
  }
}

