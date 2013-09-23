#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

int  dist_prec(cell_info_struct* cell,
               const dmy_struct    *dmy,
               filep_struct        *filep,
               out_data_file_struct *out_data_files,
               out_data_struct     *out_data,
               int                  time_step_record,
               const int            cellnum,
               char                 NEWCELL,
               const ProgramState  *state) {
/**********************************************************************
  dist_prec		Keith Cherkauer		October 9, 1997

  This subroutine calls the solution routines for a single grid cell
  for one time step.  It also controls the distribution of precipitation
  and reassembles grid cell data for output.

  The fractional coverage of precipitation over an area or grid cell, 
  mu, is estimated using the equation from Fan et. al. (1996).  The 
  coefficient, 0.6, was selected for the Arkansas - Red River Basin and
  was found using precipitation records on a 100km x 100km area.  It
  may not be applicable to all regions, please check the reference

  References:

  Modifications:
  11-30-98 Added counter to assure that a storm has been stopped
           for at least one day, before allowing the model to 
	   average soil moisture when a new precipitation event
	   arrives.                                             KAC
  03-05-01 Fixed error in which distributed precipitation accounting
           variables (DRY_TIME, STILL_STORM, ANY_SNOW) were used 
           within the vegetation loop, but did not store separate
           values for each vegetation type.                     KAC
  03-12-03 Modifed to add AboveTreeLine to soil_con_struct so that
           the model can make use of the computed treeline.     KAC
  03-27-03 Modified calculation of DRY_TIME.  Originally the check
           to see if a new storm was warranted checked if DRY_TIME
           was greater than 24/dt.  However, DRY_TIME is incremented
           by dt, so it was checking hours against time steps.  The
           division by dt has been removed, so a new storm starts if
           the cell has been drying for a full 24 hours.     RS & KAC
  04-10-03 Modified to store STILL_STORM and DRY_TIME in the model
           statefile, so that full conditions will be preserved.  KAC
  01-Nov-04 Added support for state files containing SPATIAL_FROST and
	    LAKE_MODEL state variables.					TJB
  02-Feb-05 Modified to save state file at the end of the final timestep
	    of the date indicated by STATEYEAR, STATEMONTH, and STATEDAY
	    in the global parameter file.				TJB
  2005-Mar-24 Modified parameter list of put_data() to accomodate support
	      for ALMA variables.					TJB
  2006-Sep-23 Implemented flexible output configuration; uses new out_data,
	      out_data_files, and save_data structures.			TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct.	TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2007-Apr-04 Modified to handle grid cell errors by returning to the
              main subroutine, rather than ending the simulation.	GCT/KAC
  2008-Oct-23 Modified call to put_data() to store ErrorFlag.		TJB
  2009-Mar-03 Modified routine to store put_data() error in ErrorFlag2 and 
	      return a single ERROR value if an error occurs.		KAC via TJB
  2009-Jun-19 Added T flag to indicate whether TFALLBACK occurred.	TJB
  2009-Sep-28 Added logic for initial (pre-simulation) call to put_data.TJB

**********************************************************************/

  char ANY_SNOW[MAX_VEG];
  int veg, i;
  int month;
  int ErrorFlag, ErrorFlag2;
  double Wdmax;
  double NEW_MU;

  if (state->options.DIST_PRCP) {

    /*******************************************
     Controls Distributed Precipitation Model
     *******************************************/

    NEW_MU = 1.0 - exp(-state->options.PREC_EXPT * cell->atmos[time_step_record].prec[state->NR]);
    for (veg = 0; veg <= cell->veg_con[0].vegetat_type_num; veg++) {
      ANY_SNOW[veg] = FALSE;
      for (i = 0; i < state->options.SNOW_BAND; i++)
        /* Check for snow on ground or falling */
        if (cell->prcp.snow[veg][i].swq > 0 || cell->prcp.snow[veg][i].snow_canopy > 0.)
          ANY_SNOW[veg] = TRUE;
      if (ANY_SNOW[veg] || cell->atmos[time_step_record].snowflag[state->NR]) {
        /* If snow present, mu must be set to 1. */
        NEW_MU = 1.;
        if (time_step_record == 0) {
          /* Set model variables if first time step */
          cell->prcp.mu[veg] = NEW_MU;
          if (cell->atmos[time_step_record].prec[state->NR] > 0)
            cell->init_STILL_STORM[veg] = TRUE;
          else
            cell->init_STILL_STORM[veg] = FALSE;
          cell->init_DRY_TIME[veg] = 0;
        }
        ANY_SNOW[veg] = TRUE;
      } else {
        if (time_step_record == 0) {
          if (cell->atmos[time_step_record].prec[state->NR] == 0) {
            /* If first time step has no rain, than set mu to 1. */
            cell->prcp.mu[veg] = 1.;
            NEW_MU = 1.;
            cell->init_STILL_STORM[veg] = TRUE;
            cell->init_DRY_TIME[veg] = 24;
          } else {
            /* If first time step has rain, then set mu based on intensity */
            cell->prcp.mu[veg] = NEW_MU;
            cell->init_STILL_STORM[veg] = TRUE;
            cell->init_DRY_TIME[veg] = 0;
          }
        } else if (cell->atmos[time_step_record].prec[state->NR] == 0 && cell->init_DRY_TIME[veg] >= 24.) {
          /* Check if storm has ended */
          NEW_MU = cell->prcp.mu[veg];
          cell->init_STILL_STORM[veg] = FALSE;
          cell->init_DRY_TIME[veg] = 0;
        } else if (cell->atmos[time_step_record].prec[state->NR] == 0) {
          /* May be pause in storm, keep track of pause length */
          NEW_MU = cell->prcp.mu[veg];
          cell->init_DRY_TIME[veg] += state->global_param.dt;
        }
      }

      if (!cell->init_STILL_STORM[veg]
          && (cell->atmos[time_step_record].prec[state->NR] > STORM_THRES || ANY_SNOW[veg])) {
        /** Average soil moisture before a new storm **/
        ErrorFlag = initialize_new_storm(cell->prcp.cell, cell->prcp.veg_var, veg,
            cell->veg_con[0].vegetat_type_num, time_step_record, cell->prcp.mu[veg], NEW_MU, state);
        if (ErrorFlag == ERROR)
          return (ERROR);

        cell->init_STILL_STORM[veg] = TRUE;
        cell->prcp.mu[veg] = NEW_MU;
      } else if (NEW_MU != cell->prcp.mu[veg] && cell->init_STILL_STORM[veg]) {
        /** Redistribute soil moisture during the storm if mu changes **/
        if (dmy[time_step_record].day == 1 && dmy[time_step_record].hour == 0) {
          month = dmy[time_step_record].month - 2;
          if (month < 0)
            month = 11;
        } else
          month = dmy[time_step_record].month - 1;
        if (veg < cell->veg_con[0].vegetat_type_num)
          Wdmax = state->veg_lib[cell->veg_con[veg].veg_class].Wdmax[month];
        else
          Wdmax = 0;
        redistribute_during_storm(cell->prcp.cell, cell->prcp.veg_var, veg,
            cell->veg_con[0].vegetat_type_num, time_step_record, Wdmax, cell->prcp.mu[veg], NEW_MU,
            cell->soil_con.max_moist, state);
        cell->prcp.mu[veg] = NEW_MU;
      }
    }

    /** Solve model time step **/
    ErrorFlag = full_energy(NEWCELL, cellnum, time_step_record, &cell->atmos[time_step_record], &cell->prcp, dmy,
        &cell->lake_con, &cell->soil_con, cell->veg_con, &cell->writeDebug, state);

  }

  else { /* FIXME dumb flow control */

    /**************************************************
     Controls Grid Cell Averaged Precipitation Model
     **************************************************/

    ErrorFlag = full_energy(NEWCELL, cellnum, time_step_record, &cell->atmos[time_step_record], &cell->prcp, dmy,
        &cell->lake_con, &cell->soil_con, cell->veg_con, &cell->writeDebug, state);

  }

  /**************************************************
   Write cell average values for current time step
   **************************************************/

  ErrorFlag2 = put_data(cell, out_data_files, out_data, &dmy[time_step_record],
      time_step_record, state);
  if (ErrorFlag2 == ERROR)
    ErrorFlag = ERROR;

  /************************************
   Save model state at assigned date
   (after the final time step of the assigned date)
   ************************************/

  if (filep->statefile != NULL
      && (dmy[time_step_record].year == state->global_param.stateyear
          && dmy[time_step_record].month == state->global_param.statemonth
          && dmy[time_step_record].day == state->global_param.stateday
          && (time_step_record + 1 == state->global_param.nrecs
              || dmy[time_step_record + 1].day != state->global_param.stateday)))
    write_model_state(&cell->prcp, cell->veg_con[0].vegetat_type_num,
        cell->soil_con.gridcel, filep, &cell->soil_con, cell->init_STILL_STORM, cell->init_DRY_TIME, cell->lake_con, state);

  return (ErrorFlag);

}
