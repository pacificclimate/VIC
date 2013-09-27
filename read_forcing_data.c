#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>
 
static char vcid[] = "$Id$";

double **read_forcing_data(FILE                **infile,
                           int                  *ncids,
			   global_param_struct   global_param,
                           soil_con_struct      *soil_con,
                           const ProgramState   *state)
/**********************************************************************
  read_forcing_data    Keith Cherkauer      January 10, 2000

  This subroutine controls the order and number of forcing variables
  read from the forcing data files.  Two forcing files are allowed, 
  variables, time step and file format must be defined in the global
  control file.

**********************************************************************/
{
  char                 errorstr[MAXSTRING];
  int                  i;
  double             **forcing_data;

  /** Allocate data arrays for input forcing data **/
  forcing_data = (double **)calloc(N_FORCING_TYPES,sizeof(double*));
  for(i=0;i<N_FORCING_TYPES;i++) 
    if (state->param_set.TYPE[i].SUPPLIED)
      forcing_data[i] = (double *)calloc((global_param.nrecs * state->NF),
			   sizeof(double));

  /** Read First Forcing Data File **/
  if(IS_VALID(state->param_set.FORCE_DT[0]) && state->param_set.FORCE_DT[0] > 0) {
    read_atmos_data(infile[0], ncids[0], 0, global_param.forceskip[0],
		    forcing_data, soil_con, state);
  }
  else {
    sprintf(errorstr,"ERROR: File time step must be defined for at least the first forcing file (FILE_DT).\n");
    vicerror(errorstr);
  }

  /** Read Second Forcing Data File **/
  if(IS_VALID(state->param_set.FORCE_DT[1]) && state->param_set.FORCE_DT[1] > 0) {
    read_atmos_data(infile[1], ncids[1], 1, global_param.forceskip[1],
		    forcing_data, soil_con, state);
  }

  return(forcing_data);

}
