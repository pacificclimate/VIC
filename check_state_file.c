#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vicNl.h"

static char vcid[] = "$Id$";

FILE *check_state_file(char *init_state_name, ProgramState *state)
/*********************************************************************
  check_state_file      Keith Cherkauer           April 17, 2000

  This subroutine opens a model state file and verifys that the 
  starting date, number of layers and number of thermal nodes in the 
  file agrees with what was defined in the model global control file.

  Modifications:
  04-10-03 modified to open and read from a binary state file.    KAC
  04-10-03 modified to compute record where state file starts,
           this allows VIC to read in the full array of atmospheric
           forcing data but start the simulation at the same time
           step as the state file.  This should eliminate the 
           problems associated with restarting the model with an 
           incomplete record of forcing data, which can lead to 
           differences in the interpolated sub-daily forcings.    KAC
  06-03-03 modified to handle both ASCII and BINARY state files.  KAC
  2006-08-23 Changed order of fread/fwrite statements from ...1, sizeof...
             to ...sizeof, 1,... GCT
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct. TJB

*********************************************************************/
{
  FILE   *init_state;
  char    filename[MAXSTRING];
  char    ErrStr[MAXSTRING];
  double  Nsum;
  int     tmp_Nlayer;
  int     tmp_Nnodes;
  int     startday, startmonth, startyear;

  /* open state file */
  if ( state->options.STATE_FORMAT == StateOutputFormat::BINARY_STATEFILE )
    init_state = open_file(init_state_name,"rb");
  else 
    init_state = open_file(init_state_name,"r");

  /* Check state date information */
  if ( state->options.STATE_FORMAT == StateOutputFormat::BINARY_STATEFILE ) {
    fread( &startyear, sizeof(int), 1, init_state );
    fread( &startmonth, sizeof(int), 1, init_state );
    fread( &startday, sizeof(int), 1, init_state );
  }
  else {
    fscanf(init_state,"%d %d %d\n", &startyear, &startmonth, &startday);
  }

  /* Check simulation options */
  if ( state->options.STATE_FORMAT == StateOutputFormat::BINARY_STATEFILE ) {
    fread( &tmp_Nlayer, sizeof(int), 1, init_state );
    fread( &tmp_Nnodes, sizeof(int), 1, init_state );
  }
  else {
    fscanf(init_state,"%d %d\n", &tmp_Nlayer, &tmp_Nnodes);
  }
  if ( tmp_Nlayer != state->options.Nlayer ) {
    sprintf(ErrStr,"The number of soil moisture layers in the model state file (%d) does not equal that defined in the global control file (%d).  Check your input files.", tmp_Nlayer, state->options.Nlayer);
    nrerror(ErrStr);
  }
  if ( tmp_Nnodes != state->options.Nnode ) {
    sprintf(ErrStr,"The number of soil thermal nodes in the model state file (%d) does not equal that defined in the global control file (%d).  Check your input files.", tmp_Nnodes, state->options.Nnode);
    nrerror(ErrStr);
  }

  return(init_state);

}
