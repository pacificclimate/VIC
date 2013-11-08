#include "vicNl.h"

static char vcid[] = "$Id$";

#include <string.h>


FILE *open_state_file(filenames_struct filenames, const ProgramState* state)
/*********************************************************************
  open_state_file      Keith Cherkauer           April 15, 2000

  This subroutine opens the model state file for output.

  Modifications:
  04-10-03 Modified to open and write to a binary state file.    KAC
  06-03-03 modified to handle both ASCII and BINARY state files.  KAC
  2005-11-29 SAVE_STATE is set in global param file, not in user_def.h GCT
  2005-12-06 Moved setting of statename to get_global_param     GCT
  2006-08-23 Changed order of fread/fwrite statements from ...1, sizeof...
             to ...sizeof, 1,... GCT
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      This included moving global->statename to filenames->statefile. TJB

*********************************************************************/
{
  FILE* statefile = NULL;

  StateIO* outputState = getStateIO(NULL, state);
  outputState->initializeOutput(&statefile, filenames.statefile, state);

  delete outputState;
  return statefile;
}

