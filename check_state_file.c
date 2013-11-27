#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>

#include "vicNl.h"
#include "StateIOContext.h"

static char vcid[] = "$Id$";

void check_state_file(char *init_state_name, ProgramState *state)
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

  StateIOContext context(init_state_name, StateIO::Reader, state);
  StateHeader header = context.stream->readHeader();

  if ( header.nLayer != state->options.Nlayer ) {
    std::stringstream ss;
    ss << "The number of soil moisture layers in the model state file (" << header.nLayer << ")";
    ss << "does not equal that defined in the global control file (" << state->options.Nlayer << ").  Check your input files.";
    throw VICException(ss.str());
  }
  if ( header.nNode != state->options.Nnode ) {
    std::stringstream ss;
    ss << "The number of soil thermal nodes in the model state file (" << header.nNode << ")";
    ss << "does not equal that defined in the global control file (" << state->options.Nnode << ").  Check your input files.";
    throw VICException(ss.str());
  }
}
