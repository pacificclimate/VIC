#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void vicerror(const char *error_text)
/**********************************************************************
	vicerror.c	Keith Cherkauer		April 23, 1997

  This subroutine was written to handle numerical errors within the
  VIC model.  This will flush all file buffers so that all records 
  that have been run will be written to disk before the model is exited.

  Modifications:
  2006-Sep-23 Implemented flexible output configuration; uses the new
              out_data and out_data_files structures. TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct. TJB

**********************************************************************/
{
        filenames_struct fnames;

	/* turn off compression of last set of files */
  bool compress = false;

	fprintf(stderr,"VIC model run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now writing output files...\n");
	//TODO: this is broken since global state is not updated with each cell, and this would not work with multiple threads anyways
  //close_files(&(state->Error.filep), state->Error.out_data_files, &fnames, compress);
	fprintf(stderr,"...now exiting to system...\n");
        fflush(stdout);
        fflush(stderr);
	_Exit(1);
}
