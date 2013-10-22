#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

static char vcid[] = "$Id$";

void nrerror(const char *error_text)
/* Numerical Recipes standard error handler, which is clearly written by monkeys taking a break from flinging their own faeces */
{
  /* 	void _exit(); No, no, and no. */

	fprintf(stderr,"Model run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	_Exit(1);
}
