/*
 * OutputData.h
 *
 *  Created on: Aug 30, 2016
 *      Author: mfischer
 */

#ifndef OUTPUTDATA_H_
#define OUTPUTDATA_H_

#include <stdlib.h>
#include <string>

/*******************************************************
  This class stores output information for one variable.
  It replaces the old C out_data_struct and is more memory-safe.
  *******************************************************/
class OutputData {
public:
	OutputData();
	~OutputData();
  std::string	varname;
  int		write;       /* FALSE = don't write; TRUE = write */
// FIXME: format, type, and mult members should go away once we extricate ASCII and BINARY format output code
  std::string	format; /* format, when written to an ascii file;
		                should match the desired fprintf format specifier, e.g. %.4f */
  int		type;        /* type, when written to a binary file;
		                OUT_TYPE_USINT  = unsigned short int
		                OUT_TYPE_SINT   = short int
		                OUT_TYPE_FLOAT  = single precision floating point
		                OUT_TYPE_DOUBLE = double precision floating point */
  float		mult;        /* multiplier, when written to a binary file */
  int		aggtype;     /* type of aggregation to use;
				AGG_TYPE_AVG    = take average value over agg interval
				AGG_TYPE_BEG    = take value at beginning of agg interval
				AGG_TYPE_END    = take value at end of agg interval
				AGG_TYPE_MAX    = take maximum value over agg interval
				AGG_TYPE_MIN    = take minimum value over agg interval
				AGG_TYPE_SUM    = take sum over agg interval */
  int		nelem;       /* number of data values */
  double	*data;       /* array of data values */
  double	*aggdata;    /* array of aggregated data values */
};

#endif /* OUTPUTDATA_H_ */
