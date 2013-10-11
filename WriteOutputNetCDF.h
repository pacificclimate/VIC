#ifndef WRITEOUTPUTNETCDF_H_
#define WRITEOUTPUTNETCDF_H_

#include "WriteOutputFormat.h"

class WriteOutputNetCDF: public WriteOutputFormat {
public:
  WriteOutputNetCDF(const ProgramState* state) : WriteOutputFormat(state) {}
  const char* getDescriptionOfOutputType();
  void openFile();
  void compressFiles();
  void write_data(out_data_struct *out_data, const dmy_struct *dmy, int dt, const ProgramState* state);
  void write_header(out_data_struct *out_data, const dmy_struct *dmy, const ProgramState* state);
};

#endif /* WRITEOUTPUTNETCDF_H_ */
