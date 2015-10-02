#ifndef WRITEOUTPUTNETCDF_H_
#define WRITEOUTPUTNETCDF_H_

//#include "WriteOutputFormat.h"
#include "user_def.h"
#include "vicNl.h"

#if NETCDF_OUTPUT_AVAILABLE

namespace netCDF {
  class NcFile;
}

//class WriteOutputAllCells: public WriteOutputFormat {
class WriteOutputAllCells {
public:
  WriteOutputAllCells(const ProgramState* state);
  ~WriteOutputNetCDF();
//  void write_all_cells_output(out_data_struct *all_out_data, const dmy_struct *dmy, int dt, const ProgramState* state);
  void write_all_cells_output(out_data_struct *all_out_data, out_data_file_struct* out_data_file, const dmy_struct *dmy, int dt, const ProgramState* state);

private:
  int getTimeIndex(const dmy_struct* curTime, const int timeIndexDivisor, const ProgramState* state);
  netCDF::NcFile* netCDF;
  int timeIndexDivisor;
};

#endif /* NETCDF_OUTPUT_AVAILABLE */

#endif /* WRITEOUTPUTNETCDF_H_ */
