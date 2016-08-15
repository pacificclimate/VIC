#ifndef WRITEOUTPUTNETCDF_H_
#define WRITEOUTPUTNETCDF_H_

#include <string>
//#include <map>
#include "user_def.h"
#include "WriteOutputFormat.h"

#if NETCDF_OUTPUT_AVAILABLE

namespace netCDF {
  class NcFile;
}

class WriteOutputNetCDF: public WriteOutputFormat {
public:
  WriteOutputNetCDF(const ProgramState* state);
  ~WriteOutputNetCDF();
  const char* getDescriptionOfOutputType();
  // This should only be called once per invocation of VIC. It creates a fresh netCDF output file.
  void initializeFile(const ProgramState*, const out_data_struct*);
  void openFile();
  void compressFiles();
  void write_data_one_cell(out_data_struct *out_data, const dmy_struct *dmy, int numTimeRecords, const ProgramState* state);
  void write_data_all_cells(std::vector<out_data_struct*>& all_out_data, out_data_file_struct *out_data_files_template, const dmy_struct *dmy, int dt, const ProgramState *state);
  void write_header(out_data_struct *out_data, const dmy_struct *dmy, const ProgramState* state);
  int getLengthOfTimeDimension(const ProgramState* state);
  int getTimeIndex(const dmy_struct* curTime, const int timeIndexDivisor, const ProgramState* state);
  netCDF::NcFile* netCDF;
  int timeIndexDivisor;
};

#endif /* NETCDF_OUTPUT_AVAILABLE */

#endif /* WRITEOUTPUTNETCDF_H_ */
