#ifndef WRITEOUTPUTNETCDF_H_
#define WRITEOUTPUTNETCDF_H_

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "user_def.h"
#include "vicNl.h"
//#include "WriteOutputFormat.h"

namespace netCDF {
  class NcFile;
}

//class WriteOutputNetCDF: public WriteOutputFormat {
class WriteOutputNetCDF {
public:
	WriteOutputNetCDF(const ProgramState* state);
	~WriteOutputNetCDF() ;
  const char* getDescriptionOfOutputType();
  // This should only be called once per invocation of VIC. It creates a fresh netCDF output file.
  void initializeFile(const ProgramState*, const out_data_struct*);
  void openFile();
  void write_data_one_cell(out_data_struct *out_data, const dmy_struct *dmy, int dt, const ProgramState* state);
  void write_data_all_cells(std::vector<out_data_struct*>& all_out_data, out_data_file_struct *out_data_files_template, const dmy_struct *dmy, int dt, const ProgramState *state);
  void write_header(out_data_struct *out_data, const dmy_struct *dmy, const ProgramState* state);
  int getLengthOfTimeDimension(const ProgramState* state);
  int getTimeIndex(const dmy_struct* curTime, const int timeIndexDivisor, const ProgramState* state);
  netCDF::NcFile* netCDF;
  std::string netCDFOutputFileName;
  std::vector<out_data_file_struct*> dataFiles;
  int timeIndexDivisor;
  float lat;
  float lon;
};

#endif /* WRITEOUTPUTNETCDF_H_ */
