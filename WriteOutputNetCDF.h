#ifndef WRITEOUTPUTNETCDF_H_
#define WRITEOUTPUTNETCDF_H_

#include <string>

#include "WriteOutputFormat.h"
namespace netCDF {
  class NcFile;
}

class WriteOutputNetCDF: public WriteOutputFormat {
public:
  WriteOutputNetCDF(const ProgramState* state);
  ~WriteOutputNetCDF();
  const char* getDescriptionOfOutputType();
  // This should only be called once per invocation of VIC. It creates a fresh netCDF output file.
  void initializeFile(const ProgramState*);
  void openFile();
  void compressFiles();
  void write_data(out_data_struct *out_data, const dmy_struct *dmy, int dt, const ProgramState* state);
  void write_header(out_data_struct *out_data, const dmy_struct *dmy, const ProgramState* state);
private:
  std::string netCDFOutputFileName;
  netCDF::NcFile* netCDF;
};

#endif /* WRITEOUTPUTNETCDF_H_ */
