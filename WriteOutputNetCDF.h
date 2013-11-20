#ifndef WRITEOUTPUTNETCDF_H_
#define WRITEOUTPUTNETCDF_H_

#include <string>
#include <map>
#include "user_def.h"
#include "WriteOutputFormat.h"

#if NETCDF_OUTPUT_AVAILABLE

namespace netCDF {
  class NcFile;
}

using std::string;
class VariableMetaData {
public:
  VariableMetaData() {}
  VariableMetaData(string units, string name, string standardName, string longName, string cellMethods, double scalingFactor = 1.0, double addFactor = 0.0, bool isBands = false): units(units), name(name), standardName(standardName), longName(longName), cellMethods(cellMethods), scalingFactor(scalingFactor), addFactor(addFactor), isBands(isBands) {}
  string units, name, standardName, longName, cellMethods;
  double scalingFactor;
  double addFactor;
  bool isBands;
};

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
  std::map<std::string, VariableMetaData> getMapping(bool isHourly = false);
  std::map<std::string, VariableMetaData> mapping;
  std::string netCDFOutputFileName;
  netCDF::NcFile* netCDF;
  int timeIndexDivisor;
};

#endif /* NETCDF_OUTPUT_AVAILABLE */

#endif /* WRITEOUTPUTNETCDF_H_ */
