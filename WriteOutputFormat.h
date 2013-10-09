#ifndef WRITEOUTPUTFORMAT_H_
#define WRITEOUTPUTFORMAT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

/*
 * This is the abstract class which defines the interface to the
 * methods which each different output format must implement.
 */
class WriteOutputFormat {
public:
  WriteOutputFormat() {}
  virtual ~WriteOutputFormat() {}
  virtual const char* getDescriptionOfOutputType() = 0;
  virtual FILE* openFile(const char* string) = 0;
  /***************************************************************
   Write output files using default VIC ASCII or BINARY formats
   - multiple files, all variables, no truncation

   see VIC web page for format details:
   www.hydro.washington.edu/Lettenmaier/Models/VIC/VIChome.html
   ***************************************************************/
  virtual void write_data(out_data_file_struct *out_data_files, out_data_struct *out_data, const dmy_struct *dmy, int dt, const ProgramState* state) = 0;
  virtual void write_header(out_data_file_struct *out_data_files, out_data_struct *out_data, const dmy_struct *dmy, const ProgramState* state) = 0;
};

#endif /* WRITEOUTPUTFORMAT_H_ */
