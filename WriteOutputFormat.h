#ifndef WRITEOUTPUTFORMAT_H_
#define WRITEOUTPUTFORMAT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "vicNl.h"

/*
 * This is the abstract class which defines the interface to the
 * methods which each different output format must implement.
 */
class WriteOutputFormat {
public:
  WriteOutputFormat(const ProgramState* state) {
    compressAfterWrite = state->options.COMPRESS;
  }
  virtual ~WriteOutputFormat() {
    for (unsigned int filenum = 0; filenum< dataFiles.size(); filenum++) {
      delete dataFiles[filenum];
    }
  }

  virtual const char* getDescriptionOfOutputType() = 0;
  virtual void openFile() = 0;
  virtual void cleanup() {
    // Close Output Files and free memory (if any).
    for (unsigned int filenum=0; filenum<dataFiles.size(); filenum++) {
      if (dataFiles[filenum]->fh != NULL) {
        fclose(dataFiles[filenum]->fh);
        dataFiles[filenum]->fh = NULL;
      }
    }
    if(compressAfterWrite) compressFiles();
  }
  virtual void compressFiles() = 0;
  /***************************************************************
   Write output files using default VIC ASCII or BINARY formats
   - multiple files, all variables, no truncation

   see VIC web page for format details:
   www.hydro.washington.edu/Lettenmaier/Models/VIC/VIChome.html
   ***************************************************************/
  virtual void write_data_one_cell(out_data_struct *out_data, const dmy_struct *dmy, size_t numTimeRecords, const ProgramState* state) = 0;
  virtual void write_header(out_data_struct *out_data, const dmy_struct *dmy, const ProgramState* state) = 0;

  std::vector<out_data_file_struct*> dataFiles;
  std::string netCDFOutputFileName;
  float lat;
  float lon;
private:
  bool compressAfterWrite;
};

#endif /* WRITEOUTPUTFORMAT_H_ */
