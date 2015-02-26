#ifndef WRITEOUTPUTASCII_H_
#define WRITEOUTPUTASCII_H_

#include "WriteOutputFormat.h"

class WriteOutputAscii: public WriteOutputFormat {
public:
  WriteOutputAscii(const ProgramState* state) : WriteOutputFormat(state) {}
  const char* getDescriptionOfOutputType();
  void openFile();
  void compressFiles();
  void write_data(out_data_struct *out_data, const dmy_struct *dmy, int dt, const ProgramState* state);
//  void write_data(out_data_struct *out_data, const dmy_struct *dmy, int dt, ProgramState* state);
  void write_header(out_data_struct *out_data, const dmy_struct *dmy, const ProgramState* state);
};

#endif /* WRITEOUTPUTASCII_H_ */
