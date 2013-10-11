#ifndef WRITEOUTPUTBINARY_H_
#define WRITEOUTPUTBINARY_H_

#include "WriteOutputFormat.h"

class WriteOutputBinary: public WriteOutputFormat {
public:
  WriteOutputBinary(const ProgramState* state) : WriteOutputFormat(state) {}
  const char* getDescriptionOfOutputType();
  void openFile();
  void compressFiles();
  void write_data(out_data_struct *out_data, const dmy_struct *dmy, int dt, const ProgramState* state);
  void write_header(out_data_struct *out_data, const dmy_struct *dmy, const ProgramState* state);

private:
  void prepareDataForWriting(out_data_struct* out_data);
};

#endif /* WRITEOUTPUTBINARY_H_ */
