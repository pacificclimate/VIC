#ifndef WRITEOUTPUTBINARY_H_
#define WRITEOUTPUTBINARY_H_

#include "WriteOutputFormat.h"

class WriteOutputBinary: public WriteOutputFormat {
public:
  const char* getDescriptionOfOutputType();
  FILE* openFile(const char* string);
  void write_data(out_data_file_struct *out_data_files, out_data_struct *out_data, const dmy_struct *dmy, int dt, const ProgramState* state);
  void write_header(out_data_file_struct *out_data_files, out_data_struct *out_data, const dmy_struct *dmy, const ProgramState* state);

private:
  void prepareDataForWriting(out_data_struct* out_data);
};

#endif /* WRITEOUTPUTBINARY_H_ */
