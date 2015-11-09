#ifndef WRITEOUTPUTALLCELLS_H_
#define WRITEOUTPUTALLCELLS_H_

#include "WriteOutputNetCDF.h"
#include "user_def.h"
#include "vicNl.h"

class WriteOutputAllCells: public WriteOutputNetCDF {
public:
  WriteOutputAllCells(const ProgramState* state);
  ~WriteOutputAllCells();

  void write_data(std::vector<out_data_struct*>& all_out_data, out_data_file_struct *out_data_files_template, const dmy_struct *dmy, int dt, const ProgramState *state);

  void closeFile();
};

#endif /* WRITEOUTPUTALLCELLS_H_ */
