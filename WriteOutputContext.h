#ifndef WRITEOUTPUTCONTEXT_H_
#define WRITEOUTPUTCONTEXT_H_

#include "WriteOutputFormat.h"
class WriteOutputFormat;

/* This class provides the context for the different output format types.
 * This is the strategy design pattern.
 * This class takes care of creating and cleaning up the correct output format class.
 * After setting up this object, you can call the various methods of the abstract outputFormat class
 * and through polymorphism, the correct format will be written.
 */
class WriteOutputContext {
public:
  WriteOutputContext(const ProgramState* state);
  virtual ~WriteOutputContext();
  WriteOutputFormat* outputFormat;
};

#endif /* WRITEOUTPUTCONTEXT_H_ */
