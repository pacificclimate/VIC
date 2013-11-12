#ifndef STATEIOCONTEXT_H_
#define STATEIOCONTEXT_H_

#include <string>

#include "StateIO.h"
#include "vicNl_def.h"

/*
 * This is a convenience class which just provides automatic memory management of a StateIO object.
 * It sets up a StateIO object to the right format based on the program options, and deletes it
 * when the instance goes out of scope.
 */
class StateIOContext {
public:
  StateIOContext(std::string filename, StateIO::IOType ioType, const ProgramState* state);
  virtual ~StateIOContext();
  StateIO* stream;
};

#endif /* STATEIOCONTEXT_H_ */
