#include "StateIO.h"

StateIO::StateIO(const ProgramState* state) : state(state) {
}

StateIO::~StateIO() {
}

int StateIO::writeNewline() {
  return 0;
}
