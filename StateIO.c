#include "StateIO.h"

StateIO::StateIO(std::string filename, const ProgramState* state) : filename(filename), state(state) {
}

StateIO::~StateIO() {
}

int StateIO::writeNewline() {
  return 0;
}
