#include "StateIO.h"

StateIO::StateIO(std::string filename, IOType type, const ProgramState* state) : filename(filename), state(state), ioType(type) {
}

StateIO::~StateIO() {
}

int StateIO::writeNewline() {
  return 0;
}
