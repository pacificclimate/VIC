#ifndef STATEIO_H_
#define STATEIO_H_

#include "vicNl_def.h"
#include <string>

using std::string;
class StateVariableMetaData {
public:
  StateVariableMetaData() {}
  StateVariableMetaData(string units, string name, string standardName, string longName, string cellMethods, double scalingFactor = 1.0, double addFactor = 0.0, bool isBands = false): units(units), name(name), standardName(standardName), longName(longName), cellMethods(cellMethods), scalingFactor(scalingFactor), addFactor(addFactor), isBands(isBands) {}
  string units, name, standardName, longName, cellMethods;
  double scalingFactor;
  double addFactor;
  bool isBands;
};

class StateIO {
public:
  StateIO(const ProgramState* state);
  virtual ~StateIO();
  virtual void initializeOutput(FILE** f, const char* filename, const ProgramState* state) = 0;
  virtual int write(const int* data, int numValues, const StateVariableMetaData* meta) = 0;
  virtual int write(const double* data, int numValues, const StateVariableMetaData* meta) = 0;
  virtual int write(const char* data, int numValues, const StateVariableMetaData* meta) = 0;
  virtual int writeNewline();
  virtual int read(int* data, int numValues) = 0;
  virtual int read(double* data, int numValues) = 0;
  virtual void flush() = 0;
private:
  const ProgramState* state;
};

#endif /* STATEIO_H_ */
