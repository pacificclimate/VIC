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

class StateHeader {
public:
  StateHeader(int year, int month, int day, int nLayer, int nNode) : year(year), month(month), day(day), nLayer(nLayer), nNode(nNode) {}
  int year;
  int month;
  int day;
  int nLayer;
  int nNode;
  bool isValid() { return (year > 0 && month > 0 && day > 0); }
};

class StateIO {
public:
  enum IOType { Reader, Writer };
  StateIO(std::string filename, IOType type, const ProgramState* state);
  virtual ~StateIO();
  virtual void initializeOutput() = 0;
  virtual int write(const int* data, int numValues, const StateVariableMetaData* meta) = 0;
  virtual int write(const double* data, int numValues, const StateVariableMetaData* meta) = 0;
  virtual int write(const char* data, int numValues, const StateVariableMetaData* meta) = 0;
  virtual int writeNewline();
  virtual StateHeader readHeader() = 0;
  virtual int seekToCell(int cellid, int* nVeg, int* nBand) = 0;
  virtual int read(int* data, int numValues, const StateVariableMetaData* meta) = 0;
  virtual int read(double* data, int numValues, const StateVariableMetaData* meta) = 0;
  virtual int read(char* data, int numValues, const StateVariableMetaData* meta) = 0;
  virtual void flush() = 0;
  virtual void rewindFile() = 0;
protected:
  std::string filename;
  const ProgramState* state;
  const IOType ioType;
};

#endif /* STATEIO_H_ */
