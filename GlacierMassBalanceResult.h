#ifndef GLACIERMASSBALANCERESULT_H_
#define GLACIERMASSBALANCERESULT_H_

#include <vector>
#include "vicNl_def.h"

struct GlacierMassBalanceResult {
  GlacierMassBalanceResult(const std::vector<HRU>& hruList, const dmy_struct& dmy);
  void printForDebug();
private:
  void calculateEquationFromPoints();
  double b0;
  double b1;
  double b2;
  double fitError;
  dmy_struct date;
  std::vector< std::pair<double, double> > graphPoints;
};



#endif /* GLACIERMASSBALANCERESULT_H_ */
