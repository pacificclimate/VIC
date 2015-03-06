#ifndef GLACIERMASSBALANCERESULT_H_
#define GLACIERMASSBALANCERESULT_H_

#include <vector>
#include "GraphingEquation.h"
#include "vicNl_def.h"

struct GlacierMassBalanceResult {
  GlacierMassBalanceResult(const std::vector<HRU>& hruList, const soil_con_struct* soil, const dmy_struct& dmy);
  void printForDebug();
  void calculateEquationFromPoints();
  void calculateFitError();
  GraphingEquation equation;
  dmy_struct date;
  std::vector<GraphPoint> graphPoints;
};


#endif /* GLACIERMASSBALANCERESULT_H_ */
