
#include <cstdio>
#include <sstream>
#include "GlacierMassBalanceResult.h"


GlacierMassBalanceResult::GlacierMassBalanceResult(const std::vector<HRU>& hruList, const dmy_struct& dmy) {

  date = dmy;

  for (std::vector<HRU>::const_iterator hru = hruList.begin(); hru != hruList.end(); ++hru) {
    if (hru->isGlacier) {
      //TODO: is bandIndex sufficient or do we need actual elevation?
      if (IS_VALID(hru->glacier.cum_mass_balance)) {
        // Check if there is already an entry for this elevation.
        bool found = false;
        for (std::vector<std::pair<double,double> >::iterator it = graphPoints.begin(); it != graphPoints.end(); ++it) {
          if (it->second == hru->bandIndex) {
            it->first += hru->glacier.cum_mass_balance;
            found = true;
          }
        }
        if (!found) {
          graphPoints.push_back(std::pair<double, double>(hru->glacier.cum_mass_balance, (double)hru->bandIndex));
        }
      }
    }
  }

  calculateEquationFromPoints();
}

void GlacierMassBalanceResult::calculateEquationFromPoints() {
  //TODO: actually calculate equation parameters.
  b0 = 0;
  b1 = 0;
  b2 = 0;
  fitError = -1;
}

void GlacierMassBalanceResult::printForDebug() {

  std::stringstream output;
  output << date.year << "," << date.month << "," << date.day << "," << b0 << "," << b1 << "," << b2;
  for (std::vector<std::pair<double, double> >::iterator dataPoint = graphPoints.begin(); dataPoint != graphPoints.end(); ++dataPoint) {
    output << "," << "(" << dataPoint->first << "," << dataPoint->second << ")";
  }

  // Only print results if they mean anything (if there are actually data points to consider.
  if (graphPoints.size() > 0) {
    fprintf(stderr, "Output from GlacierMassBalanceResult: %s\n", output.str().c_str());
  }
}
