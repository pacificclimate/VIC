
#include <cstdio>
#include <sstream>
#include "GlacierMassBalanceResult.h"

// Test function, used for development purposes only.
void testGlacierMassBalanceResult() {
  soil_con_struct soil;
  soil.BandElev = new float[MAX_BANDS];
  const int NUM_DATA_POINTS = 6;
  double xValues [] = { 0.2, 3, 5, 9, 14, 17, 22 };
  double yValues [] = { -6.6, -3.2, 0, 4, 5.9, 6 };
  for (int i = 0; i < NUM_DATA_POINTS; i++) {
    soil.BandElev[i] = xValues[i];
  }
  std::vector<HRU> hruList;

  for (int i = 0; i < NUM_DATA_POINTS; i++) {
    HRU hru;
    hru.isGlacier = true;
    hru.bandIndex = i;
    hru.glacier.cum_mass_balance = yValues[i];
    hruList.push_back(hru);
  }
  dmy_struct dmy;
  dmy.day = 1;
  dmy.month = 1;
  dmy.year = 2006;

  GlacierMassBalanceResult result(hruList, &soil, dmy);
  result.printForDebug();
}

GlacierMassBalanceResult::GlacierMassBalanceResult(const std::vector<HRU>& hruList, const soil_con_struct* soil, const dmy_struct& dmy) {

  date = dmy;

  for (std::vector<HRU>::const_iterator hru = hruList.begin(); hru != hruList.end(); ++hru) {
    if (hru->isGlacier) {
      if (IS_VALID(hru->glacier.cum_mass_balance)) {
        // Check if there is already an entry for this elevation.
        bool found = false;
        for (std::vector<GraphPoint>::iterator it = graphPoints.begin(); it != graphPoints.end(); ++it) {
          if (it->elevationX == soil->BandElev[hru->bandIndex]) {
            it->massBalanceY += hru->glacier.cum_mass_balance;
            found = true;
          }
        }
        if (!found) {
          GraphPoint p;
          p.elevationX = soil->BandElev[hru->bandIndex];
          p.massBalanceY = hru->glacier.cum_mass_balance;
          graphPoints.push_back(p);
        }
      }
    }
  }

  // Remove graph points which are meaningless (no data for this cell in these bands)
  double lastElevation = 0;
  for (std::vector<GraphPoint>::iterator point = graphPoints.begin(); point != graphPoints.end(); /* nothing */) {
    if (point->elevationX == 0 && point->elevationX <= lastElevation) {
      point = graphPoints.erase(point);
    } else {
      ++point;
    }
  }

  if (graphPoints.size() > 0) {
    calculateEquationFromPoints();
  }
}

void GlacierMassBalanceResult::calculateEquationFromPoints() {
  equation.calculateEquationFromNPoints(graphPoints);
  equation.calculateFitError(graphPoints);
}

// The output from printing (graph points list and function equation) can be directly copied to fooplot.com
void GlacierMassBalanceResult::printForDebug() {
  std::stringstream output;
  output << date.year << "\t" << date.month << "\t" << date.day << "\n";
  output << "    best fit equation: " << equation.b0 << " + " << equation.b1 << "*x + " << equation.b2 << "*x^2\n";
  output << "    accumulated total error of the curve: " << equation.fitError << "\n";
  output << "    elevation, mass balance points used for fit:";
  for (std::vector<GraphPoint>::iterator dataPoint = graphPoints.begin(); dataPoint != graphPoints.end(); ++dataPoint) {
    output << "\n        " << dataPoint->elevationX << "," << dataPoint->massBalanceY;
  }

  // Only print results if they mean anything (if there are actually data points to consider.
  if (graphPoints.size() > 0) {
    fprintf(stderr, "Output from GlacierMassBalanceResult for: %s\n", output.str().c_str());
  }
}
