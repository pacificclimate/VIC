#ifndef GRAPHINGEQUATION_H_
#define GRAPHINGEQUATION_H_

#include <vector>

struct GraphPoint {
  GraphPoint() : elevationX(0), massBalanceY(0) {}
  double elevationX;
  double massBalanceY;
};

struct GraphingEquation {
  GraphingEquation() : b0(0), b1(0), b2(0), fitError(-1) {}
  void calculateFitError(const std::vector<GraphPoint>& points);
  void calculateEquationFrom2Points(GraphPoint p1, GraphPoint p2);
  void calculateEquationFrom1Point(GraphPoint p1);
  void calculateEquationFromNPoints(const std::vector<GraphPoint>& points);
  double b0;  // x^0 coefficient
  double b1;  // x^1 coefficient
  double b2;  // x^2 coefficient
  double fitError;
};

#endif /* GRAPHINGEQUATION_H_ */
