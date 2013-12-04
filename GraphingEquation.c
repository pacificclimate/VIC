
#include "GraphingEquation.h"
#include "vicNl.h"

#include <cstdio>
#include <cmath>

void GraphingEquation::calculateFitError(const std::vector<GraphPoint>& points) {
  fitError = 0;
  for (std::vector<GraphPoint>::const_iterator point = points.begin(); point != points.end(); ++point) {
    double curveValue = b0 + b1*(point->elevationX) + b2*(point->elevationX * point->elevationX);
    fitError += std::abs(curveValue - point->massBalanceY);
  }
  fprintf(stderr, "Cumulative fit error: %lf\n", fitError);
}

void GraphingEquation::calculateEquationFrom2Points(GraphPoint p1, GraphPoint p2) {
  double slope = (p2.massBalanceY - p1.massBalanceY) / (p2.elevationX - p1.elevationX);
  b0 = p1.massBalanceY - slope * p1.elevationX; // Solve for b0 in the equation y = mx + b0. b0 = y - mx. Use p1 as x and y.
  b1 = slope;
  b2 = 0;
}

void GraphingEquation::calculateEquationFrom1Point(GraphPoint p1) {
  // This is arbitrarily defined as a horizontal line through this point.
  b0 = p1.massBalanceY;
  b1 = 0;
  b2 = 0;
}

// The matrix algorithm here was adapted from a stack overflow question: http://stackoverflow.com/questions/11233404/simple-curve-fitting-implimentation-in-c-svd-least-sqares-fit-or-similar
// Where one of the users posted his working code at: http://pastebin.com/tUvKmGPn (note: his code has memory leaks which have been fixed here)
void GraphingEquation::calculateEquationFromNPoints(const std::vector<GraphPoint>& points) {

  if (points.size() == 0) {
    throw VICException("Error while calculating GlacierMassBalanceResult. Cannot calculate an equation from no graph points!");
  } else if (points.size() == 1) {
    calculateEquationFrom1Point(points[0]);
    return;
  } else if (points.size() == 2) {
    calculateEquationFrom2Points(points[0], points[1]);
    return;
  }

  //determine how many points used
  int size = points.size();
  //Allocate memory for data[size][2]
  //data[3][0] -> x4, data[3][1] -> y4
  std::vector<GraphPoint> data;
  for (int i = 0; i < size; i++) {
    data.push_back(GraphPoint());
  }

  // sumx3 = (x0)^3 + (x1)^3 + ... + (x{size})^3
  double sumx4 = 0, sumx3 = 0, sumx2 = 0, sumx1 = 0, det;

  // Allocate memory for X^T
  std::vector<std::vector<double> > transpose;
  std::vector<double> filler;
  for (int i = 0; i < size; i++) {
    filler.push_back(0);
  }
  for (int i = 0; i < 3; i++)
    transpose.push_back(filler);

  //Allocate memory for (X^T *X)^-1 * X^T
  std::vector<std::vector<double> > stuff;
  for (int i = 0; i < 3; i++)
    stuff.push_back(filler);

  // y = a[0]x^2 + a[1]x + a[2]
  double a[3] = {0};

  //input data, do some precalculations at the same time to avoid making more loops
  for (int i = 0; i < size; i++) {
    data[i].elevationX = points[i].elevationX;
    data[i].massBalanceY = points[i].massBalanceY;

    for (int j = 0; j < 2; j++) {
      //these computations only happen with x values
      if (j == 0) {
        //develops the sum values needed for inverse
        sumx4 += data[i].elevationX * data[i].elevationX * data[i].elevationX * data[i].elevationX;
        sumx3 += data[i].elevationX * data[i].elevationX * data[i].elevationX;
        sumx2 += data[i].elevationX * data[i].elevationX;
        sumx1 += data[i].elevationX;

        //develops transpose matrix
        transpose[2][i] = 1;
        transpose[1][i] = data[i].elevationX;
        transpose[0][i] = data[i].elevationX * data[i].elevationX;
      }
    }
  }

  //After solving all the math
  //determinate
  det = (sumx4*sumx2*size) + (sumx3*sumx1*sumx2) + (sumx2*sumx3*sumx1) - (sumx2*sumx2*sumx2) - (sumx1*sumx1*sumx4) - (size*sumx3*sumx3);

  //precalculated the inverse matrix to avoid numerical methods which take time or lose accuracy, NOTE: does not include division of determinate
  double inverse[3][3] = {
          {size*sumx2 - sumx1*sumx1, -(size*sumx3 - sumx1*sumx2), sumx1*sumx3-sumx2*sumx2},
          {-(size*sumx3-sumx2*sumx1), size*sumx4-sumx2*sumx2, -(sumx1*sumx4-sumx3*sumx2)},
          {sumx1*sumx3 - sumx2*sumx2, -(sumx1*sumx4 - sumx2*sumx3), sumx2*sumx4 - sumx3*sumx3}
  };

  //This is matrix multiplication for this particular pair of matrices
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < size; j++) {
      stuff[i][j] = inverse[i][0] * transpose[0][j] + inverse[i][1] * transpose[1][j] + inverse[i][2] * transpose[2][j];
    }
  }

  //This is the final matrix multiplication that outputs a 1x3 matrix with our curve parameters
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < size; j++) {
      a[i] += stuff[i][j] * data[j].massBalanceY;
    }
    //don't forget to divide by determinate
    a[i] /= det;
  }

  b0 = a[2];
  b1 = a[1];
  b2 = a[0];
}
