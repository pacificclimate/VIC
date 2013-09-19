#include "root_brent.h"

#ifndef SOIL_THERMAL_EQN_H_
#define SOIL_THERMAL_EQN_H_

class SoilThermalEqn : public RootBrent {
public:
  SoilThermalEqn(double TL, double TU, double T0, double moist,
      double max_moist, double** ufwc_table, double bubble, double expt,
      double porosity, double effective_porosity, double ice0, double gamma,
      double A, double B, double C, double D, double E, int EXP_TRANS,
      int node) :
      TL(TL), TU(TU), T0(T0), moist(moist), max_moist(max_moist), ufwc_table(ufwc_table), bubble(bubble),
      expt(expt), porosity(porosity), effective_porosity(effective_porosity), ice0(ice0), gamma(gamma), A(A),
      B(B), C(C), D(D), E(E), EXP_TRANS(EXP_TRANS), node(node) {
  }

  double calculate(double);
private:
  double TL;
  double TU;
  double T0;
  double moist;
  double max_moist;
  double** ufwc_table;
  double bubble;
  double expt;
  double porosity;
  double effective_porosity;
  double ice0;
  double gamma;
  double A;
  double B;
  double C;
  double D;
  double E;
  int EXP_TRANS;
  int node;
};


#endif /* SOIL_THERMAL_EQN_H_ */
