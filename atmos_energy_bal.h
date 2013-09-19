#include "root_brent.h"

#ifndef ATMOS_ENERGY_BAL_H_
#define ATMOS_ENERGY_BAL_H_

class AtmosEnergyBal : public RootBrent {
public:
  AtmosEnergyBal(double LatentHeat, double NetRadiation, double Ra, double Tair,
      double atmos_density, double InSensible, double *SensibleHeat) :
      LatentHeat(LatentHeat), NetRadiation(NetRadiation), Ra(Ra), Tair(Tair), atmos_density(
          atmos_density), InSensible(InSensible), SensibleHeat(SensibleHeat) {
  }
  double calculate(double Tcanopy);

private:
  double  LatentHeat;
  double  NetRadiation;
  double  Ra;
  double  Tair;
  double  atmos_density;
  double  InSensible;
  double *SensibleHeat;
};


#endif /* ATMOS_ENERGY_BAL_H_ */
