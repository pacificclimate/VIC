#include "root_brent.h"

#ifndef CANOPY_ENERGY_BAL_H_
#define CANOPY_ENERGY_BAL_H_

class CanopyEnergyBal : public RootBrent {
public:
  CanopyEnergyBal(int band, int month, int rec, double delta_t,
      const double elevation, const double* Wcr, const double* Wpwp, const double* depth,
      const double* frost_fract, double AirDens, double EactAir, double Press,
      double latent_heat_Le, double Tcanopy, double Vpd,
      double precipitation_mu, double* Evap, VegConditions& Ra, AeroResistUsed& Ra_used,
      double* Rainfall, VegConditions& wind_speed, VegConditions::VegetationConditions UnderStory, int iveg, int veg_class,
      VegConditions& displacement, VegConditions& ref_height, VegConditions& roughness, const float* root,
      double IntRain, double IntSnow, double* Wdew,
      layer_data_struct* layer_wet, layer_data_struct * layer_dry,
      veg_var_struct * veg_var_wet, veg_var_struct * veg_var_dry,
      double LongOverIn, double LongUnderOut, double NetShortOver,
      double* AdvectedEnergy, double* LatentHeat, double* LatentHeatSub,
      double* LongOverOut, double* NetLongOver, double* NetRadiation,
      double* RefreezeEnergy, double* SensibleHeat, double* VaporMassFlux,
      const ProgramState* state) :
      band(band), month(month), rec(rec), delta_t(delta_t), elevation(elevation), Wcr(Wcr), Wpwp(Wpwp),
      depth(depth), frost_fract(frost_fract), AirDens(AirDens), EactAir(EactAir), Press(Press), latent_heat_Le(latent_heat_Le),
      Tcanopy(Tcanopy), Vpd(Vpd), precipitation_mu(precipitation_mu), Evap(Evap), Ra(Ra), Ra_used(Ra_used),
      Rainfall(Rainfall), wind_speed(wind_speed), UnderStory(UnderStory), iveg(iveg), veg_class(veg_class), displacement(displacement),
      ref_height(ref_height), roughness(roughness), root(root), IntRain(IntRain), IntSnow(IntSnow), Wdew(Wdew), layer_wet(layer_wet),
      layer_dry(layer_dry), veg_var_wet(veg_var_wet), veg_var_dry(veg_var_dry), LongOverIn(LongOverIn), LongUnderOut(LongUnderOut),
      NetShortOver(NetShortOver), AdvectedEnergy(AdvectedEnergy), LatentHeat(LatentHeat), LatentHeatSub(LatentHeatSub),
      LongOverOut(LongOverOut), NetLongOver(NetLongOver), NetRadiation(NetRadiation), RefreezeEnergy(RefreezeEnergy),
      SensibleHeat(SensibleHeat), VaporMassFlux(VaporMassFlux), state(state) {
  }

  double calculate(double);
private:
  int band;
  int month;
  int rec;
  double delta_t;
  const double elevation;
  const double* Wcr;
  const double* Wpwp;
  const double* depth;
  const double* frost_fract;
  double AirDens;
  double EactAir;
  double Press;
  double latent_heat_Le;
  double Tcanopy;
  double Vpd;
  double precipitation_mu;
  double* Evap;
  VegConditions& Ra;
  AeroResistUsed& Ra_used;
  double* Rainfall;
  VegConditions& wind_speed;
  VegConditions::VegetationConditions UnderStory;
  int iveg;
  int veg_class;
  VegConditions& displacement;
  VegConditions& ref_height;
  VegConditions& roughness;
  const float* root;
  double IntRain;
  double IntSnow;
  double* Wdew;
  layer_data_struct* layer_wet;
  layer_data_struct * layer_dry;
  veg_var_struct * veg_var_wet;
  veg_var_struct * veg_var_dry;
  double LongOverIn;
  double LongUnderOut;
  double NetShortOver;
  double* AdvectedEnergy;
  double* LatentHeat;
  double* LatentHeatSub;
  double* LongOverOut;
  double* NetLongOver;
  double* NetRadiation;
  double* RefreezeEnergy;
  double* SensibleHeat;
  double* VaporMassFlux;
  const ProgramState* state;


};


#endif /* CANOPY_ENERGY_BAL_H_ */
