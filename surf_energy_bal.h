#include "root_brent.h"
#include "VegConditions.h"

#ifndef SURF_ENERGY_BAL_H_
#define SURF_ENERGY_BAL_H_

class SurfEnergyBal : public RootBrent {
public:
  SurfEnergyBal(int rec, int nrecs, int month, int VEG, int veg_class,
      double delta_t, double Cs1, double Cs2, double D1, double D2,
      double T1_old, double T2, double Ts_old, double bubble, double dp,
      double expt, double ice0, double kappa1, double kappa2, double max_moist,
      double moist, const float* root, VegConditions::VegSurfType UnderStory, int overstory,
      double NetShortBare, double NetShortGrnd, double NetShortSnow,
      double Tair, double atmos_density, double atmos_pressure,
      double emissivity, double LongBareIn, double LongSnowIn,
      double precipitation_mu, double surf_atten, double vp, double vpd,
      double* Wdew, VegConditions& displacement, VegConditions& aero_resist, AeroResistUsed& aero_resist_used,
      double* rainfall, VegConditions& ref_height, VegConditions& roughness, VegConditions& wind_speed,
      double latent_heat_Le, double Advection, double OldTSurf, double TPack,
      double Tsnow_surf, double kappa_snow, double melt_energy,
      double snow_coverage, double snow_density, double snow_swq,
      double snow_water, double* deltaCC, double* refreeze_energy,
      double* vapor_flux, double* blowing_flux, double* surface_flux,
      int Nnodes, double* Cs_node, double* T_node, double* Tnew_node,
      char* Tnew_fbflag, int* Tnew_fbcount,
      double* ice_node, double* kappa_node,
      double* moist_node, const soil_con_struct* soil_con,
      layer_data_struct* layer_wet, layer_data_struct* layer_dry,
      veg_var_struct* veg_var_wet, veg_var_struct* veg_var_dry,
      int INCLUDE_SNOW, int NOFLUX, int EXP_TRANS, int SNOWING, int* FIRST_SOLN,
      double* NetLongBare, double* NetLongSnow, double* T1, double* deltaH,
      double* fusion, double* grnd_flux, double* latent_heat,
      double* latent_heat_sub, double* sensible_heat, double* snow_flux,
      double* store_error, const ProgramState* state) :

      rec(rec), nrecs(nrecs), month(month), VEG(VEG), veg_class(veg_class),
          delta_t(delta_t), Cs1(Cs1), Cs2(Cs2), D1(D1), D2(D2), T1_old(
          T1_old), T2(T2), Ts_old(Ts_old), bubble(bubble), dp(dp), expt(expt), ice0(
          ice0), kappa1(kappa1), kappa2(kappa2), max_moist(max_moist), moist(
          moist), root(root), UnderStory(UnderStory), overstory(overstory), NetShortBare(
          NetShortBare), NetShortGrnd(NetShortGrnd), NetShortSnow(NetShortSnow), Tair(
          Tair), atmos_density(atmos_density), atmos_pressure(atmos_pressure), emissivity(
          emissivity), LongBareIn(LongBareIn), LongSnowIn(LongSnowIn), precipitation_mu(
          precipitation_mu), surf_atten(surf_atten), vp(vp), vpd(vpd), Wdew(
          Wdew), displacement(displacement), aero_resist(aero_resist), aero_resist_used(aero_resist_used), rainfall(
          rainfall), ref_height(ref_height), roughness(roughness), wind_speed(wind_speed), latent_heat_Le(
          latent_heat_Le), Advection(Advection), OldTSurf(OldTSurf), TPack(
          TPack), Tsnow_surf(Tsnow_surf), kappa_snow(kappa_snow), melt_energy(
          melt_energy), snow_coverage(snow_coverage), snow_density(
          snow_density), snow_swq(snow_swq), snow_water(snow_water), deltaCC(
          deltaCC), refreeze_energy(refreeze_energy), vapor_flux(vapor_flux), blowing_flux(
          blowing_flux), surface_flux(surface_flux), Nnodes(Nnodes), Cs_node(
          Cs_node), T_node(T_node), Tnew_node(Tnew_node), Tnew_fbflag(
          Tnew_fbflag), Tnew_fbcount(Tnew_fbcount), ice_node(ice_node), kappa_node(
          kappa_node), moist_node(moist_node), soil_con(soil_con), layer_wet(
          layer_wet), layer_dry(layer_dry), veg_var_wet(veg_var_wet), veg_var_dry(
          veg_var_dry), INCLUDE_SNOW(INCLUDE_SNOW), NOFLUX(NOFLUX), EXP_TRANS(
          EXP_TRANS), SNOWING(SNOWING), FIRST_SOLN(FIRST_SOLN), NetLongBare(
          NetLongBare), NetLongSnow(NetLongSnow), T1(T1), deltaH(deltaH), fusion(
          fusion), grnd_flux(grnd_flux), latent_heat(latent_heat), latent_heat_sub(
          latent_heat_sub), sensible_heat(sensible_heat), snow_flux(snow_flux), store_error(
          store_error), state(state), error_cnt0(0), error_cnt1(0) {
  }

  double calculate(double);


  int rec;
  int nrecs;
  int month;
  int VEG;
  int veg_class;
  double delta_t;
  double Cs1;
  double Cs2;
  double D1;
  double D2;
  double T1_old;
  double T2;
  double Ts_old;
  double bubble;
  double dp;
  double expt;
  double ice0;
  double kappa1;
  double kappa2;
  double max_moist;
  double moist;
  const float* root;
  VegConditions::VegSurfType UnderStory;
  int overstory;
  double NetShortBare;
  double NetShortGrnd;
  double NetShortSnow;
  double Tair;
  double atmos_density;
  double atmos_pressure;
  double emissivity;
  double LongBareIn;
  double LongSnowIn;
  double precipitation_mu;
  double surf_atten;
  double vp;
  double vpd;
  double* Wdew;
  VegConditions& displacement;
  VegConditions& aero_resist;
  AeroResistUsed& aero_resist_used;
  double* rainfall;
  VegConditions& ref_height;
  VegConditions& roughness;
  VegConditions& wind_speed;
  double latent_heat_Le;
  double Advection;
  double OldTSurf;
  double TPack;
  double Tsnow_surf;
  double kappa_snow;
  double melt_energy;
  double snow_coverage;
  double snow_density;
  double snow_swq;
  double snow_water;
  double* deltaCC;
  double* refreeze_energy;
  double* vapor_flux;
  double* blowing_flux;
  double* surface_flux;
  int Nnodes;
  double* Cs_node;
  double* T_node;
  double* Tnew_node;
  char* Tnew_fbflag;
  int* Tnew_fbcount;
  double* ice_node;
  double* kappa_node;
  double* moist_node;
  const soil_con_struct* soil_con;
  layer_data_struct* layer_wet;
  layer_data_struct* layer_dry;
  veg_var_struct* veg_var_wet;
  veg_var_struct* veg_var_dry;
  int INCLUDE_SNOW;
  int NOFLUX;
  int EXP_TRANS;
  int SNOWING;
  int* FIRST_SOLN;
  double* NetLongBare;
  double* NetLongSnow;
  double* T1;
  double* deltaH;
  double* fusion;
  double* grnd_flux;
  double* latent_heat;
  double* latent_heat_sub;
  double* sensible_heat;
  double* snow_flux;
  double* store_error;
  const ProgramState* state;
  //error counting variables for IMPLICIT option
  int error_cnt0;
  int error_cnt1;

};

#endif /* SURF_ENERGY_BAL_H_ */
