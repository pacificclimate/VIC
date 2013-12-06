#include "root_brent.h"
#include "VegConditions.h"

#ifndef SNOWPACKENERGYBALANCE_H_
#define SNOWPACKENERGYBALANCE_H_

class SnowPackEnergyBalance : public RootBrent {
public:
  SnowPackEnergyBalance(  double Dt,      /* Model time step (sec) */
      double Ra,                          /* Aerodynamic resistance (s/m) */
      AeroResistUsed& Ra_used,            /* Aerodynamic resistance (s/m) after stability correction */
      double Displacement,                /* Displacement height (m) */
      double Z,                           /* Reference height (m) */
      VegConditions& roughness,           /* surface roughness height (m) */
      double AirDens,                     /* Density of air (kg/m3) */
      double EactAir,                     /* Actual vapor pressure of air (Pa) */
      double LongSnowIn,                  /* Incoming longwave radiation (W/m2) */
      double Lv,                          /* Latent heat of vaporization (J/kg3) */
      double Press,                       /* Air pressure (Pa) */
      double Rain,                        /* Rain fall (m/timestep) */
      double NetShortUnder,               /* Net incident shortwave radiation  (W/m2) */
      double Vpd,                         /* Vapor pressure deficit (Pa) */
      double Wind,                        /* Wind speed (m/s) */
      double OldTSurf,                    /* Surface temperature during previous time step */
      double SnowCoverFract,              /* Fraction of area covered by snow */
      double SnowDepth,                   /* Depth of snowpack (m) */
      double SnowDensity,                 /* Density of snowpack (kg/m^3) */
      double SurfaceLiquidWater,          /* Liquid water in the surface layer (m) */
      double SweSurfaceLayer,             /* Snow water equivalent in surface layer  (m) */
      double Tair,                        /* Canopy air / Air temperature (C) */
      double TGrnd,                       /* Ground surface temperature (C) */
      double* AdvectedEnergy,             /* Energy advected by precipitation (W/m2) */
      double* AdvectedSensibleHeat,       /* Sensible heat advected from snow-free area into snow covered area (W/m^2) */
      double* DeltaColdContent,           /* Change in cold content of surface  layer (W/m2) */
      double* GroundFlux,                 /* Ground Heat Flux (W/m2) */
      double* LatentHeat,                 /* Latent heat exchange at surface (W/m2) */
      double* LatentHeatSub,              /* Latent heat of sublimation exchange at  surface (W/m2) */
      double* NetLongUnder,               /* Net longwave radiation at snowpack  surface (W/m^2) */
      double* RefreezeEnergy,             /* Refreeze energy (W/m2) */
      double* SensibleHeat,               /* Sensible heat exchange at surface (W/m2) */
      double* vapor_flux,                 /* Mass flux of water vapor to or from the intercepted snow (m/timestep) */
      double* blowing_flux,               /* Mass flux of water vapor from blowing snow. (m/timestep) */
      double* surface_flux                /* Mass flux of water vapor from pack snow. (m/timestep) */
  ) :
      Dt(Dt), Ra(Ra), Ra_used(Ra_used), Displacement(Displacement), Z(Z), roughness(roughness), AirDens(AirDens),
      EactAir(EactAir), LongSnowIn(LongSnowIn), Lv(Lv), Press(Press), Rain(Rain), NetShortUnder(NetShortUnder),
      Vpd(Vpd), Wind(Wind), OldTSurf(OldTSurf), SnowCoverFract(SnowCoverFract), SnowDepth(SnowDepth),
      SnowDensity(SnowDensity), SurfaceLiquidWater(SurfaceLiquidWater), SweSurfaceLayer(SweSurfaceLayer), Tair(Tair),
      TGrnd(TGrnd), AdvectedEnergy(AdvectedEnergy), AdvectedSensibleHeat(AdvectedSensibleHeat), DeltaColdContent(DeltaColdContent),
      GroundFlux(GroundFlux), LatentHeat(LatentHeat), LatentHeatSub(LatentHeatSub), NetLongUnder(NetLongUnder),
      RefreezeEnergy(RefreezeEnergy), SensibleHeat(SensibleHeat), vapor_flux(vapor_flux), blowing_flux(blowing_flux), surface_flux(surface_flux) {
  }

  double calculate(double);
private:

  double Dt;
  double Ra;
  AeroResistUsed& Ra_used;
  double Displacement;
  double Z;
  VegConditions& roughness;
  double AirDens;
  double EactAir;
  double LongSnowIn;
  double Lv;
  double Press;
  double Rain;
  double NetShortUnder;
  double Vpd;
  double Wind;
  double OldTSurf;
  double SnowCoverFract;
  double SnowDepth;
  double SnowDensity;
  double SurfaceLiquidWater;
  double SweSurfaceLayer;
  double Tair;
  double TGrnd;
  double* AdvectedEnergy;
  double* AdvectedSensibleHeat;
  double* DeltaColdContent;
  double* GroundFlux;
  double* LatentHeat;
  double* LatentHeatSub;
  double* NetLongUnder;
  double* RefreezeEnergy;
  double* SensibleHeat;
  double* vapor_flux;
  double* blowing_flux;
  double* surface_flux;
};

#endif /* SNOWPACKENERGYBALANCE_H_ */






































