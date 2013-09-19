#include "root_brent.h"

#ifndef ICEENERGYBALANCE_H_
#define ICEENERGYBALANCE_H_

class IceEnergyBalance : public RootBrent {
public:
  IceEnergyBalance( double Dt,                  /* Model time step (hours) */
                    double Ra,                  /* Aerodynamic resistance (s/m) */
                    double* Ra_used,            /* Aerodynamic resistance (s/m) after stability correction */
                    double Z,                   /* Reference height (m) */
                    double Displacement,        /* Displacement height (m) */
                    double Z0,                  /* surface roughness height (m) */
                    double Wind,                /* Wind speed (m/s) */
                    double ShortRad,            /* Net incident shortwave radiation (W/m2) */
                    double LongRadIn,           /* Incoming longwave radiation (W/m2) */
                    double AirDens,             /* Density of air (kg/m3) */
                    double Lv,                  /* Latent heat of vaporization (J/kg3) */
                    double Tair,                /* Air temperature (C) */
                    double Press,               /* Air pressure (Pa) */
                    double Vpd,                 /* Vapor pressure deficit (Pa) */
                    double EactAir,             /* Actual vapor pressure of air (Pa) */
                    double Rain,                /* Rain fall (m/timestep) */
                    double SweSurfaceLayer,     /* Snow water equivalent in surface layer (m) */
                    double SurfaceLiquidWater,  /* Liquid water in the surface layer (m) */
                    double OldTSurf,            /* Surface temperature during previous time step */
                    double* RefreezeEnergy,     /* Refreeze energy (W/m2) */
                    double* vapor_flux,         /* Total mass flux of water vapor to or from snow (m/timestep) */
                    double* blowing_flux,       /* Mass flux of water vapor to or from blowing snow (m/timestep) */
                    double* surface_flux,       /* Mass flux of water vapor to or from snow pack (m/timestep) */
                    double* AdvectedEnergy,     /* Energy advected by precipitation (W/m2) */
                    double DeltaColdContent,    /* Change in cold content (W/m2) */
                    double Tfreeze,
                    double AvgCond,
                    double SWconducted,
                    double SnowDepth,
                    double SnowDensity,
                    double SurfAttenuation,
                    double* qf,                  /* Ground Heat Flux (W/m2) */
                    double* LatentHeat,          /* Latent heat exchange at surface (W/m2) */
                    double* LatentHeatSub,       /* Latent heat exchange at surface (W/m2) due to sublimation */
                    double* SensibleHeat,        /* Sensible heat exchange at surface (W/m2) */
                    double* LongRadOut) :
      Dt(Dt), Ra(Ra), Ra_used(Ra_used), Z(Z), Displacement(Displacement), Z0(Z0), Wind(Wind), ShortRad(ShortRad),
      LongRadIn(LongRadIn), AirDens(AirDens), Lv(Lv), Tair(Tair), Press(Press), Vpd(Vpd), EactAir(EactAir), Rain(Rain),
      SweSurfaceLayer(SweSurfaceLayer), SurfaceLiquidWater(SurfaceLiquidWater), OldTSurf(OldTSurf), RefreezeEnergy(RefreezeEnergy),
      vapor_flux(vapor_flux), blowing_flux(blowing_flux), surface_flux(surface_flux), AdvectedEnergy(AdvectedEnergy),
      DeltaColdContent(DeltaColdContent), Tfreeze(Tfreeze), AvgCond(AvgCond), SWconducted(SWconducted), SnowDepth(SnowDepth),
      SnowDensity(SnowDensity), SurfAttenuation(SurfAttenuation), qf(qf), LatentHeat(LatentHeat), LatentHeatSub(LatentHeatSub),
      SensibleHeat(SensibleHeat), LongRadOut(LongRadOut) {}

  double calculate(double);
private:
  double Dt;
  double Ra;
  double* Ra_used;
  double Z;
  double Displacement;
  double Z0;
  double Wind;
  double ShortRad;
  double LongRadIn;
  double AirDens;
  double Lv;
  double Tair;
  double Press;
  double Vpd;
  double EactAir;
  double Rain;
  double SweSurfaceLayer;
  double SurfaceLiquidWater;
  double OldTSurf;
  double* RefreezeEnergy;
  double* vapor_flux;
  double* blowing_flux;
  double* surface_flux;
  double* AdvectedEnergy;
  double DeltaColdContent;
  double Tfreeze;
  double AvgCond;
  double SWconducted;
  double SnowDepth;
  double SnowDensity;
  double SurfAttenuation;
  double* qf;
  double* LatentHeat;
  double* LatentHeatSub;
  double* SensibleHeat;
  double* LongRadOut;
};

#endif /* ICEENERGYBALANCE_H_ */








































