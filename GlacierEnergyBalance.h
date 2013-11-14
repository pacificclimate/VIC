#ifndef GLACIERENERGYBALANCE_H_
#define GLACIERENERGYBALANCE_H_

#include "root_brent.h"

class GlacierEnergyBalance : public RootBrent {
public:
  GlacierEnergyBalance(
  double Dt,                      /* Model time step (sec) */
  double Ra,                      /* Aerodynamic resistance (s/m) */
  double *Ra_used,                /* Aerodynamic resistance (s/m) after stability correction */
  /* Vegetation Parameters */
  double Displacement,            /* Displacement height (m) */
  double Z,                       /* Reference height (m) */
  double *Z0,                     /* surface roughness height (m) */
  /* Atmospheric Forcing Variables */
  double AirDens,                 /* Density of air (kg/m3) */
  double EactAir,                 /* Actual vapor pressure of air (Pa) */
  double LongSnowIn,              /* Incoming longwave radiation (W/m2) */
  double Lv,                      /* Latent heat of vaporization (J/kg3) */
  double Press,                   /* Air pressure (Pa) */
  double Rain,                    /* Rain fall (m/timestep) */
  double NetShortUnder,           /* Net incident shortwave radiation (W/m2) */
  double Vpd,                     /* Vapor pressure deficit (Pa) */
  double Wind,                    /* Wind speed (m/s) */
  /* Snowpack Variables */
  double OldTSurf,                /* Surface temperature during previous time step */
  double IceDepth,                /* Depth of glacier surface layer (m) */
  double IceWE,                   /* Liquid water in the glacier surface layer (m) */
  /* Energy Balance Components */
  double Tair,                    /* Canopy air / Air temperature (C) */
  double TGrnd,                   /* Ground surface temperature (C) */
  double *AdvectedEnergy,         /* Energy advected by precipitation (W/m2) */
  double *AdvectedSensibleHeat,   /* Sensible heat advected from snow-free area into snow covered area (W/m^2) */
  double *DeltaColdContent,       /* Change in cold content of surface layer (W/m2) */
  double *GroundFlux,             /* Ground Heat Flux (W/m2) */
  double *LatentHeat,             /* Latent heat exchange at surface (W/m2) */
  double *LatentHeatSub,          /* Latent heat of sublimation exchange at surface (W/m2) */
  double *NetLongUnder,           /* Net longwave radiation at snowpack surface (W/m^2) */
  double *SensibleHeat,           /* Sensible heat exchange at surface (W/m2) */
  double *vapor_flux,             /* Mass flux of water vapor to or from the intercepted snow (m/timestep) */
  double *blowing_flux,           /* Mass flux of water vapor from blowing snow. (m/timestep) */
  double *surface_flux            /* Mass flux of water vapor from pack snow. (m/timestep) */
  ) :
      Dt(Dt), Ra(Ra), Ra_used(Ra_used), Displacement(Displacement), Z(Z),
      Z0(Z0), AirDens(AirDens), EactAir(EactAir), LongSnowIn(LongSnowIn),
      Lv(Lv), Press(Press), Rain(Rain), NetShortUnder(NetShortUnder), Vpd(Vpd),
      Wind(Wind), OldTSurf(OldTSurf), IceDepth(IceDepth), IceWE(IceWE),
      Tair(Tair), TGrnd(TGrnd), AdvectedEnergy(AdvectedEnergy),
      AdvectedSensibleHeat(AdvectedSensibleHeat), DeltaColdContent(DeltaColdContent),
      GroundFlux(GroundFlux), LatentHeat(LatentHeat), LatentHeatSub(LatentHeatSub),
      NetLongUnder(NetLongUnder), SensibleHeat(SensibleHeat), vapor_flux(vapor_flux),
      blowing_flux(blowing_flux), surface_flux(surface_flux) {}

  double calculate(double);
private:
  double Dt;                      /* Model time step (sec) */
  double Ra;                      /* Aerodynamic resistance (s/m) */
  double *Ra_used;                /* Aerodynamic resistance (s/m) after stability correction */
  /* Vegetation Parameters */
  double Displacement;            /* Displacement height (m) */
  double Z;                       /* Reference height (m) */
  double *Z0;                     /* surface roughness height (m) */
  /* Atmospheric Forcing Variables */
  double AirDens;                 /* Density of air (kg/m3) */
  double EactAir;                 /* Actual vapor pressure of air (Pa) */
  double LongSnowIn;              /* Incoming longwave radiation (W/m2) */
  double Lv;                      /* Latent heat of vaporization (J/kg3) */
  double Press;                   /* Air pressure (Pa) */
  double Rain;                    /* Rain fall (m/timestep) */
  double NetShortUnder;           /* Net incident shortwave radiation (W/m2) */
  double Vpd;                     /* Vapor pressure deficit (Pa) */
  double Wind;                    /* Wind speed (m/s) */
  /* Snowpack Variables */
  double OldTSurf;                /* Surface temperature during previous time step */
  double IceDepth;                /* Depth of glacier surface layer (m) */
  double IceWE;                   /* Liquid water in the glacier surface layer (m) */
  /* Energy Balance Components */
  double Tair;                    /* Canopy air / Air temperature (C) */
  double TGrnd;                   /* Ground surface temperature (C) */
  double *AdvectedEnergy;         /* Energy advected by precipitation (W/m2) */
  double *AdvectedSensibleHeat;   /* Sensible heat advected from snow-free area into snow covered area (W/m^2) */
  double *DeltaColdContent;       /* Change in cold content of surface layer (W/m2) */
  double *GroundFlux;             /* Ground Heat Flux (W/m2) */
  double *LatentHeat;             /* Latent heat exchange at surface (W/m2) */
  double *LatentHeatSub;          /* Latent heat of sublimation exchange at surface (W/m2) */
  double *NetLongUnder;           /* Net longwave radiation at snowpack surface (W/m^2) */
  double *SensibleHeat;           /* Sensible heat exchange at surface (W/m2) */
  double *vapor_flux;             /* Mass flux of water vapor to or from the intercepted snow (m/timestep) */
  double *blowing_flux;           /* Mass flux of water vapor from blowing snow. (m/timestep) */
  double *surface_flux;           /* Mass flux of water vapor from pack snow. (m/timestep) */
};


#endif /* GLACIERENERGYBALANCE_H_ */
