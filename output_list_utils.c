#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include <string.h>

static char vcid[] = "$Id$";

void copy_output_data(std::vector<OutputData*>&current_output_data, OutputData *out_data_list, const ProgramState *state) {
	OutputData *cur_data = new OutputData[N_OUTVAR_TYPES];

	for (int i = 0; i < N_OUTVAR_TYPES; i++) {
	cur_data[i].varname = out_data_list[i].varname;
	cur_data[i].format = out_data_list[i].format;
    cur_data[i].nelem = out_data_list[i].nelem;
    cur_data[i].aggtype = out_data_list[i].aggtype;
    cur_data[i].mult = out_data_list[i].mult;
    cur_data[i].type = out_data_list[i].type;
    cur_data[i].write = out_data_list[i].write;
    // Allocate space for data
    cur_data[i].data = new double[cur_data[i].nelem];
    cur_data[i].aggdata = new double[cur_data[i].nelem];
    for (int element = 0; element < cur_data[i].nelem; element++) {
      cur_data[i].data[element] = out_data_list[i].data[element];
      cur_data[i].aggdata[element] = out_data_list[i].aggdata[element];
    }
  }
  current_output_data.push_back(cur_data);
}

OutputData *create_output_list(const ProgramState* state) {
/*************************************************************
  create_output_list()      Ted Bohn     September 08, 2006

  This routine creates the list of output variables.

  Modifications:
  2006-Sep-14 Implemented ALMA-compliant input and output;
              now more variables are tracked.				TJB
  2006-Sep-18 Implemented aggregation of output variables.		TJB
  2006-Oct-10 Shortened the names of variables whose names were
	      too long; fixed typos in other names; added
	      OUT_IN_LONG.						TJB
  2006-Nov-07 Changed default precision from %.1f to %.4f.		TJB
  2006-Nov-07 Added OUT_SOIL_TNODE.					TJB
  2006-Nov-30 Added OUT_DELSURFSTOR.					TJB
  2007-Feb-28 Corrected AGG_TYPE definitions for miscellaneous
	      output variables; re-organized the code to make
	      it easier to debug.					TJB
  2007-Aug-17 Added EXCESS_ICE variables to output list.        	JCA
  2007-Aug-22 Added OUTPUT_WATER_ERROR as output variable.      	JCA
  2008-Sep-09 Added SOIL_TNODE_WL as an output variable, the
	      soil temperature in the wetland fraction of the
	      grid cell.						LCB via TJB
  2009-Jan-16 Added AERO_COND1&2 and AERO_RESIST1&2 to track
	      surface and overstory values; changed AERO_COND
	      and AERO_RESIST to track "scene" values.			TJB
  2009-Feb-22 Added OUT_VPD.						TJB
  2009-May-17 Added OUT_ASAT.						TJB
  2009-Jun-09 Added OUT_PET_*, potential evap computed for
	      various landcover types.					TJB
  2009-Jun-19 Added T flag to indicate whether TFALLBACK occurred.	TJB
  2009-Jul-07 Fixed nelem assignments for some band-specific vars.	TJB
  2009-Sep-19 Changed "*_FLAG" to "*_FBFLAG".				TJB
  2009-Oct-08 Extended T fallback scheme to snow and ice T.		TJB
  2010-Feb-14 Added OUT_LAKE_AREA_FRAC.					TJB
  2010-Mar-31 Added OUT_RUNOFF_IN.					TJB
  2010-Sep-24 Renamed RUNOFF_IN and OUT_RUNOFF_IN to CHANNEL_IN and
	      OUT_LAKE_CHAN_IN, respectively.  Renamed OUT_EVAP_LAKE
	      to OUT_LAKE_EVAP.  Added other lake water balance terms
	      to set of output variables.  Added volumetric versions 
	      of these too.						TJB
  2010-Nov-02 Added OUT_LAKE_RO_IN and OUT_LAKE_RO_IN_V for reporting
	      overland runoff input to lake.  Added OUT_LAKE_RCHRG and
	      OUT_LAKE_RCHRG_V for reporting lake recharge of
	      surrounding wetland.  Added OUT_LAKE_VAPFLX and
	      OUT_LAKE_VAPFLX_V.					TJB
  2010-Nov-21 Added OUT_LAKE_DSTOR, OUT_LAKE_DSTOR_V, OUT_LAKE_DSWE,
	      OUT_LAKE_DSWE_V, OUT_LAKE_SWE, and OUT_LAKE_SWE_V.	TJB
  2010-Dec-01 Added OUT_ZWT.						TJB
  2011-Mar-01 Added OUT_ZWT2, OUT_ZWT3, and OUT_ZWTL.			TJB
  2011-Nov-04 Added OUT_TSKC.						TJB
*************************************************************/

  int v;
  OutputData *out_data;

  out_data = new OutputData[N_OUTVAR_TYPES];

  // Build the list of supported output variables

  // Water Balance Terms - state variables
  out_data[OUT_ASAT].varname = "OUT_ASAT";                       /* saturated area fraction */
  out_data[OUT_LAKE_AREA_FRAC].varname = "OUT_LAKE_AREA_FRAC";   /* lake surface area as fraction of grid cell area [fraction] */
  out_data[OUT_LAKE_DEPTH].varname = "OUT_LAKE_DEPTH";           /* lake depth [m] */
  out_data[OUT_LAKE_ICE].varname = "OUT_LAKE_ICE";               /* moisture stored as lake ice [mm] */
  out_data[OUT_LAKE_ICE_FRACT].varname = "OUT_LAKE_ICE_FRACT";   /* fractional coverage of lake ice [fraction] */
  out_data[OUT_LAKE_ICE_HEIGHT].varname = "OUT_LAKE_ICE_HEIGHT"; /* thickness of lake ice [cm] */
  out_data[OUT_LAKE_MOIST].varname = "OUT_LAKE_MOIST";           /* liquid water stored in lake [mm over lake area?] */
  out_data[OUT_LAKE_SURF_AREA].varname = "OUT_LAKE_SURF_AREA";   /* lake surface area [m2] */
  out_data[OUT_LAKE_SWE].varname = "OUT_LAKE_SWE";               /* liquid water equivalent of snow on top of lake ice [m over lake ice] */
  out_data[OUT_LAKE_SWE_V].varname = "OUT_LAKE_SWE_V";           /* volumetric liquid water equivalent of snow on top of lake ice [m3] */
  out_data[OUT_LAKE_VOLUME].varname = "OUT_LAKE_VOLUME";         /* lake volume [m3] */
  out_data[OUT_ROOTMOIST].varname = "OUT_ROOTMOIST";             /* root zone soil moisture [mm] */
  out_data[OUT_SMFROZFRAC].varname = "OUT_SMFROZFRAC";           /* fraction of soil moisture (by mass) that is ice, for each soil layer */
  out_data[OUT_SMLIQFRAC].varname = "OUT_SMLIQFRAC";             /* fraction of soil moisture (by mass) that is liquid, for each soil layer */
  out_data[OUT_SNOW_CANOPY].varname = "OUT_SNOW_CANOPY";         /* snow interception storage in canopy [mm] */
  out_data[OUT_SNOW_COVER].varname = "OUT_SNOW_COVER";           /* fractional area of snow cover [fraction] */
  out_data[OUT_SNOW_DEPTH].varname = "OUT_SNOW_DEPTH";           /* depth of snow pack [cm] */
  out_data[OUT_SOIL_ICE].varname = "OUT_SOIL_ICE";               /* soil ice content [mm] for each soil layer */
  out_data[OUT_SOIL_ICE_TOT].varname = "OUT_SOIL_ICE_TOT";       /* soil ice content [mm] for all soil layers */
  out_data[OUT_SOIL_LIQ].varname = "OUT_SOIL_LIQ";               /* soil liquid moisture content [mm] for each soil layer */
  out_data[OUT_SOIL_LIQ_TOT].varname = "OUT_SOIL_LIQ_TOT";       /* soil liquid moisture content [mm] for all soil layers */
  out_data[OUT_SOIL_MOIST].varname = "OUT_SOIL_MOIST";           /* soil total moisture content [mm] for each soil layer */
  out_data[OUT_SOIL_MOIST_TOT].varname = "OUT_SOIL_MOIST_TOT";   /* soil total moisture content [mm] for all soil layers */
  out_data[OUT_SOIL_WET].varname = "OUT_SOIL_WET";               /* vertical average of (soil moisture - wilting point)/(maximum soil moisture - wilting point) [mm/mm] */
  out_data[OUT_SURFSTOR].varname = "OUT_SURFSTOR";               /* storage of liquid water on surface (ponding) [mm] */
  out_data[OUT_SURF_FROST_FRAC].varname = "OUT_SURF_FROST_FRAC"; /* fraction of soil surface that is frozen [fraction] */
  out_data[OUT_SWE].varname = "OUT_SWE";                         /* snow water equivalent in snow pack [mm] */
  out_data[OUT_WDEW].varname = "OUT_WDEW";                       /* total moisture interception storage in canopy [mm] */
  out_data[OUT_ZWT].varname = "OUT_ZWT";                         /* water table position [cm] - method 1 (zwt within lowest unsaturated layer) */
  out_data[OUT_ZWT2].varname = "OUT_ZWT2";                       /* water table position [cm] - method 2 (zwt of total moisture across top-most N-1 layers, lumped together) */
  out_data[OUT_ZWT3].varname = "OUT_ZWT3";                       /* water table position [cm] - method 3 (zwt of total moisture across all layers, lumped together) */
  out_data[OUT_ZWTL].varname = "OUT_ZWTL";                       /* per-layer water table positions [cm] (one per soil layer) */

  // Water Balance Terms - fluxes
  out_data[OUT_BASEFLOW].varname = "OUT_BASEFLOW";               /* baseflow out of the bottom layer [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_DELINTERCEPT].varname = "OUT_DELINTERCEPT";       /* change in canopy interception storage [mm] */
  out_data[OUT_DELSOILMOIST].varname = "OUT_DELSOILMOIST";       /* change in soil water content [mm] */
  out_data[OUT_DELSWE].varname = "OUT_DELSWE";                   /* change in snow water equivalent [mm] */
  out_data[OUT_DELSURFSTOR].varname = "OUT_DELSURFSTOR";         /* change in surface liquid water storage  [mm] */
  out_data[OUT_EVAP].varname = "OUT_EVAP";                       /* total net evaporation [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_EVAP_BARE].varname = "OUT_EVAP_BARE";             /* net evaporation from bare soil [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_EVAP_CANOP].varname = "OUT_EVAP_CANOP";           /* net evaporation from canopy interception [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_INFLOW].varname = "OUT_INFLOW";                   /* moisture that reaches top of soil column [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_LAKE_BF_IN].varname = "OUT_LAKE_BF_IN";           /* incoming baseflow from lake catchment [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_LAKE_BF_IN_V].varname = "OUT_LAKE_BF_IN_V";       /* incoming volumetric baseflow from lake catchment [m3] (ALMA_OUTPUT: [m3/s]) */
  out_data[OUT_LAKE_BF_OUT].varname = "OUT_LAKE_BF_OUT";         /* outgoing baseflow lake [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_LAKE_BF_OUT_V].varname = "OUT_LAKE_BF_OUT_V";     /* outgoing volumetric baseflow from lake [m3] (ALMA_OUTPUT: [m3/s]) */
  out_data[OUT_LAKE_CHAN_IN].varname = "OUT_LAKE_CHAN_IN";       /* channel inflow into lake [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_LAKE_CHAN_IN_V].varname = "OUT_LAKE_CHAN_IN_V";   /* volumetric channel inflow into lake [m3] (ALMA_OUTPUT: [m3/s]) */
  out_data[OUT_LAKE_CHAN_OUT].varname = "OUT_LAKE_CHAN_OUT";     /* channel outflow from lake [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_LAKE_CHAN_OUT_V].varname = "OUT_LAKE_CHAN_OUT_V"; /* volumetric channel outflow from lake [m3] (ALMA_OUTPUT: [m3/s]) */
  out_data[OUT_LAKE_DSTOR].varname = "OUT_LAKE_DSTOR";           /* change in lake moisture storage (liquid plus ice cover) [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_LAKE_DSTOR_V].varname = "OUT_LAKE_DSTOR_V";       /* volumetric change in lake moisture storage (liquid plus ice cover) [m3] (ALMA_OUTPUT: [m3/s]) */
  out_data[OUT_LAKE_DSWE].varname = "OUT_LAKE_DSWE";             /* change in snowpack on top of lake ice [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_LAKE_DSWE_V].varname = "OUT_LAKE_DSWE_V";         /* volumetric change in snowpack on top of lake ice [m3] (ALMA_OUTPUT: [m3/s]) */
  out_data[OUT_LAKE_EVAP].varname = "OUT_LAKE_EVAP";             /* net evaporation from lake surface [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_LAKE_EVAP_V].varname = "OUT_LAKE_EVAP_V";         /* net volumetric evaporation from lake surface [m3] (ALMA_OUTPUT: [m3/s]) */
  out_data[OUT_LAKE_PREC_V].varname = "OUT_LAKE_PREC_V";         /* volumetric precipitation over lake surface [m3] (ALMA_OUTPUT: [m3/s]) */
  out_data[OUT_LAKE_RCHRG].varname = "OUT_LAKE_RCHRG";           /* recharge from lake to surrounding wetland [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_LAKE_RCHRG_V].varname = "OUT_LAKE_RCHRG_V";       /* volumetric recharge from lake to surrounding wetland [m3] (ALMA_OUTPUT: [m3/s]) */
  out_data[OUT_LAKE_RO_IN].varname = "OUT_LAKE_RO_IN";           /* incoming runoff from lake catchment [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_LAKE_RO_IN_V].varname = "OUT_LAKE_RO_IN_V";       /* incoming volumetric runoff from lake catchment [m3] (ALMA_OUTPUT: [m3/s]) */
  out_data[OUT_LAKE_VAPFLX].varname = "OUT_LAKE_VAPFLX";         /* sublimation from lake snow pack [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_LAKE_VAPFLX_V].varname = "OUT_LAKE_VAPFLX_V";     /* volumetric sublimation from lake snow pack [m3] (ALMA_OUTPUT: [m3/s]) */
  out_data[OUT_PET_SATSOIL].varname = "OUT_PET_SATSOIL";         /* potential evap from saturated bare soil [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_PET_H2OSURF].varname = "OUT_PET_H2OSURF";         /* potential evap from open water [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_PET_SHORT].varname = "OUT_PET_SHORT";             /* potential evap from short reference crop (grass) [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_PET_TALL].varname = "OUT_PET_TALL";               /* potential evap from tall reference crop (alfalfa) [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_PET_NATVEG].varname = "OUT_PET_NATVEG";           /* potential evap from current vegetation and current canopy resistance [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_PET_VEGNOCR].varname = "OUT_PET_VEGNOCR";         /* potential evap from current vegetation and 0 canopy resistance bare soil [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_PREC].varname = "OUT_PREC";                       /* incoming precipitation [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_RAINF].varname = "OUT_RAINF";                     /* rainfall [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_REFREEZE].varname = "OUT_REFREEZE";               /* refreezing of water in the snow [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_RUNOFF].varname = "OUT_RUNOFF";                   /* surface runoff [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_SNOW_MELT].varname = "OUT_SNOW_MELT";             /* snow melt [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_SNOWF].varname = "OUT_SNOWF";                     /* snowfall [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_SUB_BLOWING].varname = "OUT_SUB_BLOWING";         /* net sublimation of blowing snow [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_SUB_CANOP].varname = "OUT_SUB_CANOP";             /* net sublimation from snow stored in canopy [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_SUB_SNOW].varname = "OUT_SUB_SNOW";               /* net sublimation from snow pack (surface and blowing) [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_SUB_SURFACE].varname = "OUT_SUB_SURFACE";         /* net sublimation from snow pack surface [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_TRANSP_VEG].varname = "OUT_TRANSP_VEG";           /* net transpiration from vegetation [mm] (ALMA_OUTPUT: [mm/s]) */

  // Energy Balance Terms - state variables
  out_data[OUT_ALBEDO].varname = "OUT_ALBEDO";                   /* albedo [fraction] */
  out_data[OUT_BARESOILT].varname = "OUT_BARESOILT";             /* bare soil surface temperature [C] (ALMA_OUTPUT: [K]) */
  out_data[OUT_FDEPTH].varname = "OUT_FDEPTH";                   /* depth of freezing fronts [cm] (ALMA_OUTPUT: [m]) for each freezing front */
  out_data[OUT_LAKE_ICE_TEMP].varname = "OUT_LAKE_ICE_TEMP";     /* lake ice temperature [K] */
  out_data[OUT_LAKE_SURF_TEMP].varname = "OUT_LAKE_SURF_TEMP";   /* lake surface temperature [K] */
  out_data[OUT_RAD_TEMP].varname = "OUT_RAD_TEMP";               /* average radiative surface temperature [K] */
  out_data[OUT_SALBEDO].varname = "OUT_SALBEDO";                 /* snow albedo [fraction] */
  out_data[OUT_SNOW_PACK_TEMP].varname = "OUT_SNOW_PACK_TEMP";   /* snow pack temperature [C] (ALMA_OUTPUT: [K]) */
  out_data[OUT_SNOW_SURF_TEMP].varname = "OUT_SNOW_SURF_TEMP";   /* snow surface temperature [C] (ALMA_OUTPUT: [K]) */
  out_data[OUT_SNOWT_FBFLAG].varname = "OUT_SNOWT_FBFLAG";       /* snow surface temperature flag */
  out_data[OUT_SOIL_TEMP].varname = "OUT_SOIL_TEMP";             /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil layer */
  out_data[OUT_SOIL_TNODE].varname = "OUT_SOIL_TNODE";           /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil thermal node */
  out_data[OUT_SOIL_TNODE_WL].varname = "OUT_SOIL_TNODE_WL";     /* soil temperature [C] (ALMA_OUTPUT: [K]) for each soil thermal node in the wetland */
  out_data[OUT_SOILT_FBFLAG].varname = "OUT_SOILT_FBFLAG";       /* soil temperature flag for each soil thermal node */
  out_data[OUT_SURF_TEMP].varname = "OUT_SURF_TEMP";             /* average surface temperature [C] (ALMA_OUTPUT: [K]) */
  out_data[OUT_SURFT_FBFLAG].varname = "OUT_SURFT_FBFLAG";       /* surface temperature flag */
  out_data[OUT_TCAN_FBFLAG].varname = "OUT_TCAN_FBFLAG";         /* Tcanopy flag */
  out_data[OUT_TDEPTH].varname = "OUT_TDEPTH";                   /* depth of thawing fronts [cm] (ALMA_OUTPUT: [m]) for each thawing front */
  out_data[OUT_TFOL_FBFLAG].varname = "OUT_TFOL_FBFLAG";         /* Tfoliage flag */
  out_data[OUT_VEGT].varname = "OUT_VEGT";                       /* average vegetation canopy temperature [C] (ALMA_OUTPUT: [K]) */

  // Energy Balance Terms - fluxes
  out_data[OUT_ADV_SENS].varname = "OUT_ADV_SENS";               /* net sensible heat advected to snow pack [W/m2] */
  out_data[OUT_ADVECTION].varname = "OUT_ADVECTION";             /* advected energy [W/m2] */
  out_data[OUT_DELTACC].varname = "OUT_DELTACC";                 /* rate of change in cold content in snow pack [W/m2] */
  out_data[OUT_DELTAH].varname = "OUT_DELTAH";                   /* rate of change in heat storage [W/m2] */
  out_data[OUT_ENERGY_ERROR].varname = "OUT_ENERGY_ERROR";       /* energy budget error [W/m2] */
  out_data[OUT_WATER_ERROR].varname = "OUT_WATER_ERROR";         /* water budget error [mm] */
  out_data[OUT_FUSION].varname = "OUT_FUSION";                   /* net energy used to melt/freeze soil moisture [W/m2] */
  out_data[OUT_GRND_FLUX].varname = "OUT_GRND_FLUX";             /* net heat flux into ground [W/m2] */
  out_data[OUT_IN_LONG].varname = "OUT_IN_LONG";                 /* incoming longwave flux at surface (under veg) [W/m2] */
  out_data[OUT_LATENT].varname = "OUT_LATENT";                   /* net upward latent heat flux [W/m2] */
  out_data[OUT_LATENT_SUB].varname = "OUT_LATENT_SUB";           /* net upward latent heat flux from sublimation [W/m2] */
  out_data[OUT_MELT_ENERGY].varname = "OUT_MELT_ENERGY";         /* energy of fusion (melting) [W/m2] */
  out_data[OUT_NET_LONG].varname = "OUT_NET_LONG";               /* net downward longwave flux [W/m2] */
  out_data[OUT_NET_SHORT].varname = "OUT_NET_SHORT";             /* net downward shortwave flux [W/m2] */
  out_data[OUT_R_NET].varname = "OUT_R_NET";                     /* net downward radiation flux [W/m2] */
  out_data[OUT_RFRZ_ENERGY].varname = "OUT_RFRZ_ENERGY";         /* net energy used to refreeze liquid water in snowpack [W/m2] */
  out_data[OUT_SENSIBLE].varname = "OUT_SENSIBLE";               /* net upward sensible heat flux [W/m2] */
  out_data[OUT_SNOW_FLUX].varname = "OUT_SNOW_FLUX";             /* energy flux through snow pack [W/m2] */

  // Miscellaneous Terms
  out_data[OUT_AERO_COND].varname = "OUT_AERO_COND";             /* "scene" aerodynamic conductance [m/s] (tiles with overstory contribute overstory conductance; others contribue surface conductance) */
  out_data[OUT_AERO_COND1].varname = "OUT_AERO_COND1";           /* surface aerodynamic conductance [m/s] */
  out_data[OUT_AERO_COND2].varname = "OUT_AERO_COND2";           /* overstory aerodynamic conductance [m/s] */
  out_data[OUT_AERO_RESIST].varname = "OUT_AERO_RESIST";         /* "scene" aerodynamic resistance [s/m] (tiles with overstory contribute overstory resistance; others contribue surface resistance)*/
  out_data[OUT_AERO_RESIST1].varname = "OUT_AERO_RESIST1";       /* surface aerodynamic resistance [m/s] */
  out_data[OUT_AERO_RESIST2].varname = "OUT_AERO_RESIST2";       /* overstory aerodynamic resistance [m/s] */
  out_data[OUT_AIR_TEMP].varname = "OUT_AIR_TEMP";               /* air temperature [C] */
  out_data[OUT_DENSITY].varname = "OUT_DENSITY";                 /* near-surface atmospheric density [kg/m3] */
  out_data[OUT_LONGWAVE].varname = "OUT_LONGWAVE";               /* incoming longwave [W/m2] */
  out_data[OUT_PRESSURE].varname = "OUT_PRESSURE";               /* near surface atmospheric pressure [kPa] */
  out_data[OUT_QAIR].varname = "OUT_QAIR";                       /* specific humidity [kg/kg] */
  out_data[OUT_REL_HUMID].varname = "OUT_REL_HUMID";             /* relative humidity [fraction]*/
  out_data[OUT_SHORTWAVE].varname = "OUT_SHORTWAVE";             /* incoming shortwave [W/m2] */
  out_data[OUT_SURF_COND].varname = "OUT_SURF_COND";             /* surface conductance [m/s] */
  out_data[OUT_TSKC].varname = "OUT_TSKC";                       /* cloud cover fraction [fraction] */
  out_data[OUT_VP].varname = "OUT_VP";                           /* near surface vapor pressure [kPa] */
  out_data[OUT_VPD].varname = "OUT_VPD";                         /* near surface vapor pressure deficit [kPa] */
  out_data[OUT_WIND].varname = "OUT_WIND";                       /* near surface wind speed [m/s] */
  out_data[OUT_AREA].varname = "OUT_AREA";                       /* cell area [m2] */

  // Dynamic Soil Layer Terms - EXCESS_ICE option
#if EXCESS_ICE
  out_data[OUT_SOIL_DEPTH].varname = "OUT_SOIL_DEPTH";             /* soil moisture layer depths [m] */
  out_data[OUT_SUBSIDENCE].varname = "OUT_SUBSIDENCE";             /* subsidence of soil layer [mm] */
  out_data[OUT_POROSITY].varname = "OUT_POROSITY";                 /* porosity [mm/mm] */
  out_data[OUT_ZSUM_NODE].varname = "OUT_ZSUM_NODE";               /* depths of thermal nodes [m] */
#endif

  // Band-specific quantities
  out_data[OUT_ADV_SENS_BAND].varname = "OUT_ADV_SENS_BAND";               /* net sensible heat flux advected to snow pack [W/m2] */
  out_data[OUT_ADVECTION_BAND].varname = "OUT_ADVECTION_BAND";             /* advected energy [W/m2] */
  out_data[OUT_ALBEDO_BAND].varname = "OUT_ALBEDO_BAND";                   /* albedo [fraction] */
  out_data[OUT_AREA_BAND].varname = "OUT_AREA_BAND";                       /* band area [fraction] */
  out_data[OUT_DELTACC_BAND].varname = "OUT_DELTACC_BAND";                 /* change in cold content in snow pack [W/m2] */
  out_data[OUT_ELEV_BAND].varname = "OUT_ELEV_BAND";                       /* band median elevation [m] */
  out_data[OUT_GRND_FLUX_BAND].varname = "OUT_GRND_FLUX_BAND";             /* net heat flux into ground [W/m2] */
  out_data[OUT_IN_LONG_BAND].varname = "OUT_IN_LONG_BAND";                 /* incoming longwave flux at surface (under veg) [W/m2] */
  out_data[OUT_LATENT_BAND].varname = "OUT_LATENT_BAND";                   /* net upward latent heat flux [W/m2] */
  out_data[OUT_LATENT_SUB_BAND].varname = "OUT_LATENT_SUB_BAND";           /* net upward latent heat flux from sublimation [W/m2] */
  out_data[OUT_MELT_ENERGY_BAND].varname = "OUT_MELT_ENERGY_BAND";         /* energy of fusion (melting) [W/m2] */
  out_data[OUT_NET_LONG_BAND].varname = "OUT_NET_LONG_BAND";               /* net downward longwave flux [W/m2] */
  out_data[OUT_NET_SHORT_BAND].varname = "OUT_NET_SHORT_BAND";             /* net downward shortwave flux [W/m2] */
  out_data[OUT_RFRZ_ENERGY_BAND].varname = "OUT_RFRZ_ENERGY_BAND";         /* net energy used to refreeze liquid water in snowpack [W/m2] */
  out_data[OUT_SENSIBLE_BAND].varname = "OUT_SENSIBLE_BAND";               /* net upward sensible heat flux [W/m2] */
  out_data[OUT_SNOW_CANOPY_BAND].varname = "OUT_SNOW_CANOPY_BAND";         /* snow interception storage in canopy [mm] */
  out_data[OUT_SNOW_COVER_BAND].varname = "OUT_SNOW_COVER_BAND";           /* fractional area of snow cover [fraction] */
  out_data[OUT_SNOW_DEPTH_BAND].varname = "OUT_SNOW_DEPTH_BAND";           /* depth of snow pack [cm] */
  out_data[OUT_SNOW_FLUX_BAND].varname = "OUT_SNOW_FLUX_BAND";             /* energy flux through snow pack [W/m2] */
  out_data[OUT_SNOW_MELT_BAND].varname = "OUT_SNOW_MELT_BAND";             /* snow melt [mm] (ALMA_OUTPUT: [mm/s]) */
  out_data[OUT_SNOW_PACKT_BAND].varname = "OUT_SNOW_PACKT_BAND";           /* snow pack temperature [C] (ALMA_OUTPUT: [K]) */
  out_data[OUT_SNOW_SURFT_BAND].varname = "OUT_SNOW_SURFT_BAND";           /* snow surface temperature [C] (ALMA_OUTPUT: [K]) */
  out_data[OUT_SWE_BAND].varname = "OUT_SWE_BAND";                         /* snow water equivalent in snow pack [mm] */

  //Glacier Terms
  out_data[OUT_GLAC_WAT_STOR].varname =  "OUT_GLAC_WAT_STOR";              /* glacier water storage [mm] */
  out_data[OUT_GLAC_AREA].varname =  "OUT_GLAC_AREA";                      /* glacier surface area fraction */
  out_data[OUT_GLAC_MBAL].varname =  "OUT_GLAC_MBAL";                      /* glacier mass balance [mm] */
  out_data[OUT_GLAC_IMBAL].varname =  "OUT_GLAC_IMBAL";                    /* glacier ice mass balance [mm] */
  out_data[OUT_GLAC_ACCUM].varname =  "OUT_GLAC_ACCUM";                    /* glacier ice accumulation from conversion of firn to ice [mm] */
  out_data[OUT_GLAC_MELT].varname =  "OUT_GLAC_MELT";                      /* glacier ice melt [mm] */
  out_data[OUT_GLAC_SUB].varname =  "OUT_GLAC_SUB";                        /* Net sublimation of glacier ice [mm] */
  out_data[OUT_GLAC_INFLOW].varname =  "OUT_GLAC_INFLOW";                  /* glacier water inflow from snow melt, ice melt and rainfall [mm] */
  out_data[OUT_GLAC_OUTFLOW].varname =  "OUT_GLAC_OUTFLOW";                /* glacier water outflow [mm] */
  out_data[OUT_GLAC_SURF_TEMP].varname =  "OUT_GLAC_SURF_TEMP";            /* glacier surface temperature [C] */
  out_data[OUT_GLAC_TSURF_FBFLAG].varname =  "OUT_GLAC_TSURF_FBFLAG";      /* glacier surface temperature flag */
  out_data[OUT_GLAC_DELTACC].varname =  "OUT_GLAC_DELTACC";                /* rate of change of cold content in glacier surface layer [W/m2] */
  out_data[OUT_GLAC_FLUX].varname =  "OUT_GLAC_FLUX";                      /* energy flux through glacier surface layer [W/m2] */
  out_data[OUT_GLAC_OUTFLOW_COEF].varname =  "OUT_GLAC_OUTFLOW_COEF";      /* glacier outflow coefficient [fraction] */
  out_data[OUT_GLAC_DELTACC_BAND].varname =  "OUT_GLAC_DELTACC_BAND";      /* rate of change of cold content in glacier surface layer [W/m2] */
  out_data[OUT_GLAC_FLUX_BAND].varname =  "OUT_GLAC_FLUX_BAND";            /* energy flux through glacier surface layer [W/m2] */
  out_data[OUT_GLAC_WAT_STOR_BAND].varname =  "OUT_GLAC_WAT_STOR_BAND";    /* glacier water storage [mm] */
  out_data[OUT_GLAC_AREA_BAND].varname =  "OUT_GLAC_AREA_BAND";            /* glacier surface area fraction */
  out_data[OUT_GLAC_MBAL_BAND].varname =  "OUT_GLAC_MBAL_BAND";            /* glacier mass balance [mm] */
  out_data[OUT_GLAC_IMBAL_BAND].varname =  "OUT_GLAC_IMBAL_BAND";          /* glacier ice mass balance [mm] */
  out_data[OUT_GLAC_ACCUM_BAND].varname =  "OUT_GLAC_ACCUM_BAND";          /* glacier ice accumulation from conversion of firn to ice [mm] */
  out_data[OUT_GLAC_MELT_BAND].varname =  "OUT_GLAC_MELT_BAND";            /* glacier ice melt [mm] */
  out_data[OUT_GLAC_SUB_BAND].varname =  "OUT_GLAC_SUB_BAND";              /* Net sublimation of glacier ice [mm] */
  out_data[OUT_GLAC_INFLOW_BAND].varname =  "OUT_GLAC_INFLOW_BAND";        /* glacier water inflow from snow melt, ice melt and rainfall [mm] */
  out_data[OUT_GLAC_OUTFLOW_BAND].varname =  "OUT_GLAC_OUTFLOW_BAND";      /* glacier water outflow [mm] */

  // Set number of elements - default is 1
  for (v=0; v<N_OUTVAR_TYPES; v++) {
    out_data[v].nelem = 1;
  }
  if (state->options.FROZEN_SOIL) {
    out_data[OUT_FDEPTH].nelem = MAX_FRONTS;
    out_data[OUT_TDEPTH].nelem = MAX_FRONTS;
  }
  out_data[OUT_SMLIQFRAC].nelem = state->options.Nlayer;
  out_data[OUT_SMFROZFRAC].nelem = state->options.Nlayer;
  out_data[OUT_SOIL_ICE].nelem = state->options.Nlayer;
  out_data[OUT_SOIL_LIQ].nelem = state->options.Nlayer;
  out_data[OUT_SOIL_MOIST].nelem = state->options.Nlayer;
  out_data[OUT_SOIL_TEMP].nelem = state->options.Nlayer;
  out_data[OUT_ZWTL].nelem = state->options.Nlayer;
#if EXCESS_ICE
  out_data[OUT_SOIL_DEPTH].nelem = state->options.Nlayer;
  out_data[OUT_SUBSIDENCE].nelem = state->options.Nlayer;
  out_data[OUT_POROSITY].nelem = state->options.Nlayer;
  out_data[OUT_ZSUM_NODE].nelem = state->options.Nnode;
#endif
  out_data[OUT_SOIL_TNODE].nelem = state->options.Nnode;
  out_data[OUT_SOIL_TNODE_WL].nelem = state->options.Nnode;
  out_data[OUT_SOILT_FBFLAG].nelem = state->options.Nnode;
  out_data[OUT_ADV_SENS_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_ADVECTION_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_ALBEDO_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_AREA_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_DELTACC_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_ELEV_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_GRND_FLUX_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_IN_LONG_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_LATENT_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_LATENT_SUB_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_MELT_ENERGY_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_NET_LONG_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_NET_SHORT_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_RFRZ_ENERGY_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_SENSIBLE_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_SNOW_CANOPY_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_SNOW_COVER_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_SNOW_DEPTH_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_SNOW_FLUX_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_SNOW_MELT_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_SNOW_PACKT_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_SNOW_SURFT_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_SWE_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_GLAC_DELTACC_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_GLAC_FLUX_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_GLAC_WAT_STOR_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_GLAC_AREA_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_GLAC_MBAL_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_GLAC_IMBAL_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_GLAC_ACCUM_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_GLAC_MELT_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_GLAC_SUB_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_GLAC_INFLOW_BAND].nelem = state->options.SNOW_BAND;
  out_data[OUT_GLAC_OUTFLOW_BAND].nelem = state->options.SNOW_BAND;

  // Set aggregation method - default is to average over the interval
  for (v=0; v<N_OUTVAR_TYPES; v++) {
    out_data[v].aggtype = AGG_TYPE_AVG;
  }
  out_data[OUT_ASAT].aggtype = AGG_TYPE_END;
  out_data[OUT_LAKE_AREA_FRAC].aggtype = AGG_TYPE_END;
  out_data[OUT_LAKE_DEPTH].aggtype = AGG_TYPE_END;
  out_data[OUT_LAKE_ICE].aggtype = AGG_TYPE_END;
  out_data[OUT_LAKE_ICE_FRACT].aggtype = AGG_TYPE_END;
  out_data[OUT_LAKE_ICE_HEIGHT].aggtype = AGG_TYPE_END;
  out_data[OUT_LAKE_MOIST].aggtype = AGG_TYPE_END;
  out_data[OUT_LAKE_SURF_AREA].aggtype = AGG_TYPE_END;
  out_data[OUT_LAKE_SWE].aggtype = AGG_TYPE_END;
  out_data[OUT_LAKE_SWE_V].aggtype = AGG_TYPE_END;
  out_data[OUT_LAKE_VOLUME].aggtype = AGG_TYPE_END;
  out_data[OUT_ROOTMOIST].aggtype = AGG_TYPE_END;
  out_data[OUT_SMFROZFRAC].aggtype = AGG_TYPE_END;
  out_data[OUT_SMLIQFRAC].aggtype = AGG_TYPE_END;
  out_data[OUT_SNOW_CANOPY].aggtype = AGG_TYPE_END;
  out_data[OUT_SNOW_COVER].aggtype = AGG_TYPE_END;
  out_data[OUT_SNOW_DEPTH].aggtype = AGG_TYPE_END;
  out_data[OUT_SOIL_ICE].aggtype = AGG_TYPE_END;
  out_data[OUT_SOIL_ICE_TOT].aggtype = AGG_TYPE_END;
  out_data[OUT_SOIL_LIQ].aggtype = AGG_TYPE_END;
  out_data[OUT_SOIL_LIQ_TOT].aggtype = AGG_TYPE_END;
  out_data[OUT_SOIL_MOIST].aggtype = AGG_TYPE_END;
  out_data[OUT_SOIL_MOIST_TOT].aggtype = AGG_TYPE_END;
  out_data[OUT_SOIL_WET].aggtype = AGG_TYPE_END;
  out_data[OUT_SURFSTOR].aggtype = AGG_TYPE_END;
  out_data[OUT_SURF_FROST_FRAC].aggtype = AGG_TYPE_END;
  out_data[OUT_SWE].aggtype = AGG_TYPE_END;
  out_data[OUT_WDEW].aggtype = AGG_TYPE_END;
  out_data[OUT_ZWT].aggtype = AGG_TYPE_END;
  out_data[OUT_ZWT2].aggtype = AGG_TYPE_END;
  out_data[OUT_ZWT3].aggtype = AGG_TYPE_END;
  out_data[OUT_ZWTL].aggtype = AGG_TYPE_END;
  out_data[OUT_SNOW_CANOPY_BAND].aggtype = AGG_TYPE_END;
  out_data[OUT_SNOW_COVER_BAND].aggtype = AGG_TYPE_END;
  out_data[OUT_SNOW_DEPTH_BAND].aggtype = AGG_TYPE_END;
  out_data[OUT_SWE_BAND].aggtype = AGG_TYPE_END;
#if EXCESS_ICE
  out_data[OUT_SOIL_DEPTH].aggtype = AGG_TYPE_END;
  out_data[OUT_POROSITY].aggtype = AGG_TYPE_END;
  out_data[OUT_ZSUM_NODE].aggtype = AGG_TYPE_END;
  out_data[OUT_SUBSIDENCE].aggtype = AGG_TYPE_SUM;
#endif
  out_data[OUT_BASEFLOW].aggtype = AGG_TYPE_SUM;
  out_data[OUT_DELINTERCEPT].aggtype = AGG_TYPE_SUM;
  out_data[OUT_DELSOILMOIST].aggtype = AGG_TYPE_SUM;
  out_data[OUT_DELSWE].aggtype = AGG_TYPE_SUM;
  out_data[OUT_DELSURFSTOR].aggtype = AGG_TYPE_SUM;
  out_data[OUT_EVAP].aggtype = AGG_TYPE_SUM;
  out_data[OUT_EVAP_BARE].aggtype = AGG_TYPE_SUM;
  out_data[OUT_EVAP_CANOP].aggtype = AGG_TYPE_SUM;
  out_data[OUT_INFLOW].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_BF_IN].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_BF_IN_V].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_BF_OUT].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_BF_OUT_V].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_CHAN_IN].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_CHAN_IN_V].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_CHAN_OUT].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_CHAN_OUT_V].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_DSTOR].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_DSTOR_V].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_DSWE].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_DSWE_V].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_EVAP].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_EVAP_V].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_PREC_V].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_RCHRG].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_RCHRG_V].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_RO_IN].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_RO_IN_V].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_VAPFLX].aggtype = AGG_TYPE_SUM;
  out_data[OUT_LAKE_VAPFLX_V].aggtype = AGG_TYPE_SUM;
  out_data[OUT_PET_SATSOIL].aggtype = AGG_TYPE_SUM;
  out_data[OUT_PET_H2OSURF].aggtype = AGG_TYPE_SUM;
  out_data[OUT_PET_SHORT].aggtype = AGG_TYPE_SUM;
  out_data[OUT_PET_TALL].aggtype = AGG_TYPE_SUM;
  out_data[OUT_PET_NATVEG].aggtype = AGG_TYPE_SUM;
  out_data[OUT_PET_VEGNOCR].aggtype = AGG_TYPE_SUM;
  out_data[OUT_PREC].aggtype = AGG_TYPE_SUM;
  out_data[OUT_RAINF].aggtype = AGG_TYPE_SUM;
  out_data[OUT_REFREEZE].aggtype = AGG_TYPE_SUM;
  out_data[OUT_RUNOFF].aggtype = AGG_TYPE_SUM;
  out_data[OUT_SNOW_MELT].aggtype = AGG_TYPE_SUM;
  out_data[OUT_SNOWF].aggtype = AGG_TYPE_SUM;
  out_data[OUT_SUB_BLOWING].aggtype = AGG_TYPE_SUM;
  out_data[OUT_SUB_CANOP].aggtype = AGG_TYPE_SUM;
  out_data[OUT_SUB_SNOW].aggtype = AGG_TYPE_SUM;
  out_data[OUT_SUB_SURFACE].aggtype = AGG_TYPE_SUM;
  out_data[OUT_TRANSP_VEG].aggtype = AGG_TYPE_SUM;
  out_data[OUT_DELTACC_BAND].aggtype = AGG_TYPE_SUM;
  out_data[OUT_SNOW_MELT].aggtype = AGG_TYPE_SUM;
  out_data[OUT_SNOWT_FBFLAG].aggtype = AGG_TYPE_SUM;
  out_data[OUT_SOILT_FBFLAG].aggtype = AGG_TYPE_SUM;
  out_data[OUT_SURFT_FBFLAG].aggtype = AGG_TYPE_SUM;
  out_data[OUT_TCAN_FBFLAG].aggtype = AGG_TYPE_SUM;
  out_data[OUT_TFOL_FBFLAG].aggtype = AGG_TYPE_SUM;

  out_data[OUT_AREA_BAND].aggtype = AGG_TYPE_END;
  out_data[OUT_ELEV_BAND].aggtype = AGG_TYPE_END;
  out_data[OUT_AREA].aggtype = AGG_TYPE_END;

  out_data[OUT_GLAC_WAT_STOR].aggtype = AGG_TYPE_END;
  out_data[OUT_GLAC_AREA].aggtype = AGG_TYPE_END;
  out_data[OUT_GLAC_MBAL].aggtype = AGG_TYPE_SUM;
  out_data[OUT_GLAC_IMBAL].aggtype = AGG_TYPE_SUM;
  out_data[OUT_GLAC_ACCUM].aggtype = AGG_TYPE_SUM;
  out_data[OUT_GLAC_MELT].aggtype = AGG_TYPE_SUM;
  out_data[OUT_GLAC_SUB].aggtype = AGG_TYPE_SUM;
  out_data[OUT_GLAC_INFLOW].aggtype = AGG_TYPE_SUM;
  out_data[OUT_GLAC_OUTFLOW].aggtype = AGG_TYPE_SUM;
  out_data[OUT_GLAC_SURF_TEMP].aggtype = AGG_TYPE_END;
  out_data[OUT_GLAC_TSURF_FBFLAG].aggtype = AGG_TYPE_END;
  out_data[OUT_GLAC_OUTFLOW_COEF].aggtype = AGG_TYPE_END;
  out_data[OUT_GLAC_WAT_STOR_BAND].aggtype = AGG_TYPE_END;
  out_data[OUT_GLAC_AREA_BAND].aggtype = AGG_TYPE_END;
  out_data[OUT_GLAC_MBAL_BAND].aggtype = AGG_TYPE_SUM;
  out_data[OUT_GLAC_IMBAL_BAND].aggtype = AGG_TYPE_SUM;
  out_data[OUT_GLAC_ACCUM_BAND].aggtype = AGG_TYPE_SUM;
  out_data[OUT_GLAC_MELT_BAND].aggtype = AGG_TYPE_SUM;
  out_data[OUT_GLAC_SUB_BAND].aggtype = AGG_TYPE_SUM;
  out_data[OUT_GLAC_INFLOW_BAND].aggtype = AGG_TYPE_SUM;
  out_data[OUT_GLAC_OUTFLOW_BAND].aggtype = AGG_TYPE_SUM;

  // Allocate space for data
  for (v=0; v<N_OUTVAR_TYPES; v++) {
    out_data[v].data = new double[out_data[v].nelem];
    out_data[v].aggdata = new double[out_data[v].nelem];
  }

  // Initialize data values
  init_output_list(out_data, FALSE, "%.4f", OUT_TYPE_FLOAT, 1);

  return out_data;

}

out_data_file_struct::out_data_file_struct() : fh(NULL), varid(NULL) {

}

out_data_file_struct::~out_data_file_struct() {
  if (varid != NULL) {
    free(varid);
  }
}

void copy_data_file_format(const out_data_file_struct* out_template, std::vector<out_data_file_struct*>& list, const ProgramState* state) {
  for (int i = 0; i < state->options.Noutfiles; i++) {
    out_data_file_struct* curData = new out_data_file_struct();
    curData->fh = NULL;
    strncpy(curData->filename, out_template[i].filename, MAXSTRING);
    strncpy(curData->prefix, out_template[i].prefix, OUT_DATA_FILE_STRUCT_PREFIX_LENGTH);
    curData->nvars = out_template[i].nvars;
    curData->varid = (int *)calloc(curData->nvars, sizeof(int));
    for (int curVar = 0; curVar < curData->nvars; curVar++) {
      curData->varid[curVar] = out_template[i].varid[curVar];
    }
    list.push_back(curData);
  }
}

void init_output_list(OutputData *out_data, int write, const char *format, int type, float mult) {
/*************************************************************
  init_output_list()      Ted Bohn     September 08, 2006

  This routine initializes the output information for all output variables.

*************************************************************/
  int varid, i;

  for (varid=0; varid<N_OUTVAR_TYPES; varid++) {
    out_data[varid].write = write;
    out_data[varid].format = format;
    out_data[varid].type = type;
    out_data[varid].mult = mult;
    for(i=0; i<out_data[varid].nelem; i++) {
      out_data[varid].data[i] = 0;
    }
  }

}

int set_output_var(out_data_file_struct *out_data_files,
	                    int write,
	                    int filenum,
	                    OutputData *out_data,
	                    std::string varname,
	                    int varnum,
											std::string format,
	                    int type,
	                    float mult) {
/*************************************************************
  set_output_var()      Ted Bohn     September 08, 2006

  This routine updates the output information for a given output variable.

*************************************************************/
  int varid;
  int found=FALSE;
  int status=0;

  for (varid=0; varid<N_OUTVAR_TYPES; varid++) {
    if (out_data[varid].varname == varname) {
      found = TRUE;
      out_data[varid].write = write;
      if (format != "*")
        out_data[varid].format = format;
      if (type != 0)
        out_data[varid].type = type;
      if (mult != 0)
        out_data[varid].mult = mult;
      out_data_files[filenum].varid[varnum] = varid;
    }
  }
  if (!found) {
    status = -1;
    std::cout << "Error: set_output_var: " << varname << " was not found in the list of supported output variable names.  Please use the exact name listed in vicNl_def.h.\n";

  }
  return status;

}


void zero_output_list(OutputData *out_data) {
/*************************************************************
  zero_output_list()      Ted Bohn     September 08, 2006

  This routine resets the values of all output variables to 0.

*************************************************************/
  int varid, i;

  for (varid=0; varid<N_OUTVAR_TYPES; varid++) {
    for(i=0; i<out_data[varid].nelem; i++) {
      out_data[varid].data[i] = 0;
    }
  }

}
