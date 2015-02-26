/*
 * variable_mapping.c
 *
 *  Created on: Feb 24, 2015
 *      Author: mfischer
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vicNl.h"

/* Initializes mapping of input forcing variable FORCE_TYPE names in global param file to their names in the NetCDF forcing file.
 * Source of this list comes from VIC documentation at http://www.hydro.washington.edu/Lettenmaier/Models/VIC/Documentation/InputVarList.shtml
 */
void ProgramState:: build_forcing_variable_mapping() {

	forcing_mapping = { // TODO: get updates for UPPERCASE values on right from Markus
		{"AIR_TEMP", "ta"},
		{"ALBEDO", "ALBEDO"},
		{"CHANNEL_IN", "CHANNEL_IN"},
		{"CAT_M", "CAT_M"},
		{"CRAINF", "CRAINF"},
		{"CSNOWF", "CSNOWF"},
		{"DENSITY", "DENSITY"},
		{"FDIR", "FDIR"},
		{"LONGWAVE", "LONGWAVE"},
		{"LSRAINF", "LSRAINF"},
		{"LSSNOWF", "LSSNOWF"},
		{"PREC", "pr"},
		{"PRESSURE", "PRESSURE"},
		{"QAIR", "QAIR"},
		{"RAINF", "RAINF"},
		{"REL_HUMID", "REL_HUMID"},
		{"SHORTWAVE", "SHORTWAVE"},
		{"SNOWF", "SNOWF"},
		{"TMAX", "tasmax"},
		{"TMIN", "tasmin"},
		{"TSKC", "TSKC"},
		{"VEGCOVER", "VEGCOVER"},
		{"VP", "OUT_VP"},
		{"WIND", "wind"},
		{"WIND_E", "WIND_E"},
		{"WIND_N", "WIND_N"}
	};
}

// Changes value for a given forcing_mapping key. VIC will look for this variable name in the input forcing file instead.
void ProgramState:: set_forcing_variable_name(std::string variableKey, std::string newName){
	if (forcing_mapping.find(variableKey) == forcing_mapping.end()) {
		        throw VICException("Error: set_forcing_variable_name could not find variable in forcing_mapping: " + variableKey);
		}
		forcing_mapping.at(variableKey) = newName;
}


// Initializes default output variable naming and metadata mapping
void ProgramState:: build_output_variable_mapping() {
	// TODO: get updates for UPPERCASE name values from Markus
	output_mapping = {
		{"OUT_PREC", VariableMetaData("kg m-2 s-1", "PREC", "precipitation_flux", "Precipitation", "time: mean")},
		{"OUT_EVAP",       VariableMetaData("kg m-2 s-1", "EVAP", "water_evaporation_flux", "Evaporation", "time: mean")},
		{"OUT_RUNOFF",     VariableMetaData("kg m-2 s-1", "RUNOFF", "runoff_flux", "Total Runoff", "time: mean area: mean where land ")},
		{"OUT_BASEFLOW",     VariableMetaData("mm", "BASEFLOW", "", "baseflow out of the bottom layer", "")},
		{"OUT_WDEW",     VariableMetaData("mm", "WDEW", "", "total moisture interception storage in canopy", "")},
		{"OUT_SOIL_LIQ",   VariableMetaData("kg m-2", "SOIL_LIQ", "moisture_content_of_soil_layer", "Water Content of Soil Layer", "time: mean area: mean where land")},
		{"OUT_RAD_TEMP",   VariableMetaData("K s-1", "RAD_TEMP", "tendency_of_air_temperature_due_to_radiative_heating", "Tendency of Air Temperature due to Radiative Heating", "time: point")},
		{"OUT_NET_SHORT",  VariableMetaData("W m-2", "NET_SHORT", "net_downward_shortwave_flux_at_sea_water_surface", "Net Downward Shortwave Radiation at Sea Water Surface", "time: mean area: mean where sea")},
		{"OUT_R_NET",      VariableMetaData("W m-2", "R_NET", "net_downward_radiative_flux_at_top_of_atmosphere_model", "Net Downward Flux at Top of Model", "time: mean")},
		{"OUT_LATENT",     VariableMetaData("W m-2", "LATENT", "surface_upward_latent_heat_flux", "Surface Upward Latent Heat Flux", "time: mean")},
		{"OUT_EVAP_CANOP", VariableMetaData("kg m-2 s-1", "EVAP_CANOP", "water_evaporation_flux_from_canopy", "Evaporation from Canopy", "time: mean area: mean where land")},
		{"OUT_TRANSP_VEG", VariableMetaData("kg m-2 s-1", "TRANSP_VEG", "transpiration_flux", "Transpiration", "time: mean area: mean where land")},
		{"OUT_EVAP_BARE",  VariableMetaData("kg m-2 s-1", "EVAP_BARE", "water_evaporation_flux_from_soil", "Water Evaporation from Soil", "time: mean area: mean where land")},
		{"OUT_SUB_CANOP",    VariableMetaData("mm", "SUB_CANOP", "", "net sublimation from snow stored in canopy", "")},
		{"OUT_SUB_SNOW",   VariableMetaData("kg m-2 s-1", "SUB_SNOW", "surface_snow_and_ice_sublimation_flux", "Surface Snow and Ice Sublimation Flux", "time: mean")},
		{"OUT_SENSIBLE",   VariableMetaData("W m-2", "SENSIBLE", "surface_upward_sensible_heat_flux", "Surface Upward Sensible Heat Flux", "time: mean")},
		{"OUT_GRND_FLUX",  VariableMetaData("W m-2", "GRND_FLUX", "surface_downward_latent_heat_flux", "Surface Downward Latent Heat Flux", "time: mean area: mean where ice_free_sea over sea")},
		{"OUT_DELTAH",     VariableMetaData("W/m2", "DELTAH", "", "rate of change in heat storage", "")},
		{"OUT_FUSION",     VariableMetaData("W/m2", "FUSION", "", "net energy used to melt/freeze soil moisture", "")},
		{"OUT_AERO_RESIST",    VariableMetaData("s/m", "AERO_RESIST", "", "\"scene\"canopy aerodynamic resistance", "")},
		{"OUT_SURF_TEMP",  VariableMetaData("K", "SURF_TEMP", "surface_temperature", "Surface Temperature", "time: mean")},
		{"OUT_ALBEDO",     VariableMetaData("%", "ALBEDO", "", "average surface albedo", "")},
		{"OUT_REL_HUMID",  VariableMetaData("%", "REL_HUMID", "relative_humidity", "Relative Humidity", "time: mean")},
		{"OUT_IN_LONG",    VariableMetaData("W m-2", "IN_LONG", "surface_downwelling_longwave_flux_in_air", "Surface Downwelling Longwave Radiation", "time: mean")},
		{"OUT_AIR_TEMP",   VariableMetaData("K", "AIR_TEMP", "air_temperature", "Air Temperature", "time: mean")},
		{"OUT_WIND",       VariableMetaData("m s-1", "WIND", "wind_speed", "Daily-Mean Near-Surface Wind Speed", "time: mean")},
		{"OUT_SWE",        VariableMetaData("kg m-2 s-1", "SWE", "snowfall_flux", "Snowfall Flux", "time: mean")},
		{"OUT_SNOW_DEPTH", VariableMetaData("m", "SNOW_DEPTH", "surface_snow_thickness", "Snow Depth", "time: mean area: mean where land")},
		{"OUT_SNOW_CANOPY", VariableMetaData("kg m-2", "SNOW_CANOPY", "veg_snow_amount", "Canopy Snow Amount", "time: mean area: mean where land")},
		{"OUT_SNOW_COVER", VariableMetaData("%", "SNOW_COVER", "surface_snow_area_fraction", "Snow Area Fraction", "time: mean")},
		{"OUT_ADVECTION",    VariableMetaData("W/m2", "ADVECTION", "", "advected energy", "")},
		{"OUT_DELTACC",    VariableMetaData("W/m2", "DELTACC", "", "rate of change in cold content in snow pack", "")},
		{"OUT_SNOW_FLUX",  VariableMetaData("W m-2", "SNOW_FLUX", "surface_downward_heat_flux_in_snow", "Downward Heat Flux into Snow Where Land over Land", "time: mean area: mean where land")},
		{"OUT_RFRZ_ENERGY",    VariableMetaData("W/m2", "RFRZ_ENERGY", "", "net energy used to refreeze liquid water in snowpack", "")},
		{"OUT_MELT_ENERGY",    VariableMetaData("W/m2", "MELT_ENERGY", "", "energy of fusion (melting) in snowpack", "")},
		{"OUT_ADV_SENS",     VariableMetaData("W/m2", "ADV_SENS", "", "net sensible flux advected to snow pack", "")},
		{"OUT_LATENT_SUB", VariableMetaData("W m-2", "LATENT_SUB", "surface_upward_latent_heat_flux", "Surface Upward Latent Heat Flux", "time: mean")},
		{"OUT_SNOW_SURF_TEMP", VariableMetaData("K", "SNOW_SURF_TEMP", "temperature_in_surface_snow", "Snow Internal Temperature", "time: mean (with samples weighted by snow mass) area: mean where land")},
		{"OUT_SNOW_PACK_TEMP",     VariableMetaData("C", "SNOW_PACK_TEMP", "", "snow pack temperature", "")},
		{"OUT_SNOW_MELT",  VariableMetaData("kg m-2 s-1", "SNOW_MELT", "surface_snow_melt_flux", "Surface Snow Melt", "time: mean area: mean where land")},
		{"OUT_SUB_BLOWING", VariableMetaData("kg m-2 s-1", "SUB_BLOWING", "surface_snow_and_ice_sublimation_flux", "Blowing Snow Sublimation Flux", "time: mean area: mean where land")},
		{"OUT_SUB_SURFACE", VariableMetaData("kg m-2 s-1", "SUB_SURFACE", "surface_snow_and_ice_sublimation_flux", "Surface Snow and Ice Sublimation Flux", "time: mean area: mean where land")},
		{"OUT_SUB_SNOW",     VariableMetaData("mm", "SUB_SNOW", "", "total net sublimation from snow pack (surface and blowing)", "")},
		{"OUT_FDEPTH",     VariableMetaData("cm", "FDEPTH", "", "depth of freezing fronts", "")},
		{"OUT_TDEPTH",     VariableMetaData("cm", "TDEPTH", "", "depth of thawing fronts", "")},
		{"OUT_SOIL_MOIST", VariableMetaData("kg m-2", "SOIL_MOIST", "soil_moisture_content", "Total Soil Moisture Content", "time: mean area: mean where land")},
		{"OUT_SURF_FROST_FRAC",  VariableMetaData("%", "SURF_FROST_FRAC", "surface_snow_area_fraction", "Snow Area Fraction", "time: mean")},
		{"OUT_SWE_BAND",         VariableMetaData("kg m-2", "SWE_BAND", "liquid_water_content_of_snow_layer", "Liquid Water Content of Snow Layer", "time: mean area: mean where land")},
		{"OUT_SNOW_DEPTH_BAND",  VariableMetaData("m", "SNOW_DEPTH_BAND", "surface_snow_thickness", "Snow Depth", "time: mean area: mean where land", 0.01, 0, true)},
		{"OUT_SNOW_CANOPY_BAND",     VariableMetaData("mm", "SNOW_CANOPY_BAND", "", "snow interception storage in canopy", "")},
		{"OUT_ADVECTION_BAND",     VariableMetaData("W/m2", "ADVECTION_BAND", "", "advected energy", "")},
		{"OUT_DELTACC_BAND",     VariableMetaData("W/m2", "DELTACC_BAND", "", "change in cold content in snow pack", "")},
		{"OUT_SNOW_FLUX_BAND",     VariableMetaData("W/m2", "SNOW_FLUX_BAND", "", "energy flux through snow pack", "")},
		{"OUT_RFRZ_ENERGY_BAND",     VariableMetaData("W/m2", "RFRZ_ENERGY_BAND", "", "net energy used to refreeze liquid water in snowpack", "")},
		{"OUT_NET_SHORT_BAND",     VariableMetaData("W/m2", "NET_SHORT_BAND", "", "net downward shortwave flux", "")},
		{"OUT_NET_LONG_BAND",    VariableMetaData("W/m2", "NET_LONG_BAND", "", "net downward longwave flux", "")},
		{"OUT_ALBEDO_BAND",    VariableMetaData("%", "ALBEDO_BAND", "", "average surface albedo", "")},
		{"OUT_LATENT_BAND",    VariableMetaData("W/m2", "LATENT_BAND", "", "net upward latent heat flux", "")},
		{"OUT_SENSIBLE_BAND",    VariableMetaData("W/m2", "SENSIBLE_BAND", "", "net upward sensible heat flux", "")},
		{"OUT_GRND_FLUX_BAND",     VariableMetaData("W/m2", "GRND_FLUX_BAND", "", "net heat flux into ground", "")},
		{"OUT_LAKE_ICE_TEMP",    VariableMetaData("C", "LAKE_ICE_TEMP", "", "temperature of lake ice", "")},
		{"OUT_LAKE_ICE_HEIGHT",    VariableMetaData("cm", "LAKE_ICE_HEIGHT", "", "thickness of lake ice", "")},
		{"OUT_LAKE_ICE_FRACT",     VariableMetaData("%", "LAKE_ICE_FRACT", "", "fractional coverage of lake ice", "")},
		{"OUT_LAKE_DEPTH",     VariableMetaData("m", "LAKE_DEPTH", "", "lake depth (distance between surface and deepest point)", "")},
		{"OUT_LAKE_SURF_AREA",     VariableMetaData("m2", "LAKE_SURF_AREA", "", "lake surface area", "")},
		{"OUT_LAKE_VOLUME",    VariableMetaData("m3", "LAKE_VOLUME", "", "lake volume", "")},
		{"OUT_LAKE_SURF_TEMP",     VariableMetaData("C", "LAKE_SURF_TEMP", "", "lake surface temperature", "")},
		{"OUT_LAKE_EVAP",    VariableMetaData("mm", "LAKE_EVAP", "", "net evaporation from lake surface", "")},
		{"OUT_GLAC_WAT_STOR",    VariableMetaData("mm", "GLAC_WAT_STOR", "", "glacier water storage", "")},
		{"OUT_GLAC_AREA",    VariableMetaData("%", "GLAC_AREA", "", "glacier surface area fraction", "")},
		{"OUT_GLAC_MBAL",    VariableMetaData("mm", "GLAC_MBAL", "", "glacier mass balance", "")},
		{"OUT_GLAC_IMBAL",     VariableMetaData("mm", "GLAC_IMBAL", "", "glacier ice mass balance", "")},
		{"OUT_GLAC_ACCUM",     VariableMetaData("mm", "GLAC_ACCUM", "", "glacier ice accumulation from conversion of firn to ice", "")},
		{"OUT_GLAC_MELT",    VariableMetaData("mm", "GLAC_MELT", "", "glacier ice melt", "")},
		{"OUT_GLAC_SUB",     VariableMetaData("mm", "GLAC_SUB", "", "Net sublimation of glacier ice", "")},
		{"OUT_GLAC_INFLOW",    VariableMetaData("mm", "GLAC_INFLOW", "", "glacier water inflow from snow melt, ice melt and rainfall", "")},
		{"OUT_GLAC_OUTFLOW",     VariableMetaData("mm", "GLAC_OUTFLOW", "", "glacier water outflow", "")},
		{"OUT_GLAC_SURF_TEMP",     VariableMetaData("C", "GLAC_SURF_TEMP", "", "glacier surface temperature", "")},
		{"OUT_GLAC_TSURF_FBFLAG",    VariableMetaData("%", "GLAC_TSURF_FBFLAG", "", "glacier surface temperature flag", "")},
		{"OUT_GLAC_DELTACC",     VariableMetaData("W/m2", "GLAC_DELTACC", "", "rate of change of cold content in glacier surface layer", "")},
		{"OUT_GLAC_FLUX",    VariableMetaData("W/m2", "GLAC_FLUX", "", "energy flux through glacier surface layer", "")},
		{"OUT_GLAC_OUTFLOW_COEF",    VariableMetaData("%", "GLAC_OUTFLOW_COEF", "", "glacier outflow coefficient", "")},
		{"OUT_GLAC_DELTACC_BAND",    VariableMetaData("W/m2", "GLAC_DELTACC_BAND", "", "rate of change of cold content in glacier surface layer", "")},
		{"OUT_GLAC_FLUX_BAND",     VariableMetaData("W/m2", "GLAC_FLUX_BAND", "", "energy flux through glacier surface layer", "")},
		{"OUT_GLAC_WAT_STOR_BAND",     VariableMetaData("mm", "GLAC_WAT_STOR_BAND", "", "glacier water storage", "")},
		{"OUT_GLAC_AREA_BAND",     VariableMetaData("%", "GLAC_AREA_BAND", "", "glacier surface area fraction", "")},
		{"OUT_GLAC_MBAL_BAND",     VariableMetaData("mm", "GLAC_MBAL_BAND", "", "glacier mass balance", "")},
		{"OUT_GLAC_IMBAL_BAND",    VariableMetaData("mm", "GLAC_IMBAL_BAND", "", "glacier ice mass balance", "")},
		{"OUT_GLAC_ACCUM_BAND",    VariableMetaData("mm", "GLAC_ACCUM_BAND", "", "glacier ice accumulation from conversion of firn to ice", "")},
		{"OUT_GLAC_MELT_BAND",     VariableMetaData("mm", "GLAC_MELT_BAND", "", "glacier ice melt", "")},
		{"OUT_GLAC_SUB_BAND",    VariableMetaData("mm", "GLAC_SUB_BAND", "", "Net sublimation of glacier ice", "")},
		{"OUT_GLAC_INFLOW_BAND",     VariableMetaData("mm", "GLAC_INFLOW_BAND", "", "glacier water inflow from snow melt, ice melt and rainfall", "")},
		{"OUT_GLAC_OUTFLOW_BAND",    VariableMetaData("mm", "GLAC_OUTFLOW_BAND", "", "glacier water outflow", "")},
		{"OUT_SHORTWAVE",    VariableMetaData("mm", "SHORTWAVE", "", "glacier SHORTWAVE", "")},
		{"OUT_LONGWAVE",    VariableMetaData("mm", "LONGWAVE", "", "glacier LONGWAVE", "")},
		{"OUT_DENSITY",    VariableMetaData("mm", "DENSITY", "", "glacier DENSITY", "")},
		{"OUT_VP",    VariableMetaData("mm", "VP", "", "glacier VP", "")},
		{"OUT_PRESSURE",    VariableMetaData("mm", "PRESSURE", "", "glacier PRESSURE", "")},
		{"OUT_RAINF",       VariableMetaData("kg m-2 s-1", "RAINF", "rainfall_flux", "Rainfall", "time: mean")},
		{"OUT_SNOWF",       VariableMetaData("kg m-2 s-1", "SNOWF", "snowfall_flux", "Snowfall", "time: mean")},
		{"OUT_INFLOW",    VariableMetaData("mm", "INFLOW", "", "moisture that reaches top of soil column", "")},
		{"OUT_WATER_ERROR",    VariableMetaData("mm", "WATER_ERROR", "", "water budget error", "")},
		{"OUT_AERO_RESIST1",    VariableMetaData("s/m", "AERO_RESIST1", "", "surface aerodynamic resistance", "")},
		{"OUT_AERO_RESIST2",    VariableMetaData("s/m", "AERO_RESIST2", "", "overstory aerodynamic resistance", "")}
	};
}

// Changes value of name member for a given output_mapping variable key.  VIC will name the variable accordingly in the output file.
void ProgramState:: set_output_variable_name(std::string variableKey, std::string newName) {

	if (output_mapping.find(variableKey) == output_mapping.end()) {
	        throw VICException("Error: set_output_variable_name could not find variable in output_mapping: " + variableKey);
	}
	output_mapping.at(variableKey).name = newName;
}



