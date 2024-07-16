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

/* Initializes mapping of input forcing variable FORCE_TYPE names in global parameter file to their names in the NetCDF forcing file.
 * Source of this list comes from VIC documentation at http://www.hydro.washington.edu/Lettenmaier/Models/VIC/Documentation/InputVarList.shtml
 */
void ProgramState:: build_forcing_variable_mapping() {

	forcing_mapping = {
		{"AIR_TEMP", "tas"}, //air temperature
		{"ALBEDO", "ALBEDO"}, //surface albedo
		{"CHANNEL_IN", "CHANNEL_IN"}, //incoming channel flow
		{"CAT_M", "CAT_M"}, //atmospheric CO2 mixing ratio
		{"CRAINF", "CRAINF"}, //convective rainfall
		{"CSNOWF", "CSNOWF"}, //convective snowfall
		{"DENSITY", "DENSITY"}, //atmospheric density
		{"FDIR", "FDIR"}, //fraction of incoming shortwave that is direct
		{"LONGWAVE", "rlds"}, //incoming longwave radiation
		{"LSRAINF", "LSRAINF"}, //large-scale rainfall
		{"LSSNOWF", "LSSNOWF"}, //large-scale snowfall
		{"PREC", "pr"}, //total precipitation
		{"PRESSURE", "ps"}, //atmospheric pressure
		{"QAIR", "huss"}, //specific humidity
		{"RAINF", "RAINF"}, //rainfall
		{"REL_HUMID", "rhs"}, //relative humidity
		{"SHORTWAVE", "rsds"}, //incoming shortwave radiation
		{"SNOWF", "prs"}, //snowfall
		{"TMAX", "tasmax"}, //maximum daily air temperature
		{"TMIN", "tasmin"}, //minimum daily air temperature
		{"TSKC", "clt"}, //cloud cover
		{"VEGCOVER", "VEGCOVER"}, //partial vegetation cover fraction
		{"VP", "VP"}, //atmospheric vapor pressures
		{"WIND", "wind"}, //wind speed
		{"WIND_E", "uas"}, //east component of wind speed
		{"WIND_N", "vas"} //north component of wind speed
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
	output_mapping = {
		//water balance terms - states
		{"OUT_ASAT",			VariableMetaData("1", "ASAT", "", "Saturated area fraction", "time: point area: sum")},
		{"OUT_LAKE_DEPTH",     		VariableMetaData("m", "LAKE_DEPTH", "", "Lake depth (distance between surface and deepest point)", "time: point area: point")},
		{"OUT_LAKE_ICE_FRACT",     	VariableMetaData("1", "LAKE_ICE_FRACT", "", "Fractional coverage of lake ice", "time: point area: sum")},
		{"OUT_LAKE_ICE_HEIGHT",    	VariableMetaData("cm", "LAKE_ICE_HEIGHT", "", "Thickness of lake ice", "time: point area: mean")},
		{"OUT_LAKE_SURF_AREA",    	VariableMetaData("m2", "LAKE_SURF_AREA", "", "Lake surface area", "time: point area: sum")},
		{"OUT_LAKE_VOLUME",    		VariableMetaData("m3", "LAKE_VOLUME", "", "Lake volume", "time: point area: sum")},
		{"OUT_ROOTMOIST",		VariableMetaData("mm", "ROOTMOIST", "", "Root zone soil moisture", "time: point area: mean")},
		{"OUT_SMFROZFRAC",		VariableMetaData("1", "SMFROZFRAC", "", "Fraction of soil moisture (by mass) that is ice, for each soil layer", "time: point area: mean")},
		{"OUT_SMLIQFRAC",		VariableMetaData("1", "SMLIQFRAC", "", "Fraction of soil moisture (by mass) that is liquid, for each soil layer", "time: point area: mean")},
		{"OUT_SNOW_CANOPY", 		VariableMetaData("mm", "SNOW_CANOPY", "", "Snow interception storage in canopy", "time: point area: mean")},
		{"OUT_SNOW_COVER", 		VariableMetaData("1", "SNOW_COVER", "surface_snow_area_fraction", "Snow area fraction", "time: point area: sum")},
		{"OUT_SNOW_DEPTH", 		VariableMetaData("cm", "SNOW_DEPTH", "surface_snow_thickness", "Snow depth", "time: point area: mean")},
		{"OUT_SOIL_ICE",		VariableMetaData("mm", "SOIL_ICE", "lwe_thickenss_of_frozen_water_content_of_soil_layer", "Soil ice content of soil layer", "time: point area: mean")},
		{"OUT_SOIL_ICE_TOT", 		VariableMetaData("mm", "SOIL_ICE_TOT", "soil_frozen_water_content", "Total soil ice content", "time: point area: mean")},
		{"OUT_SOIL_LIQ",   		VariableMetaData("mm", "SOIL_LIQ", "lwe_thickness_of_liquid_water_content_of_soil_layer", "Water content of soil layer", "time: point area: mean")},
		{"OUT_SOIL_LIQ_TOT",		VariableMetaData("mm", "SOIL_LIQ_TOT", "", "Total water content of all soil layers", "time: point area: mean")},
		{"OUT_SOIL_MOIST", 		VariableMetaData("mm", "SOIL_MOIST", "lwe_thickness_of_soil_moisture_content", "Soil total moisture content for soil layer", "time: point area: mean")},
		{"OUT_SOIL_MOIST_TOT",		VariableMetaData("mm", "SOIL_MOIST_TOT", "soil_moisture_content", "Total soil moisture content", "time: point area: mean")},
		{"OUT_SOIL_WET", 		VariableMetaData("1", "SOIL_WET", "", "Vertical average of (soil moisture - wilting point)/(maximum soil moisture - wilting point)", "time: point area: mean")},
		{"OUT_SURFSTOR",		VariableMetaData("mm", "SURFSTOR", "", "Storage of liquid water and ice (not snow) on surface (ponding)", "time: point area: mean")},
		{"OUT_SURF_FROST_FRAC",  	VariableMetaData("1", "SURF_FROST_FRAC", "", "Fraction of soil surface that is frozen", "time: point area: sum")},
		{"OUT_SWE",       		VariableMetaData("mm", "SWE", "lwe_thickness_of_surface_snow_amount", "Snow water equivalent in snow pack (including vegetation-intercepted snow)", "time: point area: mean")},
		{"OUT_WDEW",     		VariableMetaData("mm", "WDEW", "lwe_thickness_of_canopy_water_amount", "Total moisture interception storage in canopy", "time: point area: mean")},
		{"ZWT",		 		VariableMetaData("cm", "ZWT", "", "Water table position- method 1 (zwt within lowest unsaturated layer)", "time: point area: mean")},
		{"ZWT2",		 	VariableMetaData("cm", "ZWT2", "", "Water table position- method 2 (zwt of total moisture across top-most N-1 layers, lumped together)", "time: point area: mean")},
		{"ZWT3",		 	VariableMetaData("cm", "ZWT3", "", "Water table position- method 3 (zwt of total moisture across all layers, lumped together)", "time: point area: mean")},
		{"ZWTL",		 	VariableMetaData("cm", "ZWTL", "", "Water table positions per soil layer", "time: point area: mean")},
		//water balance terms - fluxes
		{"OUT_BASEFLOW",  		VariableMetaData("mm", "BASEFLOW", "lwe_thickness_of_baseflow_amount", "Baseflow out of the bottom layer", "time: mean area: mean")},
		{"OUT_DELINTERCEPT",		VariableMetaData("mm", "DELINTERCEPT", "", "Change in canopy interception storage", "time: mean area: mean")},
		{"OUT_DELSOILMOIST",		VariableMetaData("mm", "DELSOILMOIST", "", "Change in soil water content", "time: mean area: mean")},
		{"OUT_DELSURFSTOR",		VariableMetaData("mm", "DELSURFSTOR", "", "Change in snow water equivalent", "time: mean area: mean")},
		{"OUT_DELSWE",			VariableMetaData("mm", "DELSWE", "", "Change in surface liquid water storage", "time: mean area: mean")},
		{"OUT_EVAP", 			VariableMetaData("mm", "EVAP", "lwe_thickness_of_water_evaporation_amount", "Total net evaporation", "time: mean area: mean")},
		{"OUT_EVAP_BARE",  		VariableMetaData("mm", "EVAP_BARE", "lwe_thickness_of_water_evaporation_amount_from_soil", "Net evaporation from bare soil", "time: mean area: mean")},
		{"OUT_EVAP_CANOP", 		VariableMetaData("mm", "EVAP_CANOP", "lwe_thickness_of_water_evaporation_amount_from_canopy", "Net evaporation from canopy interception", "time: mean area: mean")},
		{"OUT_INFLOW",    		VariableMetaData("mm", "INFLOW", "", "Moisture that reaches top of soil column", "time: mean area: mean")},
		{"OUT_LAKE_EVAP",    		VariableMetaData("mm", "LAKE_EVAP", "", "Net evaporation from lake surface", "time: mean area: mean")},
		{"OUT_PET_H2OSURF",		VariableMetaData("mm", "PET_H20SURF", "lwe_thickness_of_water_potential_evaporation_amount", "Potential evap from open water", "time: mean area: mean")},
		{"OUT_PET_NATVEG",		VariableMetaData("mm", "PET_NATVEG", "lwe_thickness_of_water_potential_evaporation_amount", "Potential evap (transpiration only) from current vegetation and current canopy resistance", "time: mean area: mean")},
		{"OUT_PET_SATSOIL",		VariableMetaData("mm", "PET_SATSOIL", "lwe_thickness_of_water_potential_evaporation_amount", "Potential evap from saturated bare soil", "time: mean area: mean")},
		{"OUT_PET_SHORT",		VariableMetaData("mm", "PET_SHORT", "lwe_thickness_of_water_potential_evaporation_amount", "Potential evap (transpiration only) from short reference crop (grass) ", "time: mean area: mean")},
		{"OUT_PET_TALL",		VariableMetaData("mm", "PET_TALL", "lwe_thickness_of_water_potential_evaporation_amount", "Potential evap (transpiration only) from tall reference crop (alfalfa)", "time: mean area: mean")},
		{"OUT_PET_VEGNOCR",		VariableMetaData("mm", "PET_VEGNOCR", "lwe_thickness_of_water_potential_evaporation_amount", "Potential evap (transpiration only) from current vegetation and no canopy resistance", "time: mean area: mean")},
		{"OUT_PREC",			VariableMetaData("mm", "PREC", "lwe_thickness_of_precipitation_amount", "Precipitation", "time: mean area: mean")},
		{"OUT_RAINF",       		VariableMetaData("mm", "RAINF", "thickness_of_rainfall_amount", "Rainfall", "time: mean area: mean")},
		{"OUT_REFREEZE",		VariableMetaData("mm", "REFREEZE", "", "Refreezing of water in the snow", "time: mean area: mean")},
		{"OUT_RUNOFF",     		VariableMetaData("mm", "RUNOFF", "thickness_of_surface_runoff_amount", "Surface runoff", "time: mean area: mean")},
		{"OUT_RUNOFF_SNOW",		VariableMetaData("mm", "RUNOFF_SNOW", "thickness_of_surface_runoff_amount_from_snow_and_ice_melt", "Surface runoff from snow and ice melt", "time: mean area: mean")},
		{"OUT_SNOW_MELT",  		VariableMetaData("mm", "SNOW_MELT", "thickness_of_surface_snow_melt_amount", "Snow melt", "time: mean area: mean")},
		{"OUT_SNOWF",       		VariableMetaData("mm", "SNOWF", "lwe_thickness_of_snowfall_amount", "Snowfall", "time: mean area: mean")},
		{"OUT_SUB_BLOWING", 		VariableMetaData("mm", "SUB_BLOWING", "", "Net sublimation of blowing snow", "time: mean area: mean")},
		{"OUT_SUB_CANOP",   		VariableMetaData("mm", "SUB_CANOP", "", "Net sublimation from snow stored in canopy", "time: mean area: mean")},
		{"OUT_SUB_SNOW",   		VariableMetaData("mm", "SUB_SNOW", "lwe_thickness_of_surface_snow_sublimation_amount", "Total net sublimation from snow pack (surface and blowing)", "time: mean area: mean")},
		{"OUT_SUB_SURFACE",	 	VariableMetaData("mm", "SUB_SURFACE", "", "Net sublimation from snow pack surface", "time: mean area: mean")},
		{"OUT_TRANSP_VEG", 		VariableMetaData("mm", "TRANSP_VEG", "thickness_of_transpiration_amount", "Transpiration", "time: mean area: mean")},
		{"OUT_WATER_ERROR",    		VariableMetaData("mm", "WATER_ERROR", "", "Water budget error", "time: mean area: mean")},
		//energy balance terms - state variables
		{"OUT_ALBEDO",     		VariableMetaData("1", "ALBEDO", "surface_albedo", "Average surface albedo", "time: point area: mean")},
		{"OUT_BARESOILT",		VariableMetaData("degree_Celsius", "BARESOILT", "", "Bare soil surface temperature", "time: point area: mean")},
		{"OUT_FDEPTH",     		VariableMetaData("cm", "FDEPTH", "", "Depth of freezing fronts for each freezing front", "time: point area: mean")},
		{"OUT_LAKE_ICE_TEMP",    	VariableMetaData("degree_Celsius", "LAKE_ICE_TEMP", "", "Temperature of lake ice", "time: point area: mean")},
		{"OUT_LAKE_SURF_TEMP",     	VariableMetaData("degree_Celsius", "LAKE_SURF_TEMP", "", "Lake surface temperature", "time: point area: mean")},
		{"OUT_RAD_TEMP",   		VariableMetaData("K", "RAD_TEMP", "", "Average radiative surface temperature", "time: point area: mean")},
		{"OUT_SALBEDO",     		VariableMetaData("1", "SALBEDO", "surface_albedo_assuming_deep_snow", "Snow pack albedo", "time: point area: mean")},
		{"OUT_SNOW_PACK_TEMP", 		VariableMetaData("degree_Celsius", "SNOW_PACK_TEMP", "", "Snow pack temperature", "time: point area: mean")},
		{"OUT_SNOW_SURF_TEMP", 		VariableMetaData("degree_Celsius", "SNOW_SURF_TEMP", "surface_temperature_where_snow", "Snow surface temperature", "time: point area: mean")},
		{"OUT_SNOWT_FBFLAG",		VariableMetaData("", "SNOWT_FBFLAG", "", "Count of snow surface temperature fallback occurrences", "time: point area: mean")},
		{"OUT_SOIL_TEMP",		VariableMetaData("degree_Celsius", "SOIL_TEMP", "soil_tempearture", "Soil temperature for each layer", "time: point area: mean")},
		{"OUT_SOIL_TNODE",		VariableMetaData("degree_Celsius", "SOIL_TNODE", "", "Soil temperature for each thermal node", "time: point area: mean")},
		{"OUT_SOIL_TNODE_WL",		VariableMetaData("degree_Celsius", "SOIL_TNODE_WL", "", "Soil temperature for each thermal node in the wetland", "time: point area: mean")},
		{"OUT_SOILT_FBFLAG",		VariableMetaData("", "SOILT_FBFLAG", "", "Count of soil temperature fallback occurrences", "time: point area: mean")},
		{"OUT_SURF_TEMP", 		VariableMetaData("degree_Celsius", "SURF_TEMP", "surface_temperature", "Average surface temperature", "time: point area: mean")},
		{"OUT_SURFT_FBFLAG",		VariableMetaData("", "SURFT_FBFLAG", "", "Count of surface temperature fallback occurrences", "time: point area: mean")},
		{"OUT_TCAN_FBFLAG",		VariableMetaData("", "TCAN_FBFLAG", "", "Count of canopy temperature fallback occurrences", "time: point area: mean")},
		{"OUT_TDEPTH",     		VariableMetaData("cm", "TDEPTH", "", "Depth of thawing fronts", "time: point area: mean")},
		{"OUT_TFOL_FBFLAG",		VariableMetaData("", "TFOL_FBFLAG", "", "Count of foliage temperature fallback occurrences", "time: point area: mean")},
		{"OUT_VEGT",			VariableMetaData("degree_Celsius", "VEGT", "canopy_temperature", "Average vegetation canopy temperature", "time: point area: mean")},
		//energy balance terms - fluxes
		{"OUT_ADVECTION",   		VariableMetaData("W m-2", "ADVECTION", "", "Advected energy", "time: mean area: mean")},
		{"OUT_ADV_SENS",     		VariableMetaData("W m-2", "ADV_SENS", "", "Net sensible flux advected to snow pack", "time: mean area: mean")},
		{"OUT_DELTACC",    		VariableMetaData("W m-2", "DELTACC", "", "Rate of change in cold content in snow pack", "time: mean area: mean")},
		{"OUT_DELTAH",     		VariableMetaData("W m-2", "DELTAH", "", "Rate of change in heat storage", "time: mean area: mean")},
		{"OUT_ENERGY_ERROR",		VariableMetaData("W m-2", "ENERGY_ERROR", "", "Energy budget error", "time: mean area: mean")},
		{"OUT_FUSION",     		VariableMetaData("W m-2", "FUSION", "", "Net energy used to melt/freeze soil moisture", "time: mean area: mean")},
		{"OUT_GRND_FLUX",  		VariableMetaData("W m-2", "GRND_FLUX", "downward_heat_flux_at_ground_level_in_soil", "Net heat flux in to ground", "time: mean area: mean")},
		{"OUT_IN_LONG",    		VariableMetaData("W m-2", "IN_LONG", "", "Incoming longwave at ground surface (under vegetation) ", "time: mean area: mean")},
		{"OUT_LATENT",     		VariableMetaData("W m-2", "LATENT", "surface_upward_latent_heat_flux", "Net upward latent heat flux", "time: mean area: mean")},
		{"OUT_LATENT_SUB", 		VariableMetaData("W m-2", "LATENT_SUB", "", "Net upward latent heat flux from sublimation", "time: mean area: mean")},
		{"OUT_MELT_ENERGY",    		VariableMetaData("W m-2", "MELT_ENERGY", "surface_snow_melt_heat_flux", "Energy of fusion (melting) in snowpack", "time: mean area: mean")},
		{"OUT_NET_LONG", 		VariableMetaData("W m-2", "NET_LONG", "net_downward_longwave_flux_in_air", "Net downward longwave flux", "time: mean area: mean")},
		{"OUT_NET_SHORT", 		VariableMetaData("W m-2", "NET_SHORT", "net_downward_shortwave_flux_in_air", "Net downward shortwave flux", "time: mean area: mean")},
		{"OUT_R_NET",      		VariableMetaData("W m-2", "R_NET", "surface_net_downward_radiative_flux", "Net downward radiation flux", "time: mean area: mean")},
		{"OUT_RFRZ_ENERGY",    		VariableMetaData("W m-2", "RFRZ_ENERGY", "", "Net energy used to refreeze liquid water in snowpack", "time: mean area: mean")},
		{"OUT_SENSIBLE",   		VariableMetaData("W m-2", "SENSIBLE", "surface_upward_sensible_heat_flux", "Net upward sensible heat flux", "time: mean area: mean")},
		{"OUT_SNOW_FLUX",  		VariableMetaData("W m-2", "SNOW_FLUX", "downward_heat_flux_at_ground_level_in_snow", "Energy flux through snow pack", "time: mean area: mean")},
		//Miscellaneous terms
		{"OUT_AERO_COND",		VariableMetaData("m s-1", "AERO_COND", "", "scene aerodynamic conductance (tiles with overstory contribute overstory conductance; others contribute surface conductance)", "time: mean area: mean")},
		{"OUT_AERO_COND1",		VariableMetaData("m s-1", "AERO_COND1", "", "Surface aerodynamic conductance", "time: mean area: mean")},
		{"OUT_AERO_COND2",		VariableMetaData("m s-1", "AERO_COND2", "", "Overstory aerodynamic conductance", "time: mean area: mean")},
		{"OUT_AERO_RESIST",    		VariableMetaData("s m-1", "AERO_RESIST", "", "scene canopy aerodynamic resistance (tiles with overstory contribute overstory resistance; others contribute surface resistance)", "time: mean area: mean")},
		{"OUT_AERO_RESIST1",    	VariableMetaData("s m-1", "AERO_RESIST1", "", "Surface aerodynamic resistance", "time: mean area: mean")},
		{"OUT_AERO_RESIST2",    	VariableMetaData("s m-1", "AERO_RESIST2", "", "Overstory aerodynamic resistance", "time: mean area: mean")},
		{"OUT_AIR_TEMP",   		VariableMetaData("degree_Celsius", "AIR_TEMP", "air_temperature", "Air temperature", "time: mean area: mean")},
		{"OUT_DENSITY",    		VariableMetaData("kg m-3", "DENSITY", "air_density", "Near surface atmospheric density", "time: mean area: mean")},
		{"OUT_LONGWAVE",    		VariableMetaData("W m-2", "LONGWAVE", "downwelling_longwave_flux_in_air", "Incoming longwave", "time: mean area: mean")},
		{"OUT_PRESSURE",    		VariableMetaData("kPa", "PRESSURE", "surface_air_pressure", "Near surface atmospheric pressure", "time: mean area: mean")},
		{"OUT_QAIR",			VariableMetaData("1", "QAIR", "surface_specific_humidity", "Specific humidity", "time: mean area: mean")},
		{"OUT_REL_HUMID",  		VariableMetaData("1", "REL_HUMID", "relative_humidity", "Relative humidity", "time: mean area: mean")},
		{"OUT_SHORTWAVE",    		VariableMetaData("W m-2", "SHORTWAVE", "downwelling_shortwave_flux_in_air", "Incoming shortwave", "time: mean area: mean")},
		{"OUT_SURF_COND",		VariableMetaData("m s-1", "SURF_COND", "", "Surface conductance", "time: mean area: mean")},
		{"OUT_TSKC",    		VariableMetaData("1", "TSKC", "cloud_area_fraction", "Cloud cover fraction", "time: mean area: mean")},
		{"OUT_VP",    			VariableMetaData("kPa", "VP", "water_vapor_partial_pressure_in_air", "Near surface vapor pressure", "time: mean area: mean")},
		{"OUT_VPD",    			VariableMetaData("kPa", "VPD", "water_vapor_saturation_deficit_in_air", "Near surface vapor pressure deficit", "time: mean area: mean")},
		{"OUT_WIND",       		VariableMetaData("m s-1", "WIND", "wind_speed", "Near surface wind speed", "time: mean area: mean")},
		{"OUT_AREA",       		VariableMetaData("m2", "AREA", "area", "Grid cell area", "time: point area: sum")},
		//Band-specific terms
		{"OUT_ADV_SENS_BAND", 		VariableMetaData("W m-2", "ADV_SENS_BAND", "", "Net sensible flux advected to snow pack", "time: mean area: mean")},
		{"OUT_ADVECTION_BAND",     	VariableMetaData("W m-2", "ADVECTION_BAND", "", "advected energy", "time: mean area: mean")},
		{"OUT_ALBEDO_BAND",    		VariableMetaData("1", "ALBEDO_BAND", "surface_albedo", "Average surface albedo", "time: point area: mean")},
		{"OUT_AREA_BAND",		VariableMetaData("1", "AREA_BAND", "", "Band area", "time: point area: sum")},
		{"OUT_DELTACC_BAND",     	VariableMetaData("W m-2", "DELTACC_BAND", "", "Change in cold content in snow pack", "time: mean area: mean")},
		{"OUT_ELEV_BAND",     		VariableMetaData("m", "ELEV_BAND", "", "Band elevation", "time: point area: median")},
		{"OUT_GRND_FLUX_BAND",     	VariableMetaData("W m-2", "GRND_FLUX_BAND", "downward_heat_flux_at_ground_level_in_soil", "Net heat flux into ground", "time: mean area: mean")},
		{"OUT_IN_LONG_BAND",    	VariableMetaData("W m-2", "IN_LONG_BAND", "", "Incoming longwave at ground surface (under vegetation) ", "time: mean area: mean")},
		{"OUT_LATENT_BAND",    		VariableMetaData("W m-2", "LATENT_BAND", "surface_upward_latent_heat_flux", "net upward latent heat flux", "time: mean area: mean")},
		{"OUT_LATENT_SUB_BAND",  	VariableMetaData("W m-2", "LATENT_SUB_BAND", "", "Net upward latent heat flux from sublimation", "time: mean area: mean")},
		{"OUT_MELT_ENERGY_BAND", 	VariableMetaData("W m-2", "MELT_ENERGY_BAND", "surface_snow_melt_heat_flux", "Energy of fusion (melting) in snowpack", "time: mean area: mean")},
		{"OUT_NET_LONG_BAND",    	VariableMetaData("W m-2", "NET_LONG_BAND", "net_downward_longwave_flux_in_air", "Net downward longwave flux", "time: mean area: mean")},
		{"OUT_NET_SHORT_BAND", 		VariableMetaData("W m-2", "NET_SHORT_BAND", "net_downward_shortwave_flux_in_air", "Net downward shortwave flux", "time: mean area: mean")},
		{"OUT_RFRZ_ENERGY_BAND",	VariableMetaData("W m-2", "RFRZ_ENERGY_BAND", "", "Net energy used to refreeze liquid water in snowpack", "")},
		{"OUT_SENSIBLE_BAND",    	VariableMetaData("W m-2", "SENSIBLE_BAND", "surface_upward_sensible_heat_flux", "Net upward sensible heat flux", "time: mean area: mean")},
		{"OUT_SNOW_CANOPY_BAND",	VariableMetaData("mm", "SNOW_CANOPY_BAND", "", "Snow interception storage in canopy", "time: point area: mean")},
		{"OUT_SNOW_COVER_BAND",		VariableMetaData("1", "SNOW_COVER_BAND", "surface_snow_area_fraction", "Snow area fraction", "time: point area: sum")},
		{"OUT_SNOW_DEPTH_BAND",		VariableMetaData("cm", "SNOW_DEPTH_BAND", "surface_snow_thickness", "Snow Depth", "time: point area: mean")},
		{"OUT_SNOW_FLUX_BAND",		VariableMetaData("W m-2", "SNOW_FLUX_BAND", "downward_heat_flux_at_ground_level_in_snow", "Energy flux through snow pack", "time: mean area: mean")},
		{"OUT_SNOW_MELT_BAND",		VariableMetaData("mm", "SNOW_MELT_BAND", "thickness_of_surface_snow_melt_amount", "Snow melt", "time: mean area: mean")},
		{"OUT_SNOW_PACKT_BAND",		VariableMetaData("degree_Celsius", "SNOW_PACKT_BAND", "", "Snow pack temperature", "time: point area: mean")},
		{"OUT_SNOW_SURFT_BAND",		VariableMetaData("degree_Celsius", "SNOW_SURFT_BAND", "surface_temperature_where_snow", "Snow surface temperature", "time: point area: mean")},
		{"OUT_SWE_BAND",		VariableMetaData("mm", "SWE_BAND", "liquid_water_content_of_snow_layer", "Snow water equivalent", "time: mean area: mean")},
		//Glacier water balance terms - state variables
		{"OUT_GLAC_AREA",    		VariableMetaData("1", "GLAC_AREA", "", "Glacier surface area fraction", "time: point area: sum")},
		{"OUT_GLAC_WAT_STOR",    	VariableMetaData("mm", "GLAC_WAT_STOR", "", "glacier water storage", "time: point area: mean")},
		//Glacier water balance terms - fluxes
		{"OUT_GLAC_ACCUM",     		VariableMetaData("mm", "GLAC_ACCUM", "", "Glacier ice accumulation from conversion of firn to ice", "time: mean area: mean")},
		{"OUT_GLAC_IMBAL",     		VariableMetaData("mm", "GLAC_IMBAL", "", "Glacier ice mass balance", "time: mean area: mean")},
		{"OUT_GLAC_INFLOW",    		VariableMetaData("mm", "GLAC_INFLOW", "", "Glacier water inflow from snow melt, ice melt and rainfall", "time: mean area: mean")},
		{"OUT_GLAC_MBAL",   		VariableMetaData("mm", "GLAC_MBAL", "", "Glacier mass balance", "time: mean area: mean")},
		{"OUT_GLAC_MELT",    		VariableMetaData("mm", "GLAC_MELT", "", "Glacier ice melt", "time: mean area: mean")},
		{"OUT_GLAC_OUTFLOW",     	VariableMetaData("mm", "GLAC_OUTFLOW", "", "Glacier water outflow", "time: mean area: mean")},
		{"OUT_GLAC_OUTFLOW_COEF", 	VariableMetaData("1", "GLAC_OUTFLOW_COEF", "", "Glacier outflow coefficient", "time: mean area: mean")},
		{"OUT_GLAC_SUB",     		VariableMetaData("mm", "GLAC_SUB", "", "Net sublimation of glacier ice", "time: mean area: mean")},
		//Glacier energy balance terms - state variables
		{"OUT_GLAC_SURF_TEMP",     	VariableMetaData("degree_Celsius", "GLAC_SURF_TEMP", "land_ice_temperature", "Glacier surface temperature", "time: point area: mean")},
		{"OUT_GLAC_TSURF_FBFLAG",   VariableMetaData("", "GLAC_TSURF_FBFLAG", "", "Glacier surface temperature fallback count", "time: point area: mean")},
		//Glacier energy balance terms - fluxes
		{"OUT_GLAC_DELTACC",     	VariableMetaData("W m-2", "GLAC_DELTACC", "", "rate of change of cold content in glacier surface layer", "time: mean area: mean")},
		{"OUT_GLAC_FLUX",    		VariableMetaData("W m-2", "GLAC_FLUX", "", "energy flux through glacier surface layer", "time: mean area: mean")},
		{"OUT_GLAC_MELT_ENERGY",	VariableMetaData("W m-2", "GLAC_MELT_ENERGY", "", "energy used to thaw glacier ice", "time: mean area: mean")},
		//Glacier band-specific quantities
		{"OUT_GLAC_AREA_BAND",     	VariableMetaData("1", "GLAC_AREA_BAND", "", "Glacier surface area fraction", "time: point area: sum")},
		{"OUT_GLAC_WAT_STOR_BAND",	VariableMetaData("mm", "GLAC_WAT_STOR_BAND", "", "Glacier water storage", "")},
		{"OUT_GLAC_ACCUM_BAND",    	VariableMetaData("mm", "GLAC_ACCUM_BAND", "", "Glacier ice accumulation from conversion of firn to ice", "time: mean area: mean")},
		{"OUT_GLAC_IMBAL_BAND",   	VariableMetaData("mm", "GLAC_IMBAL_BAND", "", "Glacier ice mass balance", "time: mean area: mean")},
		{"OUT_GLAC_INFLOW_BAND",   	VariableMetaData("mm", "GLAC_INFLOW_BAND", "", "Glacier water inflow from snow melt, ice melt and rainfall", "time: mean area: mean")},
		{"OUT_GLAC_MBAL_BAND",     	VariableMetaData("mm", "GLAC_MBAL_BAND", "", "Glacier mass balance", "time: mean area: mean")},
		{"OUT_GLAC_MELT_BAND",     	VariableMetaData("mm", "GLAC_MELT_BAND", "", "Glacier ice melt", "time: mean area: mean")},
		{"OUT_GLAC_OUTFLOW_BAND",	VariableMetaData("mm", "GLAC_OUTFLOW_BAND", "", "Glacier water outflow", "time: mean area: mean")},
		{"OUT_GLAC_SUB_BAND",    	VariableMetaData("mm", "GLAC_SUB_BAND", "", "Net sublimation of glacier ice", "time: mean area: mean")},
		{"OUT_GLAC_DELTACC_BAND", 	VariableMetaData("W m-2", "GLAC_DELTACC_BAND", "", "Rate of change of cold content in glacier surface layer", "time: mean area: mean")},
		{"OUT_GLAC_FLUX_BAND",     	VariableMetaData("W m-2", "GLAC_FLUX_BAND", "", "Energy flux through glacier surface layer", "time: mean area: mean")}
	};
}

// Changes value of name member for a given output_mapping variable key.  VIC will name the variable accordingly in the output file.
void ProgramState:: set_output_variable_name(std::string variableKey, std::string newName) {

	if (output_mapping.find(variableKey) == output_mapping.end()) {
	        throw VICException("Error: set_output_variable_name could not find variable in output_mapping: " + variableKey);
	}
	output_mapping.at(variableKey).name = newName;
}



