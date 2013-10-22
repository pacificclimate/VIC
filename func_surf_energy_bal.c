#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vicNl.h"
#include "surf_energy_bal.h"

static char vcid[] = "$Id$";

double SurfEnergyBal::calculate(double Ts)
/**********************************************************************
	func_surf_energy_bal	Keith Cherkauer		January 3, 1996

  This subroutine computes the surface energy balance for bare soil
  and vegetation uncovered by snow.  It computes outgoing longwave,
  sensible heat flux, ground heat flux, and storage of heat in the thin
  upper layer, based on the given surface temperature.

  The Energy Balance Equation used comes from Xu Liang's Paper 
  "Insights of the Ground Heat Flux in Land Surface Parameterization
  Schemes."

  Modifications:
  04-14-98 modified to compute evapotranspiration within this routine
           in the hopes of reducing the number of iteration 
	  needed to find a solution surface temperature.       KAC
  07-13-98 modified to include elevation bands for vegetation 
           and snow                                             KAC
  01-20-00 modified to work with the updated radiation estimation
           routines, as well as the simplified frozen soil moisture
           storage                                              KAC
  01-08-01 sensible heat is now set to 0 W/m^2 when the ground is
           fully covered by snow.                                KAC
  04-12-01 fixed error where sensible heat flux for partial bare
           ground was not multiplied by the snow cover fraction. KAC
  04-29-02 moved calculation of sensible heat so that it is computed
           even in water balance mode.  This assures that it is set
           to 0 in water balance mode and does not yield the 
           cumulative sum of sensible heat from the snowpack.    KAC
  11-18-02 modified to compute the effects of blowing snow on the
           surface energy balance.                               LCB
  10-May-04 Added check that both FS_ACTIVE and FROZEN_SOIL are true
	    before computing *fusion.  This is just a safety measure;
	    ice and ice0 should both be 0 if FS_ACTIVE is FALSE.TJB
  16-Jul-04 Renamed VaporMassFlux, BlowingMassFlux, and SurfaceMassFlux
	    to vapor_flux, blowing_flux, and surface_flux, respectively,
	    to denote fact that their units are m/timestep rather than
	    kg/m2s.  Created new variables VaporMassFlux, BlowingMassFlux,
	    and SurfaceMassFlux with units of kg/m2s.  The addresses of
	    the *MassFlux variables are passed to latent_heat_from_snow()
	    where values for the variables are computed.  After these
	    values are computed, vapor_flux, blowing_flux and surface_flux
	    are derived from them by unit conversion.  vapor_flux,
	    blowing_flux, and surface_flux are the variables that are
	    passed in/out of this function.			TJB
  16-Jul-04 Changed the type of the last few variables (lag_one, Nveg,
	    etc) in the va_list to be double.  For some reason, passing
	    them as float or int caused them to become garbage.  This may
	    have to do with the fact that they followed variables of type
	    (double *) in va_list, which may have caused memory alignment
	    problems.						TJB
  05-Aug-04 Removed iveg, LastSnow, dt, SnowDepth, lag_one, sigma_slope,
	    fetch, and Nveg from the this function's argument list,
	    since these variables were only used in the call to
	    latent_heat_from_snow() which no longer needs them.	TJB
  28-Sep-04 Added Ra_used to store the aerodynamic resistance used in
	    flux calculations.					TJB
  2007-Apr-11 Modified to handle grid cell errors by returning to the
              main subroutine, rather than ending the simulation.	GCT
  2007-Apr-24 Removed (1.-snow_coverage) from three equations where it did not
              belong: for calculating LongBareOut in the second two cases and for
              calculating NetBareRad in the third case.			JCA
  2007-Apr-24 Features included for IMPLICIT frozen soils option.	JCA
              (including passing in nrec, nrecs, and iveg)
              (including passing in bulk_density, soil_density, and quartz)
              (including counting cases when IMPLICIT fails and involes EXPLICIT)
  2007-Apr-24 Features included for EXP_TRANS frozen soils option.	JCA
  2007-Apr-24 Passing in Zsum_node.					JCA
  2007-Aug-07 Moved Implicit error counting above call for
              solve_T_profile.						JCA
  2008-Aug-08 Added option for EXCESS_ICE.				JCA
              including: passing entire soil_con structure to
              calc_surf_energy_bal
  2007-Aug-24 Modified to use arno_evap rather than canopy_evap if LAI
              is 0, e.g. winter cropland.				KAC via TJB
  2008-Mar-01 Fixed typo in declaration of ufwc_table_layer.		TJB
  2009-Feb-09 Modified to remove dz_node.				KAC via TJB
  2009-May-20 Corrected the deltaH and fusion terms to account for
	      surf_atten, as in Liang et al, 1999.
	      Added options.GRND_FLUX_TYPE to allow backwards-compatibility
	      with versions 4.0.6 and 4.1.0.				TJB
  2009-Jun-19 Added T flag to indicate whether TFALLBACK occurred.	TJB
  2009-Jul-26 Removed the special logic for the water balance mode, in
	      which net longwave is stored in the "longwave" variable.	TJB
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.		TJB
  2009-Nov-15 Changed definitions of D1 and D2 to work for arbitrary
	      node spacing.						TJB
  2010-Apr-24 Replaced ra_under with Ra_used[0].			TJB
  2010-Apr-28 Removed net_short, displacement, roughness, and ref_height
	      from arg list of arno_evap() as they are no longer used.	TJB
  2011-May-31 Removed options.GRND_FLUX.  Soil temperatures and ground
	      flux are now always computed.				TJB
  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.				TJB
  2011-Aug-09 Now method used for estimating soil temperatures depends only
	      on QUICK_FLUX setting.					TJB

**********************************************************************/
{
  /* define routine input variables */

  /* general model terms */
  int Error;

  /* Define internal routine variables */
  double Evap;		/** Total evap in m/s **/
  double LongBareOut; // outgoing LW from snow-free ground
  double NetBareRad;
  double TMean;
  double Tmp;
  double error;
  double ice;
  double temp_latent_heat;
  double temp_latent_heat_sub;
  double VaporMassFlux;
  double BlowingMassFlux;
  double SurfaceMassFlux;

#if EXCESS_ICE
  porosity = soil_con->porosity[0];
  effective_porosity = soil_con->effective_porosity[0];
  porosity_node = soil_con->porosity_node;
  effective_porosity_node = soil_con->effective_porosity_node;
#endif

  /***************
    MAIN ROUTINE
  ***************/

  Error = 0;

  TMean = Ts;
  Tmp = TMean + KELVIN;

  /**********************************************
    Compute Surface Temperature at Half Time Step
  **********************************************/

  if ( snow_coverage > 0 && !INCLUDE_SNOW ) {

    /****************************************
      Compute energy flux through snow pack
    ****************************************/

    *snow_flux = ( kappa_snow * (Tsnow_surf - TMean) );

  }
  else if ( INCLUDE_SNOW ) {
    *snow_flux = 0;
    Tsnow_surf = TMean;
  }
  else *snow_flux = 0;

  /***************************************************************
    Estimate soil temperatures for ground heat flux calculations
  ***************************************************************/

  if ( state->options.QUICK_FLUX ) {
    /**************************************************************
      Use Liang et al. 1999 Equations to Calculate Ground Heat Flux
    **************************************************************/
    *T1 = estimate_T1(TMean, T1_old, T2, D1, D2, kappa1, kappa2, Cs1, Cs2, dp, delta_t);
    
    /*****************************************************
      Compute the Ground Heat Flux from the Top Soil Layer
    *****************************************************/
    if (state->options.GRND_FLUX_TYPE == GF_406) {
      *grnd_flux = (snow_coverage + (1. - snow_coverage) * surf_atten) * (kappa1 / D1 * ((*T1) - TMean));
    }
    else {
      *grnd_flux = (snow_coverage + (1. - snow_coverage) * surf_atten) * (kappa1 / D1 * ((*T1) - TMean)
				       + (kappa2 / D2 * ( 1. - exp( -D1 / dp )) * (T2 - (*T1)))) / 2.;
    }
      
  }
  else {
    /*************************************************************
      Use Finite Difference Method to Solve Ground Heat
      Flux at Soil Thermal Nodes (Cherkauer and Lettenmaier, 1999)
    *************************************************************/
    T_node[0] = TMean;
      
    /* IMPLICIT Solution */
    if(state->options.IMPLICIT) {
      Error = solve_T_profile_implicit(Tnew_node, T_node, kappa_node, Cs_node,
          moist_node, delta_t, ice_node, dp, Nnodes, FIRST_SOLN, NOFLUX,
          EXP_TRANS, veg_class, soil_con, state);
      
      /* print out error information for IMPLICIT solution */
      if(Error==0)
        error_cnt0++;
      else
        error_cnt1++;
      if(FIRST_SOLN[1]){
        FIRST_SOLN[1] = FALSE;
#if VERBOSE
        if ( iveg == 0 && rec == nrecs - 1) 
          fprintf(stderr,"The implicit scheme failed %d instances (%.1f%c of attempts).\n",error_cnt1,100.0*(float)error_cnt1/((float)error_cnt0+(float)error_cnt1),'%');
#endif
      }
    }

    /* EXPLICIT Solution, or if IMPLICIT Solution Failed */
    if(!state->options.IMPLICIT || Error == 1) {
      if(state->options.IMPLICIT)
        FIRST_SOLN[0] = TRUE;


      Error = solve_T_profile(Tnew_node, T_node, Tnew_fbflag, Tnew_fbcount,
          kappa_node, Cs_node, moist_node, delta_t, ice_node, dp,
          soil_con->ufwc_table_node,Nnodes,
          FIRST_SOLN, NOFLUX, EXP_TRANS, veg_class, soil_con, state);
    }
      
    if ( (int)Error == ERROR ) {
      fprintf(stderr, "ERROR: func_surf_energy_bal calling solve_T_profile\n");
      return( ERROR ); 
    }
    *T1 = Tnew_node[1];
      

    /*****************************************************
      Compute the Ground Heat Flux between nodes 0 and 1
    *****************************************************/
    if (state->options.GRND_FLUX_TYPE == GF_406) {
      *grnd_flux = (snow_coverage + (1. - snow_coverage) * surf_atten) * (kappa1 / D1 * ((*T1) - TMean));
    }
    else {
      *grnd_flux = (snow_coverage + (1. - snow_coverage) * surf_atten) * (kappa1 / D1 * ((*T1) - TMean) 
		     + (kappa2 / D2 * (Tnew_node[2] - (*T1)))) / 2.;
    }
      
  }

  /******************************************************
    Compute the change in heat storage in the region between nodes 0 and 1
    (this will correspond to top soil layer for the default (non-exponential) node spacing)
  ******************************************************/
  if (state->options.GRND_FLUX_TYPE == GF_FULL) {
    *deltaH = (snow_coverage + (1. - snow_coverage) * surf_atten)
              * (Cs1 * ((Ts_old + T1_old) - (TMean + *T1)) * D1 / delta_t / 2.);
  }
  else {
    *deltaH = (Cs1 * ((Ts_old + T1_old) - (TMean + *T1)) * D1 / delta_t / 2.);
  }

  /******************************************************
    Compute the change in heat due to solid-liquid phase changes in the region between nodes 0 and 1
    (this will correspond to top soil layer for the default (non-exponential) node spacing)
  ******************************************************/
  if (soil_con->FS_ACTIVE && state->options.FROZEN_SOIL) {
    if((TMean+ *T1)/2.<0.) {
      ice = moist - maximum_unfrozen_water((TMean+ *T1)/2.,
					   max_moist,bubble,expt);
      if(ice<0.) ice=0.;
    }
    else ice=0.;

    if (state->options.GRND_FLUX_TYPE == GF_FULL) {
      *fusion = (snow_coverage + (1. - snow_coverage) * surf_atten)
                * (-ice_density * Lf * (ice0 - ice) * D1 / delta_t);
    }
    else {
      *fusion = (-ice_density * Lf * (ice0 - ice) * D1 / delta_t);
    }

  }
    
  /* if thin snowpack, compute the change in energy stored in the pack */
  if ( INCLUDE_SNOW ) {
    if ( TMean > 0 )
      *deltaCC = CH_ICE * (snow_swq - snow_water) * (0 - OldTSurf) / delta_t;
    else
      *deltaCC = CH_ICE * (snow_swq - snow_water) * (TMean - OldTSurf) / delta_t;
    *refreeze_energy = (snow_water * Lf * snow_density) / delta_t;
    *deltaCC *= snow_coverage; // adjust for snow cover fraction
    *refreeze_energy *= snow_coverage; // adjust for snow cover fraction
  }
    
  /** Compute net surface radiation of snow-free area for evaporation estimates **/
  LongBareOut = STEFAN_B * Tmp * Tmp * Tmp * Tmp;
  if ( INCLUDE_SNOW ) { // compute net LW at snow surface
    (*NetLongSnow) = (LongSnowIn - snow_coverage * LongBareOut);
  }
  (*NetLongBare)   = (LongBareIn - (1. - snow_coverage) * LongBareOut); // net LW snow-free area
  NetBareRad = (NetShortBare + (*NetLongBare) + *grnd_flux + *deltaH + *fusion);
    
  /** Compute atmospheric stability correction **/
  /** CHECK THAT THIS WORKS FOR ALL SUMMER SITUATIONS **/
  if ( wind[UnderStory] > 0.0 && overstory && SNOWING )
    Ra_used[0] = ra[UnderStory] 
      / StabilityCorrection(ref_height[UnderStory], 0.f, TMean, Tair, 
			    wind[UnderStory], roughness[UnderStory]);
  else if ( wind[UnderStory] > 0.0 )
    Ra_used[0] = ra[UnderStory] 
      / StabilityCorrection(ref_height[UnderStory], displacement[UnderStory], 
			    TMean, Tair, wind[UnderStory], roughness[UnderStory]);
  else
    Ra_used[0] = HUGE_RESIST;
  
  /*************************************************
    Compute Evapotranspiration if not snow covered

    Should evapotranspiration be active when the 
    ground is only partially covered with snow????

    Use Arno Evap if LAI is set to zero (e.g. no
    winter crop planted).
  *************************************************/
  if ( VEG && !SNOWING && state->veg_lib[veg_class].LAI[month-1] > 0 ) {
    Evap = canopy_evap(layer_wet, layer_dry, veg_var_wet, veg_var_dry, TRUE, 
		       veg_class, month, precipitation_mu, Wdew, delta_t, NetBareRad, vpd, 
		       NetShortBare, Tair, Ra_used[1], 
		       displacement[1], roughness[1], ref_height[1], 
		       (double)soil_con->elevation, rainfall, soil_con->depth, soil_con->Wcr, soil_con->Wpwp,
		       soil_con->frost_fract, root, state);
  }
  else if(!SNOWING) {
    Evap = arno_evap(layer_wet, layer_dry, NetBareRad, Tair, vpd, 
		     soil_con->depth[0], max_moist * soil_con->depth[0] * 1000., 
		     (double)soil_con->elevation, soil_con->b_infilt, Ra_used[0], delta_t, precipitation_mu,
		     soil_con->resid_moist[0], soil_con->frost_fract, state);

  }
  else Evap = 0.;
  
  /**********************************************************************
    Compute the Latent Heat Flux from the Surface and Covering Vegetation
  **********************************************************************/
  *latent_heat  = -RHO_W * latent_heat_Le * Evap;
  *latent_heat_sub = 0.;

  /** Compute the latent heat flux from a thin snowpack if present **/
  if (INCLUDE_SNOW) {

    /* Convert sublimation terms from m/timestep to kg/m2s */
    VaporMassFlux = *vapor_flux * ice_density / delta_t;
    BlowingMassFlux = *blowing_flux * ice_density / delta_t;
    SurfaceMassFlux = *surface_flux * ice_density / delta_t;

    latent_heat_from_snow(atmos_density, vp, latent_heat_Le, atmos_pressure, 
			  Ra_used[0], TMean, vpd, &temp_latent_heat, 
			  &temp_latent_heat_sub, &VaporMassFlux,
                          &BlowingMassFlux, &SurfaceMassFlux);
    *latent_heat += temp_latent_heat * snow_coverage;
    *latent_heat_sub = temp_latent_heat_sub * snow_coverage;

    /* Convert sublimation terms from kg/m2s to m/timestep */
    *vapor_flux = VaporMassFlux * delta_t / ice_density;
    *blowing_flux = BlowingMassFlux * delta_t / ice_density;
    *surface_flux = SurfaceMassFlux * delta_t / ice_density;

  }
  else *latent_heat *= (1. - snow_coverage);

  /************************************************
    Compute the Sensible Heat Flux from the Surface
  ************************************************/
  if ( snow_coverage < 1 || INCLUDE_SNOW ) {
    *sensible_heat = atmos_density * Cp * (Tair - (TMean)) / Ra_used[0];
    if ( !INCLUDE_SNOW ) (*sensible_heat) *= (1. - snow_coverage);
  }
  else *sensible_heat = 0.;

  /*************************************
    Compute Surface Energy Balance Error
  *************************************/
  error = (NetBareRad // net radiation on snow-free area
	   + NetShortGrnd + NetShortSnow // net surface SW 
	   + emissivity * (*NetLongSnow)) // net surface LW 
	   + *sensible_heat // surface sensible heat
	   + (*latent_heat + *latent_heat_sub) // surface latent heats
	   /* heat flux through snowpack - for snow covered fraction */
	   + *snow_flux * snow_coverage
	   /* energy used in reducing snow coverage area */
	   + melt_energy 
	   /* snow energy terms - values are 0 unless INCLUDE_SNOW */
	   + Advection - *deltaCC;
    
  if ( INCLUDE_SNOW ) {
    if (Tsnow_surf == 0.0 && error > -(*refreeze_energy)) {
      *refreeze_energy = -error;
      error = 0.0;
    }
    else {
      error += *refreeze_energy;
    }
  }

  *store_error = error;

  return error;

}

