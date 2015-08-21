#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

static char vcid[] = "$Id$";

int  put_data(cell_info_struct  *cell,
              WriteOutputFormat *output,
              out_data_struct   *out_data,
              const dmy_struct  *dmy,
              int                rec,
              const ProgramState* state)
/**********************************************************************
	put_data.c	Dag Lohmann		January 1996

  This routine converts data units, and stores finalized values
  in an array for later output to the output files.

  modifications:
  06-24-98  modified for new distributed precipitation data structures KAC
  01-20-00 modified to deal with simplified frozen soil moisture layers
           and frost depth / thaw depth accounting                 KAC
  03-08-00 modified to eliminate extra lines for storing bare
           soil variables.                                         KAC
  6-8-2000 modified to handle spatially distribute frozen soil     KAC
  10-6-2000 modified to handle partial snow cover                  KAC
  02-27-01 modified to output lake model variables                 KAC
  11-18-02 updated output of lake variables to reflect algorithm
           changes.  Also added output variables for blowing snow
           algorithm.                                              LCB
  03-12-03 modified to add additional energy balance variable storage
           when output of snow bands is selected.                  KAC
  03-12-03 Modified to add AboveTreeLine to soil_con_struct so that
           the model can make use of the computed treeline.     KAC
  30-Oct-03 Snow_flux was incorrectly set to Tcanopy.  Fixed.   TJB
  25-Aug-04 Sub_snow was incorrectly set to blowing_flux.  Now it is
            set to vapor_flux.                                  TJB
  28-Sep-04 Now out_data->aero_resist stores the aerodynamic resistance
            used in flux calculations.                          TJB
  2005-Mar-24 Modified to compute ALMA output variables.        TJB
  2005-Apr-23 Now aero_cond is aggregated instead of aero_resist.       TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
              out_data and out_data_files structures; removed the
              OPTIMIZE and LDAS_OUTPUT options; uses the
	      new save_data structure; implemented aggregation.		TJB
  2006-Oct-10 Shortened the names of output variables whose names were
	      too long; fixed typos in others; created new OUT_IN_LONG
	      variable.							TJB
  2006-Nov-07 Added OUT_SOIL_TNODE.					TJB
  2006-Nov-07 Assigned value to overstory.				TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2006-Nov-30 Added OUT_DELSURFSTOR.					TJB
  2006-Nov-30 Convert pressure and vapor pressure to kPa for output.	TJB
  2006-Dec-20 Changed OUT_SURF_TEMP from average of T[0] and T[1] to
	      direct assignment of T[0].				TJB
  2007-Apr-21 Moved initialization of tmp_fract to immediately before the
	        #if SPATIAL_FROST
	      block, so that it would be initialized in all cases.	TJB
  2007-Aug-17 Added EXCESS_ICE output variables.			JCA
  2007-Aug-22 Added OUT_WATER_ERROR as output variable.			JCA
  2007-Nov-06 Lake area is now the larger of lake.areai and lake.sarea.
	      Added wetland canopyevap and canopy_vapor_flux to grid
	      cell flux aggregation.					LCB via TJB
  2008-Apr-21 Made computation of out_data[OUT_SURFSTOR] more robust.	TJB
  2008-Sep-09 Calculate sarea in order to output lake surface area at
	      the end of the time step.  The stored variable
	      lake->sarea represents the sarea from the beginning of
	      the time step, not the updated value from the end of the
	      time step.						LCB via TJB
  2008-Sep-09 Added SOIL_TNODE_WL as an output variable, the soil
	      temperature in the wetland fraction of the grid cell.	LCB via TJB
  2008-Sep-09 Allow output of wetland frost/thaw depths even if Clake
	      is 1.0 since wetland energy balance is always computed.	LCB via TJB
  2008-Sep-09 Lake depth assignment moved up to precede sarea
	      assignment.						LCB via TJB
  2008-Sep-09 Check to make sure that area > 0.0 when checking to see
	      if ice area > sarea.					LCB via TJB
  2008-Oct-23 Changed data type of put_data() to be int so that it
	      can return ErrorFlag.					TJB
  2009-Jan-12 Added a final return of (0) since the data type of put_data()
	      is int rather than void.					TJB
  2009-Jan-16 Modified aero_resist_used and Ra_used to become arrays of
	      two elements (surface and overstory); added
	      options.AERO_RESIST_CANSNOW.				TJB
  2009-Jan-16 Added AERO_COND1&2 and AERO_RESIST1&2 to track
	      surface and overstory values; changed AERO_COND
	      and AERO_RESIST to track "scene" values.			TJB
  2009-Feb-09 Removed checks on PRT_SNOW_BAND option.			TJB
  2009-Feb-22 Added OUT_VPD.						TJB
  2009-May-17 Added OUT_ASAT.						TJB
  2009-Jun-09 Modified to use extension of veg_lib structure to contain
	      bare soil information.					TJB
  2009-Jun-09 Added OUT_PET_*, potential evap computed for various
	      reference land cover types.				TJB
  2009-Jun-09 Cell_data structure now only stores final aero_resist
	      values (called "aero_resist").  Preliminary uncorrected
	      aerodynamic resistances for current vegetation and various
	      reference land cover types for use in potential evap
	      calculations is stored in temporary array aero_resist.	TJB
  2009-Jun-19 Added T flag to indicate whether TFALLBACK occurred.	TJB
  2009-Jul-31 Modified so that wetland veg is now included in main loop
	      over veg tiles and aggregated the same way as all other
	      veg tiles.						TJB
  2009-Aug-28 OUT_LAKE_ICE_TEMP and OUT_LAKE_SURF_TEMP are [C].		TJB
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.		TJB
  2009-Sep-28 Created collect_wb_terms and collect_eb_terms to handle
	      summing of storages and fluxes from upland veg tiles,
	      wetland veg tile, and lake.  Added logic to handle an
	      initial (pre-simulation) call for purpose of initializing
	      water and energy balance checks.				TJB
  2009-Sep-30 Miscellaneous fixes for lake model.			TJB
  2009-Oct-05 Modifications for taking changes in lake area into account.	TJB
  2009-Oct-08 Extended T fallback scheme to snow and ice T.		TJB
  2009-Nov-09 Changed definition of sarea to include ice extent.	LCB via TJB
  2009-Nov-15 Redirected T fallback messages to stderr.			TJB
  2009-Dec-11 Added logic for initialization of save_data structure.	TJB
  2010-Feb-14 Added OUT_LAKE_AREA_FRAC.					TJB
  2010-Mar-31 Added OUT_RUNOFF_IN.					TJB
  2010-Sep-24 Renamed RUNOFF_IN and OUT_RUNOFF_IN to CHANNEL_IN and
	      OUT_LAKE_CHAN_IN, respectively.  Renamed OUT_EVAP_LAKE
	      to OUT_LAKE_EVAP.  Added other lake water balance terms
	      to set of output variables.  Added volumetric versions 
	      of these too.						TJB
  2010-Nov-02 Added OUT_LAKE_RCHRG, OUT_LAKE_RCHRG_V, OUT_RO_IN,
	      OUT_LAKE_RO_IN_V, OUT_LAKE_VAPFLX, and OUT_LAKE_VAPFLX_V.	TJB
  2010-Nov-02 Changed units of lake_var moisture fluxes to volume (m3).	TJB
  2010-Nov-11 Moved assignment of all OUT_LAKE* variables outside
	      collect_wb_terms().  Added lakefactor to collect_wb_terms()
	      arg list.							TJB
  2010-Nov-21 Added OUT_LAKE_DSTOR, OUT_LAKE_DSTOR_V, OUT_LAKE_DSWE,
	      OUT_LAKE_DSWE_V, OUT_LAKE_SWE, and OUT_LAKE_SWE_V.	TJB
  2010-Nov-26 Changed += to = in assignment of OUT_LAKE_* variables.	TJB
  2010-Dec-01 Added OUT_ZWT.						TJB
  2011-Mar-01 Added OUT_ZWT2, OUT_ZWT3, and OUT_ZWTL.			TJB
  2011-Mar-31 Added frost_fract to collect_wb_terms() arglist.		TJB
  2011-Nov-04 Added OUT_TSKC.						TJB
**********************************************************************/
{
  int                     Ndist;
  int                     Nbands;
  int                     overstory;
  bool                    HasVeg;
  bool                    IsWet;
  bool                    HasGlac;
  double                 *frost_fract;
  double                  frost_slope;
  double                  dp;
  int                     skipyear;
  double                  Cv;
  double                  Clake;
  double                  Cv_save;
  double                  precipitation_mu;
  double                  cv_baresoil;
  double                  cv_veg;
  double                  cv_overstory;
  double                  cv_snow;
  double                  cv_glacier;
  double                  inflow;
  double                  outflow;
  double                  glac_icebal;
  double                  storage;
  double                  TreeAdjustFactor[MAX_BANDS];
  double                  ThisAreaFract;
  double                  ThisTreeAdjust;
  int                     dt_sec;
  int                     out_dt_sec;
  int                     out_step_ratio;
  int                     ErrorFlag;


  frost_fract = cell->soil_con.frost_fract;
  frost_slope = cell->soil_con.frost_slope;

  dp = cell->soil_con.dp;
  skipyear = state->global_param.skipyear;
  dt_sec = state->global_param.dt*SECPHOUR;
  out_dt_sec = state->global_param.out_dt*SECPHOUR;
  out_step_ratio = (int)(out_dt_sec/dt_sec);
  if (rec >= 0) (cell->fallBackStats.step_count)++;

  if(state->options.DIST_PRCP)
    Ndist = 2;
  else 
    Ndist = 1;

  double bandCv [MAX_BANDS];
  // Initialize to zero.
  for (int band = 0; band < MAX_BANDS; band++) {
    bandCv[band] = 0;
  }

  for (std::vector<HRU>::iterator hru = cell->prcp.hruList.begin(); hru != cell->prcp.hruList.end(); ++hru) {
    if (state->veg_lib[hru->veg_con.vegIndex].overstory) {
      if (state->options.LAKES && hru->veg_con.LAKE) {
        // Fraction of tile that is flooded
        Clake = cell->prcp.lake_var.sarea / cell->lake_con.basin[0];
        bandCv[hru->bandIndex] += hru->veg_con.Cv * (1 - Clake);
      } else {
        bandCv[hru->bandIndex] += hru->veg_con.Cv;
      }
    }
  }

  // Compute treeline adjustment factors
  for (int band = 0; band < state->options.SNOW_BAND; band++) {
    if (cell->soil_con.AboveTreeLine[band]) {
      TreeAdjustFactor[band] = 1. / (1. - bandCv[band]);
    } else {
      TreeAdjustFactor[band] = 1.;
    }
    if (TreeAdjustFactor[band] != 1 && rec == 0) {
      fprintf(stderr, "WARNING: Tree adjust factor for band %i is equal to %f.\n", band, TreeAdjustFactor[band]);
    }
  }

  cv_baresoil = 0;
  cv_veg = 0;
  cv_overstory = 0;
  cv_snow = 0;
  cv_glacier = 0;

  // Initialize output data to zero
  zero_output_list(out_data);

  /* MPN */
  // Set output versions of input forcings
  if (rec >= 0) {
    out_data[OUT_AIR_TEMP].data[0] = cell->atmos[rec].air_temp[state->NR];
    out_data[OUT_DENSITY].data[0] = cell->atmos[rec].density[state->NR];
    out_data[OUT_LONGWAVE].data[0] = cell->atmos[rec].longwave[state->NR];
    out_data[OUT_PREC].data[0] = cell->atmos[rec].out_prec; // mm over grid cell
    out_data[OUT_PRESSURE].data[0] = cell->atmos[rec].pressure[state->NR]
        / kPa2Pa;
    out_data[OUT_QAIR].data[0] = EPS * cell->atmos[rec].vp[state->NR]
        / cell->atmos[rec].pressure[state->NR];
    out_data[OUT_RAINF].data[0] = cell->atmos[rec].out_rain; // mm over grid cell
    out_data[OUT_REL_HUMID].data[0] = 100. * cell->atmos[rec].vp[state->NR]
        / (cell->atmos[rec].vp[state->NR] + cell->atmos[rec].vpd[state->NR]);
    if (state->options.LAKES && cell->lake_con.Cl[0] > 0)
      out_data[OUT_LAKE_CHAN_IN].data[0] =
          cell->atmos[rec].channel_in[state->NR]; // mm over grid cell
    else
      out_data[OUT_LAKE_CHAN_IN].data[0] = 0;
    out_data[OUT_SHORTWAVE].data[0] = cell->atmos[rec].shortwave[state->NR];
    out_data[OUT_SNOWF].data[0] = cell->atmos[rec].out_snow; // mm over grid cell
    out_data[OUT_TSKC].data[0] = cell->atmos[rec].tskc[state->NR];
    out_data[OUT_VP].data[0] = cell->atmos[rec].vp[state->NR] / kPa2Pa;
    out_data[OUT_VPD].data[0] = cell->atmos[rec].vpd[state->NR] / kPa2Pa;
    out_data[OUT_WIND].data[0] = cell->atmos[rec].wind[state->NR];
  }
  /****************************************
   Store Output for all Vegetation Types (except lakes)
   ****************************************/
  for (std::vector<HRU>::iterator hru = cell->prcp.hruList.begin(); hru != cell->prcp.hruList.end(); ++hru) {

    Cv = hru->veg_con.Cv;
    Clake = 0;
    Nbands = state->options.SNOW_BAND;
    IsWet = false;

    if (hru->isArtificialBareSoil == true || hru->isGlacier)
      HasVeg = false;
    else
      HasVeg = true;

    HasGlac = hru->isGlacier;

    if ( Cv > 0) {

      // Check if this is lake/wetland tile
      if (state->options.LAKES && hru->veg_con.LAKE) {
        Clake = cell->prcp.lake_var.sarea/cell->lake_con.basin[0];
        Nbands = 1;
        IsWet = true;
      }

      overstory = state->veg_lib[hru->veg_con.vegIndex].overstory;

      /*********************************
        Store Output for all Bands 
      *********************************/
        int band = hru->bandIndex;
        ThisAreaFract = cell->soil_con.AreaFract[hru->bandIndex];
        ThisTreeAdjust = TreeAdjustFactor[band];
        if (IsWet) {
          ThisAreaFract = 1;
          ThisTreeAdjust = 1;
        }

        if(ThisAreaFract > 0. && ( hru->isArtificialBareSoil
           || ( !cell->soil_con.AboveTreeLine[band] || (cell->soil_con.AboveTreeLine[band] && !overstory)))) {

          /*******************************************************
            Store Output for Wet and Dry Fractions
          *******************************************************/
          for (int dist = 0; dist < Ndist; dist++ ) {
            if(dist==0) 
              precipitation_mu = hru->mu;
            else 
              precipitation_mu = 1. - hru->mu;

            /** compute running totals of various landcovers **/
            if (HasVeg)
              cv_veg += Cv * precipitation_mu * ThisTreeAdjust;
            else
              cv_baresoil += Cv * precipitation_mu * ThisTreeAdjust;
            if (overstory)
              cv_overstory += Cv * precipitation_mu * ThisTreeAdjust;
            if (hru->snow.swq > 0.0)
              cv_snow += Cv * precipitation_mu * ThisTreeAdjust;
            if (HasGlac) {
              cv_glacier += Cv * precipitation_mu * ThisTreeAdjust;
            }

	    /*********************************
              Record Water Balance Terms 
	    *********************************/
            collect_wb_terms(hru->cell[dist],
                             hru->veg_var[dist],
                             hru->snow,
                             hru->glacier,
                             cell->prcp.lake_var,
                             precipitation_mu,
                             Cv,
                             ThisTreeAdjust,
                             HasVeg,
                             HasGlac,
                             false,
                             (1-Clake),
                             overstory,
                             cell->soil_con.depth,
                             frost_fract,
                             out_data, state);

	  } // End wet/dry loop

	  /**********************************
	    Record Energy Balance Terms
	  **********************************/
          collect_eb_terms(hru->energy,
                           hru->snow,
                           hru->glacier,
                           hru->cell[WET],
                           &cell->fallBackStats,
                           Cv,
                           ThisAreaFract,
                           ThisTreeAdjust,
                           HasVeg,
                           HasGlac,
                           false,
                           (1-Clake),
                           overstory,
                           band,
                           cell->soil_con.depth,
                           cell->soil_con.dz_node,
                           frost_fract,
                           frost_slope,
                           out_data, state);

          // Store Wetland-Specific Variables

          if (IsWet) {
            // Wetland soil temperatures
            for(int i=0;i < state->options.Nnode; i++) {
              out_data[OUT_SOIL_TNODE_WL].data[i] = hru->energy.T[i];
            }
          }

          /**********************************
            Record Lake Variables
          **********************************/
          if (IsWet) {

            // Override some variables of soil under lake with those of wetland
            // This is for those variables whose lake values shouldn't be included
            // in grid cell average
            // Note: doing this for eb terms will lead to reporting of eb errors 
            // this should be fixed when we implement full thermal solution beneath lake
          for (int i = 0; i < MAX_FRONTS; i++) {
            cell->prcp.lake_var.energy.fdepth[i] = hru->energy.fdepth[i];
            cell->prcp.lake_var.energy.tdepth[i] = hru->energy.fdepth[i];
          }
          for (int i = 0; i < state->options.Nnode; i++) {
            cell->prcp.lake_var.energy.ice_content[i] = hru->energy.ice_content[i];
            cell->prcp.lake_var.energy.T[i] = hru->energy.T[i];
          }
          for (int i = 0; i < N_PET_TYPES; i++) {
            cell->prcp.lake_var.soil.pot_evap[i] = hru->cell[WET].pot_evap[i];
          }
          cell->prcp.lake_var.soil.rootmoist = hru->cell[WET].rootmoist;
          cell->prcp.lake_var.energy.deltaH = hru->energy.deltaH;
          cell->prcp.lake_var.energy.fusion = hru->energy.fusion;
          cell->prcp.lake_var.energy.grnd_flux = hru->energy.grnd_flux;

          glac_data_struct invalidGlacier;  // Placeholder with internal variables initialized to INVALID by default.
          /*********************************
           Record Water Balance Terms
           *********************************/
            collect_wb_terms(cell->prcp.lake_var.soil,
                             hru->veg_var[WET],
                             cell->prcp.lake_var.snow,
                             invalidGlacier,
                             cell->prcp.lake_var,
                             1.0,
                             Cv,
                             ThisTreeAdjust,
                             false,
                             false,
                             true,
                             Clake,
                             overstory,
                             cell->soil_con.depth,
                             frost_fract,
                             out_data, state);

          /**********************************
           Record Energy Balance Terms
           **********************************/
            collect_eb_terms(cell->prcp.lake_var.energy,
                             cell->prcp.lake_var.snow,
                             invalidGlacier,
                             cell->prcp.lake_var.soil,
                             &(cell->fallBackStats),
                             Cv,
                             ThisAreaFract,
                             ThisTreeAdjust,
                             false,
                             false,
                             true,
                             Clake,
                             overstory,
                             band,
                             cell->soil_con.depth,
                             cell->soil_con.dz_node,
                             frost_fract,
                             frost_slope,
                             out_data, state);

            // Store Lake-Specific Variables
            //TODO: move this lake stuff outside of the for loop since it doesn't depend on veg or band (for improved efficiency)
            // Lake ice
            if (cell->prcp.lake_var.new_ice_area > 0.0) {
              out_data[OUT_LAKE_ICE].data[0]   = (cell->prcp.lake_var.ice_water_eq/cell->prcp.lake_var.new_ice_area) * ice_density / RHO_W;
              out_data[OUT_LAKE_ICE_TEMP].data[0]   = cell->prcp.lake_var.tempi;
              out_data[OUT_LAKE_ICE_HEIGHT].data[0] = cell->prcp.lake_var.hice;
              out_data[OUT_LAKE_SWE].data[0] = cell->prcp.lake_var.swe/cell->prcp.lake_var.areai; // m over lake ice
              out_data[OUT_LAKE_SWE_V].data[0] = cell->prcp.lake_var.swe; // m3
            }
            else {
              out_data[OUT_LAKE_ICE].data[0]   = 0.0;
              out_data[OUT_LAKE_ICE_TEMP].data[0]   = 0.0;
              out_data[OUT_LAKE_ICE_HEIGHT].data[0]   = 0.0;
              out_data[OUT_LAKE_SWE].data[0]   = 0.0;
              out_data[OUT_LAKE_SWE_V].data[0]   = 0.0;
            }
            out_data[OUT_LAKE_DSWE_V].data[0] = cell->prcp.lake_var.swe - cell->prcp.lake_var.swe_save; // m3
            out_data[OUT_LAKE_DSWE].data[0] = (cell->prcp.lake_var.swe - cell->prcp.lake_var.swe_save)*1000/cell->soil_con.cell_area; // mm over gridcell

          // Lake dimensions
          out_data[OUT_LAKE_AREA_FRAC].data[0] = Cv * Clake;
          out_data[OUT_LAKE_DEPTH].data[0] = cell->prcp.lake_var.ldepth;
          out_data[OUT_LAKE_SURF_AREA].data[0] = cell->prcp.lake_var.sarea;
          if (out_data[OUT_LAKE_SURF_AREA].data[0] > 0)
            out_data[OUT_LAKE_ICE_FRACT].data[0] = cell->prcp.lake_var.new_ice_area / out_data[OUT_LAKE_SURF_AREA].data[0];
          else
            out_data[OUT_LAKE_ICE_FRACT].data[0] = 0.;
          out_data[OUT_LAKE_VOLUME].data[0] = cell->prcp.lake_var.volume;
          out_data[OUT_LAKE_DSTOR_V].data[0] = cell->prcp.lake_var.volume - cell->prcp.lake_var.volume_save;
          out_data[OUT_LAKE_DSTOR].data[0] = (cell->prcp.lake_var.volume - cell->prcp.lake_var.volume_save)*1000/cell->soil_con.cell_area; // mm over gridcell

            // Other lake characteristics
            out_data[OUT_LAKE_SURF_TEMP].data[0]  = cell->prcp.lake_var.temp[0];
            if (out_data[OUT_LAKE_SURF_AREA].data[0] > 0) {
              out_data[OUT_LAKE_MOIST].data[0]      = (cell->prcp.lake_var.volume / cell->soil_con.cell_area) * 1000.; // mm over gridcell
              out_data[OUT_SURFSTOR].data[0]        = (cell->prcp.lake_var.volume / cell->soil_con.cell_area) * 1000.; // same as OUT_LAKE_MOIST
            }
            else {
              out_data[OUT_LAKE_MOIST].data[0] = 0;
              out_data[OUT_SURFSTOR].data[0] = 0;
            }

            // Lake moisture fluxes
            out_data[OUT_LAKE_BF_IN_V].data[0] = cell->prcp.lake_var.baseflow_in; // m3
            out_data[OUT_LAKE_BF_OUT_V].data[0] = cell->prcp.lake_var.baseflow_out; // m3
            out_data[OUT_LAKE_CHAN_IN_V].data[0] = cell->prcp.lake_var.channel_in; // m3
            out_data[OUT_LAKE_CHAN_OUT_V].data[0] = cell->prcp.lake_var.runoff_out; // m3
            out_data[OUT_LAKE_EVAP_V].data[0] = cell->prcp.lake_var.evapw; // m3
            out_data[OUT_LAKE_PREC_V].data[0] = cell->prcp.lake_var.prec; // m3
            out_data[OUT_LAKE_RCHRG_V].data[0] = cell->prcp.lake_var.recharge; // m3
            out_data[OUT_LAKE_RO_IN_V].data[0] = cell->prcp.lake_var.runoff_in; // m3
            out_data[OUT_LAKE_VAPFLX_V].data[0] = cell->prcp.lake_var.vapor_flux; // m3
            out_data[OUT_LAKE_BF_IN].data[0] = cell->prcp.lake_var.baseflow_in*1000./cell->soil_con.cell_area; // mm over gridcell
            out_data[OUT_LAKE_BF_OUT].data[0] = cell->prcp.lake_var.baseflow_out*1000./cell->soil_con.cell_area; // mm over gridcell
            out_data[OUT_LAKE_CHAN_OUT].data[0] = cell->prcp.lake_var.runoff_out*1000./cell->soil_con.cell_area; // mm over gridcell
            out_data[OUT_LAKE_EVAP].data[0] = cell->prcp.lake_var.evapw*1000./cell->soil_con.cell_area; // mm over gridcell
            out_data[OUT_LAKE_RCHRG].data[0] = cell->prcp.lake_var.recharge*1000./cell->soil_con.cell_area; // mm over gridcell
            out_data[OUT_LAKE_RO_IN].data[0] = cell->prcp.lake_var.runoff_in*1000./cell->soil_con.cell_area; // mm over gridcell
            out_data[OUT_LAKE_VAPFLX].data[0] = cell->prcp.lake_var.vapor_flux*1000./cell->soil_con.cell_area; // mm over gridcell

          } // End if options.LAKES etc.

      } // End if ThisAreaFract etc.

    } // End if Cv > 0

  } // End loop over veg
 

  /*****************************************
    Aggregation of Dynamic Soil Properties      
   *****************************************/
#if EXCESS_ICE
  for(index=0;index<state->options.Nlayer;index++) {
    out_data[OUT_SOIL_DEPTH].data[index]  = cell->soil_con.depth[index];
    out_data[OUT_SUBSIDENCE].data[index]  = cell->soil_con.subsidence[index];
    out_data[OUT_POROSITY].data[index]    = cell->soil_con.effective_porosity[index];
  }  
  for(index=0;index<state->options.Nnode;index++)
    out_data[OUT_ZSUM_NODE].data[index]   = cell->soil_con.Zsum_node[index];
#endif // EXCESS_ICE

  /*****************************************
    Finish aggregation of special-case variables
   *****************************************/
  // Normalize quantities that aren't present over entire grid cell
  if (cv_baresoil > 0) {
    out_data[OUT_BARESOILT].data[0] /= cv_baresoil;
  }
  if (cv_veg > 0) {
    out_data[OUT_VEGT].data[0] /= cv_veg;
  }
  if (cv_overstory > 0) {
    out_data[OUT_AERO_COND2].data[0] /= cv_overstory;
  }
  if (cv_snow > 0) {
    out_data[OUT_SALBEDO].data[0] /= cv_snow;
    out_data[OUT_SNOW_SURF_TEMP].data[0] /= cv_snow;
    out_data[OUT_SNOW_PACK_TEMP].data[0] /= cv_snow;
  }
  if (cv_glacier > 0) {
    out_data[OUT_GLAC_SURF_TEMP].data[0] /= cv_glacier;
  }

  // Radiative temperature
  out_data[OUT_RAD_TEMP].data[0] = pow(out_data[OUT_RAD_TEMP].data[0],0.25);

  // Aerodynamic conductance and resistance
  if (out_data[OUT_AERO_COND1].data[0] > SMALL) {
    out_data[OUT_AERO_RESIST1].data[0] = 1 / out_data[OUT_AERO_COND1].data[0];
  }
  else {
    out_data[OUT_AERO_RESIST1].data[0] = HUGE_RESIST;
  }
  if (out_data[OUT_AERO_COND2].data[0] > SMALL) {
    out_data[OUT_AERO_RESIST2].data[0] = 1 / out_data[OUT_AERO_COND2].data[0];
  }
  else {
    out_data[OUT_AERO_RESIST2].data[0] = HUGE_RESIST;
  }
  if (out_data[OUT_AERO_COND].data[0] > SMALL) {
    out_data[OUT_AERO_RESIST].data[0] = 1 / out_data[OUT_AERO_COND].data[0];
  }
  else {
    out_data[OUT_AERO_RESIST].data[0] = HUGE_RESIST;
  }

  /*****************************************
    Compute derived variables
   *****************************************/
  // Water balance terms
  out_data[OUT_DELSOILMOIST].data[0] = 0;
  for (int index=0; index<state->options.Nlayer; index++) {
	out_data[OUT_SOIL_LIQ_TOT].data[0] += out_data[OUT_SOIL_LIQ].data[index];
	out_data[OUT_SOIL_ICE_TOT].data[0] += out_data[OUT_SOIL_ICE].data[index];
	out_data[OUT_SOIL_MOIST].data[index] = out_data[OUT_SOIL_LIQ].data[index]+out_data[OUT_SOIL_ICE].data[index];
    out_data[OUT_DELSOILMOIST].data[0] += out_data[OUT_SOIL_MOIST].data[index];
    out_data[OUT_SMLIQFRAC].data[index] = out_data[OUT_SOIL_LIQ].data[index]/out_data[OUT_SOIL_MOIST].data[index];
    out_data[OUT_SMFROZFRAC].data[index] = 1 - out_data[OUT_SMLIQFRAC].data[index];
  }
  if (rec >= 0) {
    out_data[OUT_DELSOILMOIST].data[0] -= cell->save_data.total_soil_moist;
    out_data[OUT_DELSWE].data[0] = out_data[OUT_SWE].data[0] + out_data[OUT_SNOW_CANOPY].data[0] - cell->save_data.swe;
    out_data[OUT_DELINTERCEPT].data[0] = out_data[OUT_WDEW].data[0] - cell->save_data.wdew;
    out_data[OUT_DELSURFSTOR].data[0] = out_data[OUT_SURFSTOR].data[0] - cell->save_data.surfstor;
  }

  // Energy terms
  out_data[OUT_REFREEZE].data[0] = (out_data[OUT_RFRZ_ENERGY].data[0]/Lf)*dt_sec;
  out_data[OUT_R_NET].data[0] = out_data[OUT_NET_SHORT].data[0] + out_data[OUT_NET_LONG].data[0];

  // Save current moisture state for use in next time step
  cell->save_data.total_soil_moist = 0;
  for (int index=0; index<state->options.Nlayer; index++) {
    cell->save_data.total_soil_moist += out_data[OUT_SOIL_MOIST].data[index];
  }
  out_data[OUT_SOIL_MOIST_TOT].data[0] = cell->save_data.total_soil_moist;
  cell->save_data.surfstor = out_data[OUT_SURFSTOR].data[0];
  cell->save_data.swe = out_data[OUT_SWE].data[0] + out_data[OUT_SNOW_CANOPY].data[0];
  cell->save_data.wdew = out_data[OUT_WDEW].data[0];

  /********************
    Check Water Balance 
    ********************/
  inflow  = out_data[OUT_PREC].data[0] + out_data[OUT_LAKE_CHAN_IN].data[0]; // mm over grid cell
  outflow = out_data[OUT_EVAP].data[0] + out_data[OUT_RUNOFF].data[0] + out_data[OUT_BASEFLOW].data[0]; // mm over grid cell
  glac_icebal = out_data[OUT_GLAC_IMBAL].data[0];
  storage = 0.;
  for(int index=0;index<state->options.Nlayer;index++)
    if(state->options.MOISTFRACT)
      storage += (out_data[OUT_SOIL_LIQ].data[index] + out_data[OUT_SOIL_ICE].data[index]) 
	* cell->soil_con.depth[index] * 1000;
    else
      storage += out_data[OUT_SOIL_LIQ].data[index] + out_data[OUT_SOIL_ICE].data[index];
  storage += out_data[OUT_SWE].data[0] + out_data[OUT_SNOW_CANOPY].data[0] + out_data[OUT_WDEW].data[0] + out_data[OUT_SURFSTOR].data[0] + out_data[OUT_GLAC_WAT_STOR].data[0];
  out_data[OUT_WATER_ERROR].data[0] = calc_water_balance_error(rec,inflow,outflow,glac_icebal,storage, state->global_param.nrecs, &cell->cellErrors);
  
  /********************
    Check Energy Balance
  ********************/
  if(state->options.FULL_ENERGY)
    calc_energy_balance_error(rec, out_data[OUT_NET_SHORT].data[0] + out_data[OUT_NET_LONG].data[0],
			      out_data[OUT_LATENT].data[0]+out_data[OUT_LATENT_SUB].data[0],
			      out_data[OUT_SENSIBLE].data[0]+out_data[OUT_ADV_SENS].data[0],
			      out_data[OUT_GRND_FLUX].data[0]+out_data[OUT_DELTAH].data[0]+out_data[OUT_FUSION].data[0],
			      out_data[OUT_ADVECTION].data[0] - out_data[OUT_DELTACC].data[0]
			      - out_data[OUT_SNOW_FLUX].data[0] + out_data[OUT_RFRZ_ENERGY].data[0],
			      state->global_param.nrecs, &cell->cellErrors);



  /******************************************************************************************
    Return to parent function if this was just an initialization of wb and eb storage terms
  ******************************************************************************************/
  if (rec < 0) return(0);



  /********************
    Report T Fallback Occurrences
  ********************/
  if (rec == state->global_param.nrecs-1) {
    fprintf(stderr,"Total number of fallbacks in Tfoliage: %d\n", cell->fallBackStats.Tfoliage_fbcount_total);
    fprintf(stderr,"Total number of fallbacks in Tcanopy: %d\n", cell->fallBackStats.Tcanopy_fbcount_total);
    fprintf(stderr,"Total number of fallbacks in Tsnowsurf: %d\n", cell->fallBackStats.Tsnowsurf_fbcount_total);
    fprintf(stderr,"Total number of fallbacks in Tsurf: %d\n", cell->fallBackStats.Tsurf_fbcount_total);
    fprintf(stderr,"Total number of fallbacks in soil T profile: %d\n", cell->fallBackStats.Tsoil_fbcount_total);
    fprintf(stderr,"Total number of fallbacks in Tglac_surf: %d\n", cell->fallBackStats.Tglacsurf_fbcount_total);
  }

  /********************
    Temporal Aggregation 
    ********************/
  for (int v=0; v<N_OUTVAR_TYPES; v++) {
    if (out_data[v].aggtype == AGG_TYPE_END) {
      for (int i=0; i<out_data[v].nelem; i++) {
        out_data[v].aggdata[i] = out_data[v].data[i];
      }
    }
    else if (out_data[v].aggtype == AGG_TYPE_SUM) {
      for (int i=0; i<out_data[v].nelem; i++) {
        out_data[v].aggdata[i] += out_data[v].data[i];
      }
    }
    else if (out_data[v].aggtype == AGG_TYPE_AVG) {
      for (int i=0; i<out_data[v].nelem; i++) {
        out_data[v].aggdata[i] += out_data[v].data[i]/out_step_ratio;
      }
    }
  }
  out_data[OUT_AERO_RESIST].aggdata[0] = 1/out_data[OUT_AERO_COND].aggdata[0];
  out_data[OUT_AERO_RESIST1].aggdata[0] = 1/out_data[OUT_AERO_COND1].aggdata[0];
  out_data[OUT_AERO_RESIST2].aggdata[0] = 1/out_data[OUT_AERO_COND2].aggdata[0];

  /********************
    Output procedure
    (only execute when we've completed an output interval)
  ********************/
  if (cell->fallBackStats.step_count == out_step_ratio) {

    /***********************************************
      Change of units for ALMA-compliant output
    ***********************************************/
    if (state->options.ALMA_OUTPUT) {
      out_data[OUT_BASEFLOW].aggdata[0] /= out_dt_sec;
      out_data[OUT_EVAP].aggdata[0] /= out_dt_sec;
      out_data[OUT_EVAP_BARE].aggdata[0] /= out_dt_sec;
      out_data[OUT_EVAP_CANOP].aggdata[0] /= out_dt_sec;
      out_data[OUT_INFLOW].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_BF_IN].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_BF_IN_V].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_BF_OUT].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_BF_OUT_V].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_CHAN_IN].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_CHAN_IN_V].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_CHAN_OUT].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_CHAN_OUT_V].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_DSTOR].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_DSTOR_V].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_DSWE].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_DSWE_V].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_EVAP].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_EVAP_V].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_ICE_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_LAKE_PREC_V].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_RCHRG].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_RCHRG_V].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_RO_IN].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_RO_IN_V].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_VAPFLX].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_VAPFLX_V].aggdata[0] /= out_dt_sec;
      out_data[OUT_LAKE_SURF_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_PREC].aggdata[0] /= out_dt_sec;
      out_data[OUT_RAINF].aggdata[0] /= out_dt_sec;
      out_data[OUT_REFREEZE].aggdata[0] /= out_dt_sec;
      out_data[OUT_RUNOFF].aggdata[0] /= out_dt_sec;
      out_data[OUT_SNOW_MELT].aggdata[0] /= out_dt_sec;
      out_data[OUT_SNOWF].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_BLOWING].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_CANOP].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_SNOW].aggdata[0] /= out_dt_sec;
      out_data[OUT_SUB_SNOW].aggdata[0] += out_data[OUT_SUB_CANOP].aggdata[0];
      out_data[OUT_SUB_SURFACE].aggdata[0] /= out_dt_sec;
      out_data[OUT_TRANSP_VEG].aggdata[0] /= out_dt_sec;
      out_data[OUT_BARESOILT].aggdata[0] += KELVIN;
      out_data[OUT_SNOW_PACK_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_SNOW_SURF_TEMP].aggdata[0] += KELVIN;
      for (int index=0; index<state->options.Nlayer; index++) {
        out_data[OUT_SOIL_TEMP].aggdata[index] += KELVIN;
      }
      for (int index=0; index<state->options.Nnode; index++) {
        out_data[OUT_SOIL_TNODE].aggdata[index] += KELVIN;
        out_data[OUT_SOIL_TNODE_WL].aggdata[index] += KELVIN;
      }
      out_data[OUT_SURF_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_VEGT].aggdata[0] += KELVIN;
      out_data[OUT_FDEPTH].aggdata[0] /= 100;
      out_data[OUT_TDEPTH].aggdata[0] /= 100;
      out_data[OUT_DELTACC].aggdata[0] *= out_dt_sec;
      out_data[OUT_DELTAH].aggdata[0] *= out_dt_sec;
      out_data[OUT_AIR_TEMP].aggdata[0] += KELVIN;
      out_data[OUT_PRESSURE].aggdata[0] *= 1000;
      out_data[OUT_VP].aggdata[0] *= 1000;
      out_data[OUT_VPD].aggdata[0] *= 1000;
    }

    /*************
      Write Data
    *************/
    if(rec >= skipyear) {
      output->write_data(out_data, dmy, state->global_param.out_dt, state);

    }

    // Reset the step count
    cell->fallBackStats.step_count = 0;

    // Reset the agg data
    for (int v=0; v<N_OUTVAR_TYPES; v++) {
      for (int i=0; i<out_data[v].nelem; i++) {
        out_data[v].aggdata[i] = 0;
      }
    }

  } // End of output procedure

  return (0);

}

void collect_wb_terms(const hru_data_struct&  cell,
                      const veg_var_struct&    veg_var,
                      const snow_data_struct&  snow,
                      const glac_data_struct& glacier,
                      const lake_var_struct&   lake_var,
                      double            precipitation_mu,
                      double            Cv,
                      double            TreeAdjustFactor,
                      bool              HasVeg,
                      bool              HasGlac,
                      bool              IsWet,
                      double            lakefactor,
                      int               overstory,
                      double           *depth,
                      double           *frost_fract,
                      out_data_struct  *out_data,
                      const ProgramState* state)
{
  double AreaFactor;
  double tmp_evap;
  double tmp_cond1;
  double tmp_cond2;
  double tmp_moist;
  double tmp_ice;
  int index;
  int                     frost_area;

  AreaFactor = Cv * precipitation_mu * TreeAdjustFactor * lakefactor;

  /** record evaporation components **/
  tmp_evap = 0.0;
  for(index=0;index<state->options.Nlayer;index++)
    tmp_evap += cell.layer[index].evap;
  if (HasVeg)
    out_data[OUT_TRANSP_VEG].data[0] += tmp_evap * AreaFactor;
  else 
    out_data[OUT_EVAP_BARE].data[0] += tmp_evap * AreaFactor;
  tmp_evap += snow.vapor_flux * 1000.;
  out_data[OUT_SUB_SNOW].data[0] += snow.vapor_flux * 1000. * AreaFactor; 
  out_data[OUT_SUB_SURFACE].data[0] += snow.surface_flux * 1000. * AreaFactor; 
  out_data[OUT_SUB_BLOWING].data[0] += snow.blowing_flux * 1000. * AreaFactor; 
  if (HasVeg) {
    tmp_evap += snow.canopy_vapor_flux * 1000.;
    out_data[OUT_SUB_CANOP].data[0] += snow.canopy_vapor_flux * 1000. * AreaFactor; 
  }
  if (HasVeg) {
    tmp_evap += veg_var.canopyevap;
    out_data[OUT_EVAP_CANOP].data[0] += veg_var.canopyevap * AreaFactor; 
  }
  if (HasGlac) {
    tmp_evap += glacier.vapor_flux * 1000.;
  }
  out_data[OUT_EVAP].data[0] += tmp_evap * AreaFactor; // mm over gridcell

  /** record potential evap **/
  out_data[OUT_PET_SATSOIL].data[0] += cell.pot_evap[0] * AreaFactor;
  out_data[OUT_PET_H2OSURF].data[0] += cell.pot_evap[1] * AreaFactor;
  out_data[OUT_PET_SHORT].data[0] += cell.pot_evap[2] * AreaFactor;
  out_data[OUT_PET_TALL].data[0] += cell.pot_evap[3] * AreaFactor;
  out_data[OUT_PET_NATVEG].data[0] += cell.pot_evap[4] * AreaFactor;
  out_data[OUT_PET_VEGNOCR].data[0] += cell.pot_evap[5] * AreaFactor;

  /** record saturated area fraction **/
  out_data[OUT_ASAT].data[0] += cell.asat * AreaFactor; 

  /** record runoff **/
  out_data[OUT_RUNOFF].data[0]   += cell.runoff * AreaFactor;

  /** record baseflow **/
  out_data[OUT_BASEFLOW].data[0] += cell.baseflow * AreaFactor; 

  /** record inflow **/
  out_data[OUT_INFLOW].data[0] += (cell.inflow) * AreaFactor;
 
  /** record canopy interception **/
  if (HasVeg) 
    out_data[OUT_WDEW].data[0] += veg_var.Wdew * AreaFactor;

  /** record aerodynamic conductance and resistance **/
  if (cell.aero_resist.surface > SMALL) {
    tmp_cond1 = (1/cell.aero_resist.surface) * AreaFactor;
  } else {
    tmp_cond1 = HUGE_RESIST;
  }
  out_data[OUT_AERO_COND1].data[0] += tmp_cond1;
  if (overstory) {
    if (cell.aero_resist.overstory > SMALL) {
      tmp_cond2 = (1/cell.aero_resist.overstory) * AreaFactor;
    } else {
      tmp_cond2 = HUGE_RESIST;
    }
  } else {
	tmp_cond2 = HUGE_RESIST;
  }
  out_data[OUT_AERO_COND2].data[0] += tmp_cond2;
  if (overstory) {
    out_data[OUT_AERO_COND].data[0] += tmp_cond2;
  } else {
    out_data[OUT_AERO_COND].data[0] += tmp_cond1;
  }

  /** record layer moistures **/
  for(index=0;index<state->options.Nlayer;index++) {
    tmp_moist = cell.layer[index].moist;
#if SPATIAL_FROST
    tmp_ice = 0;
    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
      tmp_ice  += (cell.layer[index].soil_ice[frost_area] * frost_fract[frost_area]);
#else
    tmp_ice   = cell.layer[index].soil_ice;
#endif
    tmp_moist -= tmp_ice;
    if(state->options.MOISTFRACT) {
      tmp_moist /= depth[index] * 1000.;
      tmp_ice /= depth[index] * 1000.;
    }
    out_data[OUT_SOIL_LIQ].data[index] += tmp_moist * AreaFactor;
    out_data[OUT_SOIL_ICE].data[index] += tmp_ice * AreaFactor;
  }
  out_data[OUT_SOIL_WET].data[0] += cell.wetness * AreaFactor;
  out_data[OUT_ROOTMOIST].data[0] += cell.rootmoist * AreaFactor;

  /** record water table position **/
  out_data[OUT_ZWT].data[0] += cell.zwt * AreaFactor;
  out_data[OUT_ZWT2].data[0] += cell.zwt2 * AreaFactor;
  out_data[OUT_ZWT3].data[0] += cell.zwt3 * AreaFactor;
  for(index=0;index<state->options.Nlayer;index++) {
    out_data[OUT_ZWTL].data[index] += cell.layer[index].zwt * AreaFactor;
  }

  /** record layer temperatures **/
  for(index=0;index<state->options.Nlayer;index++) {
    out_data[OUT_SOIL_TEMP].data[index] += cell.layer[index].T * AreaFactor;
  }

  /*****************************
    Record Snow Pack Variables 
  *****************************/
  
  /** record snow water equivalence **/
  out_data[OUT_SWE].data[0] += snow.swq * AreaFactor * 1000.;
  
  /** record snowpack depth **/
  out_data[OUT_SNOW_DEPTH].data[0] += snow.depth * AreaFactor * 100.;
  
  /** record snowpack albedo, temperature **/
  if (snow.swq> 0.0) {
    out_data[OUT_SALBEDO].data[0] += snow.albedo * AreaFactor;
    out_data[OUT_SNOW_SURF_TEMP].data[0] += snow.surf_temp * AreaFactor;
    out_data[OUT_SNOW_PACK_TEMP].data[0] += snow.pack_temp * AreaFactor;
  }

  /** record canopy intercepted snow **/
  if (HasVeg)
    out_data[OUT_SNOW_CANOPY].data[0] += (snow.snow_canopy) * AreaFactor * 1000.;

  /** record snowpack melt **/
  out_data[OUT_SNOW_MELT].data[0] += snow.melt * AreaFactor * 1000.;

  /** record snow cover fraction **/
  out_data[OUT_SNOW_COVER].data[0] += snow.coverage * AreaFactor;

  /************************
  Record Glacier Variables
  ************************/
  if(HasGlac){
    out_data[OUT_GLAC_WAT_STOR].data[0] += glacier.water_storage * AreaFactor * 1000.;
    out_data[OUT_GLAC_AREA].data[0] += AreaFactor;
    out_data[OUT_GLAC_MBAL].data[0] += glacier.mass_balance * AreaFactor * 1000.;
    out_data[OUT_GLAC_IMBAL].data[0] += glacier.ice_mass_balance * AreaFactor * 1000.;
    out_data[OUT_GLAC_ACCUM].data[0] += glacier.accumulation * AreaFactor * 1000.;
    out_data[OUT_GLAC_MELT].data[0] += glacier.melt * AreaFactor * 1000.;
    out_data[OUT_GLAC_SUB].data[0] += glacier.vapor_flux * AreaFactor * 1000.;
    out_data[OUT_GLAC_INFLOW].data[0] += glacier.inflow * AreaFactor * 1000.;
    out_data[OUT_GLAC_OUTFLOW].data[0] += glacier.outflow * AreaFactor * 1000.;
    out_data[OUT_GLAC_OUTFLOW_COEF].data[0] += glacier.outflow_coef * AreaFactor;
  }
}

void collect_eb_terms(const energy_bal_struct& energy,
                      const snow_data_struct&  snow,
                      const glac_data_struct& glacier,
                      const hru_data_struct&  cell_wet,
                      FallBackStats     *fallBackStats,
                      double            Cv,
                      double            AreaFract,
                      double            TreeAdjustFactor,
                      bool              HasVeg,
                      bool              HasGlac,
                      bool              IsWet,
                      double            lakefactor,
                      int               overstory,
                      int               band,
                      double           *depth,
                      double           *dz,
                      double           *frost_fract,
                      double            frost_slope,
                      out_data_struct  *out_data,
                      const ProgramState* state)
{

  double AreaFactor;
  double tmp_fract;
  double rad_temp;
  double surf_temp;
  int index;
  int    frost_area;

  AreaFactor = Cv * TreeAdjustFactor * lakefactor;

  /**********************************
    Record Frozen Soil Variables
  **********************************/

  /** record freezing and thawing front depths **/
  if(state->options.FROZEN_SOIL) {
    for(index = 0; index < MAX_FRONTS; index++) {
      if(IS_VALID(energy.fdepth[index]))
        out_data[OUT_FDEPTH].data[index] += energy.fdepth[index] * AreaFactor * 100.;
      if(IS_VALID(energy.tdepth[index]))
        out_data[OUT_TDEPTH].data[index] += energy.tdepth[index] * AreaFactor * 100.;
    }
  }

  tmp_fract = 0;
#if SPATIAL_FROST
  for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
    if ( cell_wet.layer[0].soil_ice[frost_area] )
      tmp_fract  += frost_fract[frost_area];
#else
  if ( cell_wet.layer[0].soil_ice > 0 )
    tmp_fract   = 1.;
#endif
  out_data[OUT_SURF_FROST_FRAC].data[0] += tmp_fract * AreaFactor;

  tmp_fract = 0;
#if SPATIAL_FROST
  if ( (energy.T[0] + frost_slope / 2.) > 0 ) {
    if ( (energy.T[0] - frost_slope / 2.) <= 0 )
      tmp_fract += linear_interp( 0, (energy.T[0] + frost_slope / 2.), (energy.T[0] - frost_slope / 2.), 1, 0) * AreaFactor;
  }
  else
    tmp_fract += 1 * AreaFactor;
#else
  if ( energy.T[0] <= 0 )
    tmp_fract = 1 * AreaFactor;
#endif

  /**********************************
    Record Energy Balance Variables
  **********************************/

  /** record surface radiative temperature **/
  if ( overstory && snow.snow && !(state->options.LAKES && IsWet)) {
    rad_temp = energy.Tcanopy + KELVIN;
  }
  else
    rad_temp = energy.Tsurf + KELVIN;

  /** record surface skin temperature **/
  surf_temp = energy.Tsurf;

  /** record landcover temperature **/
  if(HasVeg) {
    // landcover is bare soil
    out_data[OUT_BARESOILT].data[0] += (rad_temp-KELVIN) * AreaFactor;
  }
  else {
    // landcover is vegetation
    if ( overstory && !snow.snow )
      // here, rad_temp will be wrong since it will pick the understory temperature
      out_data[OUT_VEGT].data[0] += energy.Tfoliage * AreaFactor;
    else
      out_data[OUT_VEGT].data[0] += (rad_temp-KELVIN) * AreaFactor;
  }

  /** record mean surface temperature [C]  **/
  out_data[OUT_SURF_TEMP].data[0] += surf_temp * AreaFactor;
  
  /** record thermal node temperatures **/
  for(index=0;index<state->options.Nnode;index++) {
    out_data[OUT_SOIL_TNODE].data[index] += energy.T[index] * AreaFactor;
  }
  if (IsWet) {
    for(index=0;index<state->options.Nnode;index++) {
      out_data[OUT_SOIL_TNODE_WL].data[index] = energy.T[index];
    }
  }

  /** record temperature flags  **/
  out_data[OUT_SURFT_FBFLAG].data[0] += energy.Tsurf_fbflag * AreaFactor;
  fallBackStats->Tsurf_fbcount_total += energy.Tsurf_fbcount;
  for (index=0; index<state->options.Nnode; index++) {
    out_data[OUT_SOILT_FBFLAG].data[index] += energy.T_fbflag[index] * AreaFactor;
    fallBackStats->Tsoil_fbcount_total += energy.T_fbcount[index];
  }
  out_data[OUT_SNOWT_FBFLAG].data[0] += snow.surf_temp_fbflag * AreaFactor;
  fallBackStats->Tsnowsurf_fbcount_total += snow.surf_temp_fbcount;
  out_data[OUT_TFOL_FBFLAG].data[0] += energy.Tfoliage_fbflag * AreaFactor;
  fallBackStats->Tfoliage_fbcount_total += energy.Tfoliage_fbcount;
  out_data[OUT_TCAN_FBFLAG].data[0] += energy.Tcanopy_fbflag * AreaFactor;
  fallBackStats->Tcanopy_fbcount_total += energy.Tcanopy_fbcount;
  out_data[OUT_GLAC_TSURF_FBFLAG].data[0] += glacier.surf_temp_fbflag * AreaFactor;
  fallBackStats->Tglacsurf_fbcount_total += glacier.surf_temp_fbcount;

  /** record net shortwave radiation **/
  out_data[OUT_NET_SHORT].data[0] += energy.NetShortAtmos * AreaFactor;

  /** record net longwave radiation **/
  out_data[OUT_NET_LONG].data[0]  += energy.NetLongAtmos * AreaFactor;

  /** record incoming longwave radiation at ground surface (under veg) **/
  if ( snow.snow && overstory )
    out_data[OUT_IN_LONG].data[0] += energy.LongOverIn * AreaFactor;
  else
    out_data[OUT_IN_LONG].data[0] += energy.LongUnderIn * AreaFactor;

  /** record albedo **/
  if ( snow.snow && overstory )
    out_data[OUT_ALBEDO].data[0]    += energy.AlbedoOver * AreaFactor;
  else
    out_data[OUT_ALBEDO].data[0]    += energy.AlbedoUnder * AreaFactor;

  /** record latent heat flux **/
  out_data[OUT_LATENT].data[0]    -= energy.AtmosLatent * AreaFactor;

  /** record latent heat flux from sublimation **/
  out_data[OUT_LATENT_SUB].data[0] -= energy.AtmosLatentSub * AreaFactor;

  /** record sensible heat flux **/
  out_data[OUT_SENSIBLE].data[0]  -= energy.AtmosSensible * AreaFactor;

  /** record ground heat flux (+ heat storage) **/
  out_data[OUT_GRND_FLUX].data[0] -= energy.grnd_flux * AreaFactor;

  /** record heat storage **/
  out_data[OUT_DELTAH].data[0]    -= energy.deltaH * AreaFactor;

  /** record heat of fusion **/
  out_data[OUT_FUSION].data[0]    -= energy.fusion * AreaFactor;

  /** record energy balance error **/
  out_data[OUT_ENERGY_ERROR].data[0] += energy.error * AreaFactor;

  /** record radiative effective temperature [K], 
      emissivities set = 1.0  **/
  out_data[OUT_RAD_TEMP].data[0] += ((rad_temp) * (rad_temp) * (rad_temp) * (rad_temp)) * AreaFactor;
  
  /** record snowpack cold content **/
  out_data[OUT_DELTACC].data[0] += energy.deltaCC * AreaFactor;
  
  /** record snowpack advection **/
  if (snow.snow && overstory)
    out_data[OUT_ADVECTION].data[0] += energy.canopy_advection * AreaFactor;
  out_data[OUT_ADVECTION].data[0] += energy.advection * AreaFactor;
  
  /** record snow energy flux **/
  out_data[OUT_SNOW_FLUX].data[0] += energy.snow_flux * AreaFactor;
  
  /** record refreeze energy **/
  if (snow.snow && overstory)
    out_data[OUT_RFRZ_ENERGY].data[0] += energy.canopy_refreeze * AreaFactor;
  out_data[OUT_RFRZ_ENERGY].data[0] += energy.refreeze_energy * AreaFactor;

  /** record melt energy **/
  out_data[OUT_MELT_ENERGY].data[0] += energy.melt_energy * AreaFactor;

  /** record advected sensible heat energy **/
  if ( !overstory )
    out_data[OUT_ADV_SENS].data[0] -= energy.advected_sensible * AreaFactor;

  /** Record Glacier Variables **/
  if(HasGlac){
    out_data[OUT_GLAC_SURF_TEMP].data[0] += glacier.surf_temp * AreaFactor;
    out_data[OUT_GLAC_DELTACC].data[0] += energy.deltaCC_glac * AreaFactor;
    out_data[OUT_GLAC_FLUX].data[0] += energy.glacier_flux * AreaFactor;
  }
 
  /**********************************
    Record Band-Specific Variables
  **********************************/
  if (AreaFract == 0) {
    throw VICException("Error: AreaFract is zero! cannot divide by 0 (put_data.c)");
  }
  double bandFactor = Cv * lakefactor / AreaFract;

  /** record band snow water equivalent **/
  out_data[OUT_SWE_BAND].data[band] += snow.swq * bandFactor * 1000.;

  /** record band snowpack depth **/
  out_data[OUT_SNOW_DEPTH_BAND].data[band] += snow.depth * bandFactor * 100.;

  /** record band canopy intercepted snow **/
  if (HasVeg)
    out_data[OUT_SNOW_CANOPY_BAND].data[band] += (snow.snow_canopy) * bandFactor * 1000.;

  /** record band snow melt **/
  out_data[OUT_SNOW_MELT_BAND].data[band] += snow.melt * bandFactor;

  /** record band snow coverage **/
  out_data[OUT_SNOW_COVER_BAND].data[band] += snow.coverage * bandFactor;

  /** record band cold content **/
  out_data[OUT_DELTACC_BAND].data[band] += energy.deltaCC * bandFactor;
    
  /** record band advection **/
  out_data[OUT_ADVECTION_BAND].data[band] += energy.advection * bandFactor;
    
  /** record band snow flux **/
  out_data[OUT_SNOW_FLUX_BAND].data[band] += energy.snow_flux * bandFactor;
    
  /** record band refreeze energy **/
  out_data[OUT_RFRZ_ENERGY_BAND].data[band] += energy.refreeze_energy * bandFactor;
    
  /** record band melt energy **/
  out_data[OUT_MELT_ENERGY_BAND].data[band] += energy.melt_energy * bandFactor;

  /** record band advected sensble heat **/
  out_data[OUT_ADV_SENS_BAND].data[band] -= energy.advected_sensible * bandFactor;

  /** record surface layer temperature **/
  out_data[OUT_SNOW_SURFT_BAND].data[band] += snow.surf_temp * bandFactor;

  /** record pack layer temperature **/
  out_data[OUT_SNOW_PACKT_BAND].data[band] += snow.pack_temp * bandFactor;

  /** record latent heat of sublimation **/
  out_data[OUT_LATENT_SUB_BAND].data[band] += energy.latent_sub * bandFactor;

  /** record band net downwards shortwave radiation **/
  out_data[OUT_NET_SHORT_BAND].data[band] += energy.NetShortAtmos * bandFactor;

  /** record band net downwards longwave radiation **/
  out_data[OUT_NET_LONG_BAND].data[band] += energy.NetLongAtmos * bandFactor;

  /** record band albedo **/
  if (snow.snow && overstory)
    out_data[OUT_ALBEDO_BAND].data[band] += energy.AlbedoOver * bandFactor;
  else
    out_data[OUT_ALBEDO_BAND].data[band] += energy.AlbedoUnder * bandFactor;

  /** record band net latent heat flux **/
  out_data[OUT_LATENT_BAND].data[band] -= energy.latent * bandFactor;

  /** record band net sensible heat flux **/
  out_data[OUT_SENSIBLE_BAND].data[band] -= energy.sensible * bandFactor;

  /** record band net ground heat flux **/
  out_data[OUT_GRND_FLUX_BAND].data[band] -= energy.grnd_flux * bandFactor;

  /** record glacier band variables **/
  /** as glaciers only occupy a portion of each band, don't multiply by band factor **/
  if(HasGlac){
    out_data[OUT_GLAC_DELTACC_BAND].data[band] += energy.deltaCC_glac;
    out_data[OUT_GLAC_FLUX_BAND].data[band] += energy.glacier_flux;
    out_data[OUT_GLAC_WAT_STOR_BAND].data[band] += glacier.water_storage * 1000.;
    out_data[OUT_GLAC_AREA_BAND].data[band] += bandFactor;
    out_data[OUT_GLAC_MBAL_BAND].data[band] += glacier.mass_balance * 1000.;
    out_data[OUT_GLAC_IMBAL_BAND].data[band] +=  glacier.ice_mass_balance * 1000.;
    out_data[OUT_GLAC_ACCUM_BAND].data[band] +=  glacier.accumulation * 1000.;
    out_data[OUT_GLAC_MELT_BAND].data[band] +=  glacier.melt * 1000.;
    out_data[OUT_GLAC_SUB_BAND].data[band] +=  glacier.vapor_flux * 1000.;
    out_data[OUT_GLAC_INFLOW_BAND].data[band] +=  glacier.inflow * 1000.;
    out_data[OUT_GLAC_OUTFLOW_BAND].data[band] +=  glacier.outflow * 1000.;
  }

}
