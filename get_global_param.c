#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include <string.h>
 
static char vcid[] = "$Id: get_global_param.c,v 5.22.2.35 2011/12/23 06:57:55 vicadmin Exp $";

void finilizeForceTypes(ProgramState* state) {
  if (state->param_set.TYPE[RAINF].SUPPLIED
      && state->param_set.TYPE[SNOWF].SUPPLIED) {
    /* rainfall and snowfall supplied */
    state->param_set.TYPE[PREC].SUPPLIED =
        state->param_set.TYPE[RAINF].SUPPLIED;
  } else if (state->param_set.TYPE[CRAINF].SUPPLIED
      && state->param_set.TYPE[LSRAINF].SUPPLIED
      && state->param_set.TYPE[CSNOWF].SUPPLIED
      && state->param_set.TYPE[LSSNOWF].SUPPLIED) {
    /* convective and large-scale rainfall and snowfall supplied */
    state->param_set.TYPE[PREC].SUPPLIED =
        state->param_set.TYPE[LSRAINF].SUPPLIED;
  }

  if (state->param_set.TYPE[WIND_E].SUPPLIED
      && state->param_set.TYPE[WIND_N].SUPPLIED) {
    /* specific wind_e and wind_n supplied */
    state->param_set.TYPE[WIND].SUPPLIED =
        state->param_set.TYPE[WIND_E].SUPPLIED;
  }

  if (state->param_set.TYPE[QAIR].SUPPLIED
      && state->param_set.TYPE[PRESSURE].SUPPLIED) {
    /* specific humidity and atm. pressure supplied */
    state->param_set.TYPE[VP].SUPPLIED = state->param_set.TYPE[QAIR].SUPPLIED;
  } else if (state->param_set.TYPE[REL_HUMID].SUPPLIED
      && state->param_set.TYPE[AIR_TEMP].SUPPLIED) {
    /* relative humidity and atm. pressure supplied */
    state->param_set.TYPE[VP].SUPPLIED =
        state->param_set.TYPE[REL_HUMID].SUPPLIED;
  }

  if (!state->param_set.TYPE[VP].SUPPLIED) {
    if (state->param_set.TYPE[QAIR].SUPPLIED) {
      state->param_set.TYPE[VP].SUPPLIED = state->param_set.TYPE[QAIR].SUPPLIED;
    } else if (state->param_set.TYPE[REL_HUMID].SUPPLIED) {
      state->param_set.TYPE[VP].SUPPLIED =
          state->param_set.TYPE[REL_HUMID].SUPPLIED;
    } // end if REL_HUMID supplied
  } // end if VP not supplied

}

int latitudeToIndex(double lat, const ProgramState* state) {
  return std::abs(lat - state->global_param.gridStartLat) / state->global_param.gridStepLat;
}

int longitudeToIndex(double lon, const ProgramState* state) {
  return std::abs(lon - state->global_param.gridStartLon) / state->global_param.gridStepLon;
}

// Calculate the grid cell parameters. This is used for NetCDF outputs
void ProgramState::initGrid(const std::vector<cell_info_struct>& cells) {
  if (cells.size() == 0) {
    throw VICException("Error: cannot run the model with no cells! Make sure that some cells are enabled.");
  } else if (cells.size() == 1) {
    global_param.gridNumLatDivisions = 1;
    global_param.gridNumLonDivisions = 1;
    global_param.gridStartLat = cells[0].soil_con.lat;
    global_param.gridStartLon = cells[0].soil_con.lng;
    global_param.gridStepLat = 0;
    global_param.gridStepLon = 0;
  } else {
    global_param.gridStartLat = cells[0].soil_con.lat;
    global_param.gridStartLon = cells[0].soil_con.lng;
    global_param.gridStepLat = INT_MAX;
    global_param.gridStepLon = INT_MAX;
    double largestLat = cells[0].soil_con.lat;
    double largestLon = cells[0].soil_con.lng;

    for (unsigned int i = 0; i < cells.size(); i++) {
      // Find the minimum lat and long in the grid, use this as start of indices.
      if (cells[i].soil_con.lat < global_param.gridStartLat) {
        global_param.gridStartLat = cells[i].soil_con.lat;
      }
      if (cells[i].soil_con.lng < global_param.gridStartLon) {
        global_param.gridStartLon = cells[i].soil_con.lng;
      }
      if (cells[i].soil_con.lat > largestLat) {
        largestLat = cells[i].soil_con.lat;
      }
      if (cells[i].soil_con.lng > largestLon) {
        largestLon = cells[i].soil_con.lng;
      }
      // Find the smallest distance between cells for lat and lon
      for (unsigned int n = 0; n < cells.size(); n++) {
        double distLat = std::abs(cells[i].soil_con.lat - cells[n].soil_con.lat);
        double distLon = std::abs(cells[i].soil_con.lng - cells[n].soil_con.lng);
        if (distLat != 0 && distLat < global_param.gridStepLat) {
          global_param.gridStepLat = distLat;
        }
        if (distLon != 0 && distLon < global_param.gridStepLon) {
          global_param.gridStepLon = distLon;
        }
      }
    }
    global_param.gridNumLatDivisions = ((largestLat - global_param.gridStartLat) / global_param.gridStepLat) + 1;
    global_param.gridNumLonDivisions = ((largestLon - global_param.gridStartLon) / global_param.gridStepLon) + 1;
  }
}

void ProgramState::init_global_param(filenames_struct *names, const char* global_file_name)
/**********************************************************************
  get_global_param	Keith Cherkauer	            March 1998

  This routine reads the VIC model global control file, getting
  values for global parameters, model options, and debugging controls.

  NOTE: any additions or removals of parameters in this file must also
  be made in display_current_settings.c.

  Modifications:
  7-19-96 Modified to read time step		        KAC
  4-5-98  Modified to read model options and debugging
          controls from a single file                   KAC
  01-20-00 modified to work with new radiation estimation routines,
           new simplified frozen soil moisture, and new new open
           format forcing file rad routines.              KAC
  02-27-01 added reads for lake model parameters          KAC
  04-21-03 added parameters for blowing snow algorithm, printing
           lake variables during debugging and reading Bart's 
           new Arno parameters.                           KAC
  11-May-04 Modified to display compile-time and run-time options
	    if VERBOSE is set to TRUE.						TJB
  13-Oct-04 Added validation for GRND_FLUX option.              		TJB
  01-Nov-04 Added validation for Nnodes with QUICK_FLUX option, as
	    part of fix for QUICK_FLUX state file compatibility.		TJB
  2005-03-08 Added EQUAL_AREA option.						TJB
  2005-03-24 Added ALMA_OUTPUT option.						TJB
  2005-04-07 Fixed state file warning check.					TJB
  2005-Apr-13 Added logic for OUTPUT_FORCE option.				TJB
  2005-Apr-23 Changed ARNO_PARAMS to NIJSSEN2001_BASEFLOW			TJB
  2005-11-29 SAVE_STATE is set in global_param (not at compile time)		GCT
  2005-12-06 Moved setting of statename from open_state_file to here.		GCT
  2005-12-07 Added checks for range of STATEMONTH and STATEDAY			GCT
  2005-12-07 Allow user to use NO_FLUX in addition to NOFLUX for NOFLUX in
             global_param.param.file							GCT
  2006-09-13 Replaced NIJSSEN2001_BASEFLOW with BASEFLOW option.		TJB/GCT
  2006-Sep-23 Implemented flexible output configuration; removed the
              OPTIMIZE and LDAS_OUTPUT options; implemented aggregation of
	      output variables.							TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      This included moving global->statename to filenames->statefile;
	      also added f_path_pfx to store forcing file path and prefix.	TJB
  2006-Nov-07 Removed LAKE_MODEL option.					TJB
  2007-Jan-03 Added ALMA_INPUT option.						TJB
  2007-Jan-15 Added PRT_HEADER option.						TJB
  2007-Apr-03 Added CONTINUEONERROR option.					GCT
  2007-Apr-24 Added EXP_TRANS option.						JCA
  2007-Apr-24 Added IMPLICIT option.						JCA
  2007-Apr-23 Added initialization of global parameters.			TJB
  2007-Apr-23 Added check for FULL_ENERGY if lake model is run.			TJB
  2007-Aug-08 Added EXCESS_ICE option.						JCA
  2007-Sep-14 Added initialization of names->soil_dir.				TJB
  2007-Oct-10 Added validation of dt, start date, end date, and nrecs.		TJB
  2007-Oct-31 Added validation of input/output files.				TJB
  2008-Jan-25 Removed setting of SNOW_STEP = global_param.dt for
	      OUTPUT_FORCE == TRUE.						TJB
  2008-Jan-28 Added check that end date falls AFTER start date.			TJB
  2008-Mar-12 Relocated code validating IMPLICIT and EXCESS_ICE options.	TJB
  2008-Apr-21 Added SNOW_ALBEDO option.						KAC via TJB
  2008-Apr-21 Added SNOW_DENSITY option.					TJB
  2009-Jan-12 Added COMPUTE_TREELINE and JULY_TAVG_SUPPLIED options.		TJB
  2009-Jan-16 Modified aero_resist_used and Ra_used to become arrays of
	      two elements (surface and overstory); added
	      options.AERO_RESIST_CANSNOW.					TJB
  2009-May-17 Added AR_406_LS to options.AERO_RESIST_CANSNOW.			TJB
  2009-May-17 Added options.MIN_LIQ.						TJB
  2009-May-18 Added options.PLAPSE.						TJB
  2009-May-20 Added options.GRND_FLUX_TYPE.					TJB
  2009-May-22 Added TFALLBACK value to options.CONTINUEONERROR.			TJB
  2009-Jun-15 Changed order of options to match global parameter file.		TJB
  2009-Aug-29 Now handles commented lines that are indented.			TJB
  2009-Sep-19 Moved TFALLBACK to its own separate option.			TJB
  2009-Nov-15 Added prohibition of QUICK_FLUX=TRUE when IMPLICIT=TRUE or
	      EXP_TRANS=TRUE.							TJB
  2009-Dec-11 Removed min_liq and options.MIN_LIQ.				TJB
  2010-Apr-28 Replaced GLOBAL_LAI with VEGPARAM_LAI and LAI_SRC.		TJB
  2011-May-31 Removed options.GRND_FLUX.					TJB
  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.					TJB
  2011-Jul-05 Now QUICK_FLUX=FALSE is prohibited in water balance mode
	      (when both FULL_ENERGY and FROZEN_SOIL are FALSE).		TJB
  2011-Nov-04 Added options to access new forcing estimation features.		TJB
**********************************************************************/
{

  char cmdstr[MAXSTRING];
  char optstr[MAXSTRING];
  char flgstr[MAXSTRING];
  char ErrStr[MAXSTRING];
  int  file_num;
  int  field;
  int  i;
  int  tmpstartdate;
  int  tmpenddate;
  int  lastvalidday;
  int  lastday[] = {
            31, /* JANUARY */
            28, /* FEBRUARY */
            31, /* MARCH */
            30, /* APRIL */
            31, /* MAY */
            30, /* JUNE */
            31, /* JULY */
            31, /* AUGUST */
            30, /* SEPTEMBER */
            31, /* OCTOBER */
            30, /* NOVEMBER */
            31, /* DECEMBER */
        } ;
  bool outputTypeSet = false;
  bool outputStateTypeSet = false;

  /** Initialize global parameters (that aren't part of the options struct) **/
  global_param.dt            = INVALID_INT;
  global_param.nrecs         = INVALID_INT;
  global_param.startyear     = INVALID_INT;
  global_param.startmonth    = INVALID_INT;
  global_param.startday      = INVALID_INT;
  global_param.starthour     = INVALID_INT;
  global_param.endyear       = INVALID_INT;
  global_param.endmonth      = INVALID_INT;
  global_param.endday        = INVALID_INT;
  global_param.resolution    = INVALID;
  global_param.measure_h     = 2.0;
  global_param.wind_h        = 10.0;
  global_param.glacierAccumStartYear  = INVALID_INT;
  global_param.glacierAccumStartMonth = INVALID_INT;
  global_param.glacierAccumStartDay   = INVALID_INT;
  global_param.glacierAccumInterval   = INVALID_INT;
  for(i = 0; i < 2; i++) {
    global_param.forceyear[i]  = INVALID_INT;
    global_param.forcemonth[i] = 1;
    global_param.forceday[i]   = 1;
    global_param.forcehour[i]  = 0;
    global_param.forceskip[i]  = 0;
    strcpy(names->f_path_pfx[i],"MISSING");
  }
  file_num             = 0;
  global_param.skipyear      = 0;
  strcpy(names->init_state,   "MISSING");
  global_param.stateyear     = INVALID_INT;
  global_param.statemonth    = INVALID_INT;
  global_param.stateday      = INVALID_INT;
  strcpy(names->statefile,    "MISSING");
  strcpy(names->soil,         "MISSING");
  strcpy(names->soil_dir,     "MISSING");
  strcpy(names->veg,          "MISSING");
  strcpy(names->veglib,       "MISSING");
  strcpy(names->snowband,     "MISSING");
  strcpy(names->lakeparam,    "MISSING");
  strcpy(names->result_dir,   "MISSING");
  strcpy(names->netCDFOutputFileName, "results.nc");
  global_param.out_dt        = INVALID_INT;


  // Open the file
  FILE* gp = open_file(global_file_name, "r");

/** Read through global control file to find parameters **/
  fgets(cmdstr,MAXSTRING,gp);

  while(!feof(gp)) {
    if(cmdstr[0]!='#' && cmdstr[0]!='\n' && cmdstr[0]!='\0') {

      sscanf(cmdstr,"%s",optstr);

      /* Handle case of comment line in which '#' is indented */
      if (optstr[0] == '#') {
        fgets(cmdstr,MAXSTRING,gp);
        continue;
      }

      /*************************************
       Get Model Global Parameters
      *************************************/
      if(strcasecmp("NLAYER",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&options.Nlayer);
      }
      else if(strcasecmp("NODES",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&options.Nnode);
      }
      else if(strcasecmp("TIME_STEP",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.dt);
      }
      else if(strcasecmp("SNOW_STEP",optstr)==0) {
	sscanf(cmdstr,"%*s %d",&options.SNOW_STEP);
      }
      else if(strcasecmp("STARTYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.startyear);
      }
      else if(strcasecmp("STARTMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.startmonth);
      }
      else if(strcasecmp("STARTDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.startday);
      }
      else if(strcasecmp("STARTHOUR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.starthour);
      }
      else if (strcasecmp("GLACIER_ACCUM_START_YEAR", optstr)==0) {
        sscanf(cmdstr,"%*s %d", &global_param.glacierAccumStartYear);
      }
      else if (strcasecmp("GLACIER_ACCUM_START_MONTH", optstr)==0) {
        sscanf(cmdstr,"%*s %d", &global_param.glacierAccumStartMonth);
      }
      else if (strcasecmp("GLACIER_ACCUM_START_DAY", optstr)==0) {
        sscanf(cmdstr,"%*s %d", &global_param.glacierAccumStartDay);
      }
      else if (strcasecmp("GLACIER_ACCUM_INTERVAL", optstr)==0) {
        sscanf(cmdstr,"%*s %d", &global_param.glacierAccumInterval);
      }
      else if(strcasecmp("NRECS",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.nrecs);
      }
      else if(strcasecmp("ENDYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.endyear);
      }
      else if(strcasecmp("ENDMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.endmonth);
      }
      else if(strcasecmp("ENDDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.endday);
      }
      else if(strcasecmp("FULL_ENERGY",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) {
	  options.FULL_ENERGY=TRUE;
	}
	else options.FULL_ENERGY = FALSE;
      }
      else if(strcasecmp("FROZEN_SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) {
	  options.FROZEN_SOIL=TRUE;
          options.QUICK_FLUX = FALSE;
	}
        else options.FROZEN_SOIL = FALSE;
      }
      else if(strcasecmp("QUICK_FLUX",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.QUICK_FLUX=TRUE;
        else options.QUICK_FLUX = FALSE;
      }
      else if(strcasecmp("QUICK_SOLVE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.QUICK_SOLVE=TRUE;
        else options.QUICK_SOLVE = FALSE;
      }
      else if( (strcasecmp("NOFLUX",optstr)==0) || (strcasecmp("NO_FLUX",optstr)==0) ) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.NOFLUX=TRUE;
        else options.NOFLUX = FALSE;
      }
      else if(strcasecmp("IMPLICIT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.IMPLICIT=TRUE;
        else options.IMPLICIT = FALSE;
      }
      else if(strcasecmp("EXP_TRANS",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.EXP_TRANS=TRUE;
        else options.EXP_TRANS = FALSE;
      }
      else if (strcasecmp("SNOW_ALBEDO", optstr)==0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if(strcasecmp("SUN1999",flgstr)==0) options.SNOW_ALBEDO=SUN1999;
        else options.SNOW_ALBEDO = USACE;
      }
      else if (strcasecmp("SNOW_DENSITY", optstr)==0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if(strcasecmp("DENS_SNTHRM",flgstr)==0) options.SNOW_DENSITY=DENS_SNTHRM;
        else options.SNOW_DENSITY = DENS_BRAS;
      }
      else if(strcasecmp("BLOWING",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.BLOWING=TRUE;
        else options.BLOWING = FALSE;
      }
      else if(strcasecmp("DIST_PRCP",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.DIST_PRCP=TRUE;
        else options.DIST_PRCP = FALSE;
      }
      else if(strcasecmp("PREC_EXPT",optstr)==0) {
	sscanf(cmdstr,"%*s %f",&options.PREC_EXPT);
      }
      else if(strcasecmp("CORRPREC",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.CORRPREC=TRUE;
        else options.CORRPREC = FALSE;
      }
      else if(strcasecmp("MIN_WIND_SPEED",optstr)==0) {
	sscanf(cmdstr,"%*s %f",&options.MIN_WIND_SPEED);
      }
      else if(strcasecmp("CONTINUEONERROR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.CONTINUEONERROR=TRUE;
        else options.CONTINUEONERROR = FALSE;
      }
      else if(strcasecmp("COMPUTE_TREELINE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("FALSE",flgstr)==0) options.COMPUTE_TREELINE=FALSE;
        else {
          options.COMPUTE_TREELINE = TRUE;
          options.AboveTreelineVeg = atoi( flgstr );
        }
      }
      else if(strcasecmp("EQUAL_AREA",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.EQUAL_AREA=TRUE;
        else options.EQUAL_AREA = FALSE;
      }
      else if(strcasecmp("RESOLUTION",optstr)==0) {
        sscanf(cmdstr,"%*s %f",&global_param.resolution);
      }
      else if (strcasecmp("AERO_RESIST_CANSNOW", optstr)==0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if(strcasecmp("AR_406",flgstr)==0) options.AERO_RESIST_CANSNOW=AR_406;
        else if(strcasecmp("AR_406_LS",flgstr)==0) options.AERO_RESIST_CANSNOW=AR_406_LS;
        else if(strcasecmp("AR_406_FULL",flgstr)==0) options.AERO_RESIST_CANSNOW=AR_406_FULL;
        else if(strcasecmp("AR_410",flgstr)==0) options.AERO_RESIST_CANSNOW=AR_410;
        else if(strcasecmp("AR_COMBO",flgstr)==0) options.AERO_RESIST_CANSNOW=AR_COMBO;
      }
      else if (strcasecmp("GRND_FLUX_TYPE", optstr)==0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if(strcasecmp("GF_406",flgstr)==0) options.GRND_FLUX_TYPE=GF_406;
        else if(strcasecmp("GF_410",flgstr)==0) options.GRND_FLUX_TYPE=GF_410;
        else if(strcasecmp("GF_FULL",flgstr)==0) options.GRND_FLUX_TYPE=GF_FULL;
      }
      else if(strcasecmp("LW_TYPE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("LW_TVA",flgstr)==0) options.LW_TYPE=LW_TVA;
        else if(strcasecmp("LW_ANDERSON",flgstr)==0) options.LW_TYPE=LW_ANDERSON;
        else if(strcasecmp("LW_BRUTSAERT",flgstr)==0) options.LW_TYPE=LW_BRUTSAERT;
        else if(strcasecmp("LW_SATTERLUND",flgstr)==0) options.LW_TYPE=LW_SATTERLUND;
        else if(strcasecmp("LW_IDSO",flgstr)==0) options.LW_TYPE=LW_IDSO;
        else if(strcasecmp("LW_PRATA",flgstr)==0) options.LW_TYPE=LW_PRATA;
      }
      else if(strcasecmp("LW_CLOUD",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("LW_CLOUD_DEARDORFF",flgstr)==0) options.LW_CLOUD=LW_CLOUD_DEARDORFF;
        else options.LW_CLOUD = LW_CLOUD_BRAS;
      }
      else if(strcasecmp("MTCLIM_SWE_CORR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.MTCLIM_SWE_CORR=TRUE;
        else options.MTCLIM_SWE_CORR = FALSE;
      }
      else if(strcasecmp("PLAPSE",optstr)==0) {
        sscanf(cmdstr,"%*s %s", flgstr);
        if(strcasecmp("FALSE", flgstr) == 0) options.PLAPSE = FALSE;
        else {
	  options.PLAPSE = TRUE;
	}
      }
      else if(strcasecmp("SW_PREC_THRESH",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        options.SW_PREC_THRESH = atof(flgstr);
      }
      else if(strcasecmp("TFALLBACK",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.TFALLBACK=TRUE;
        else options.TFALLBACK = FALSE;
      }
      else if(strcasecmp("VP_INTERP",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.VP_INTERP=TRUE;
        else options.VP_INTERP = FALSE;
      }
      else if(strcasecmp("VP_ITER",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("VP_ITER_NONE",flgstr)==0) options.VP_ITER=VP_ITER_NONE;
        if(strcasecmp("VP_ITER_ALWAYS",flgstr)==0) options.VP_ITER=VP_ITER_ALWAYS;
        if(strcasecmp("VP_ITER_ANNUAL",flgstr)==0) options.VP_ITER=VP_ITER_ANNUAL;
        if(strcasecmp("VP_ITER_CONVERGE",flgstr)==0) options.VP_ITER=VP_ITER_CONVERGE;
        else options.VP_INTERP = VP_ITER_ALWAYS;
      }

      /*************************************
       Define state files
      *************************************/
      else if (strcasecmp("INIT_STATE", optstr) == 0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if (strcasecmp("FALSE", flgstr) == 0)
          options.INIT_STATE = FALSE;
        else {
          options.INIT_STATE = TRUE;
          strcpy(names->init_state, flgstr);
        }
      } else if (strcasecmp("STATENAME", optstr) == 0) {
        sscanf(cmdstr, "%*s %s", names->statefile);
        options.SAVE_STATE = TRUE;
      } else if (strcasecmp("STATEYEAR", optstr) == 0) {
        sscanf(cmdstr, "%*s %d", &global_param.stateyear);
      } else if (strcasecmp("STATEMONTH", optstr) == 0) {
        sscanf(cmdstr, "%*s %d", &global_param.statemonth);
      } else if (strcasecmp("STATEDAY", optstr) == 0) {
        sscanf(cmdstr, "%*s %d", &global_param.stateday);
      } else if (strcasecmp("BINARY_STATE_FILE", optstr) == 0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if (strcasecmp("FALSE", flgstr) == 0)
          options.STATE_FORMAT = StateOutputFormat::ASCII_STATEFILE;
        else
          options.STATE_FORMAT = StateOutputFormat::BINARY_STATEFILE;

        if (outputStateTypeSet) {
          throw VICException(
                             "ERROR: state output format specified more than once. Check for multiples of BINARY_STATE_FILE and STATE_OUTPUT_FORMAT in the global options file.\n");
        }
        outputStateTypeSet = true;
      } else if (strcasecmp("STATE_FORMAT", optstr) == 0) {
        if (outputStateTypeSet) {
          throw VICException(
                             "ERROR: state output format specified more than once. Check for multiples of BINARY_STATE_FILE and STATE_OUTPUT_FORMAT in the global options file.\n");
        }
        outputStateTypeSet = true;
        sscanf(cmdstr, "%*s %s", flgstr);
        if (strcasecmp("BINARY", flgstr) == 0) {
          options.STATE_FORMAT = StateOutputFormat::BINARY_STATEFILE;
        } else if (strcasecmp("ASCII", flgstr) == 0) {
          options.STATE_FORMAT = StateOutputFormat::ASCII_STATEFILE;
        } else if (strcasecmp("NETCDF", flgstr) == 0) {
          options.STATE_FORMAT = StateOutputFormat::NETCDF_STATEFILE;
#if !NETCDF_OUTPUT_AVAILABLE
          throw VICException("If the NETCDF output is enabled, then VIC must be built with the NETCDF_OUTPUT_AVAILABLE define set to true for this support (In user_def.h)! Please recompile with NETCDF_OUTPUT_AVAILABLE TRUE or change the STATE_FORMAT type (in the global input file) to BINARY or ASCII");
#endif
        } else {
          fprintf(stderr,
                  "Warning, input for option STATE_FORMAT was expecting either BINARY, ASCII, or NETCDF, but received: \"%s\"\n",
                  optstr);
          fprintf(stderr, "STATE_FORMAT will default to ASCII\n");
          options.STATE_FORMAT = StateOutputFormat::ASCII_STATEFILE;
        }


      } else if (strcasecmp("MAX_MEMORY", optstr) == 0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        options.MAX_MEMORY = atof(flgstr);  // Conversion atof defaults to 0.0 on error (which is consistent with our default value).
      }

      /*************************************
       Define forcing files
      *************************************/
      else if (strcasecmp("FORCING1", optstr) == 0) {
        if (strcmp(names->f_path_pfx[0], "MISSING") != 0)
          nrerror(
                  "Tried to define FORCING1 twice, if you want to use two forcing files, the second must be defined as FORCING2");
        sscanf(cmdstr, "%*s %s", names->f_path_pfx[0]);
        file_num = 0;
        field = 0;
      } else if (strcasecmp("FORCING2", optstr) == 0) {
        sscanf(cmdstr, "%*s %s", names->f_path_pfx[1]);
        if (strcasecmp("FALSE", names->f_path_pfx[1]) == 0)
          strcpy(names->f_path_pfx[1], "MISSING");
        file_num = 1;
        field = 0;
      } else if (strcasecmp("FORCE_FORMAT", optstr) == 0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if (strcasecmp(flgstr, "BINARY") == 0)
          param_set.FORCE_FORMAT[file_num] = BINARY;
        else if (strcasecmp(flgstr, "ASCII") == 0)
          param_set.FORCE_FORMAT[file_num] = ASCII;
        else if (strcasecmp(flgstr, "NETCDF") == 0)
          param_set.FORCE_FORMAT[file_num] = NETCDF;
        else
          nrerror("FORCE_FORMAT must be \"NETCDF\", \"ASCII\", or \"BINARY\".");
      } else if (strcasecmp("FORCE_ENDIAN", optstr) == 0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if (strcasecmp(flgstr, "LITTLE") == 0)
          param_set.FORCE_ENDIAN[file_num] = LITTLE;
        else if (strcasecmp(flgstr, "BIG") == 0)
          param_set.FORCE_ENDIAN[file_num] = BIG;
        else
          nrerror("FORCE_ENDIAN must be either BIG or LITTLE.");
      } else if (strcasecmp("N_TYPES", optstr) == 0) {
        sscanf(cmdstr,"%*s %d",&param_set.N_TYPES[file_num]);
      }
      else if(strcasecmp("FORCE_TYPE",optstr)==0) {
        get_force_type(cmdstr,file_num,&field, this);
      }
      else if(strcasecmp("FORCE_DT",optstr)==0) {
        sscanf(cmdstr,"%*s %d ", &param_set.FORCE_DT[file_num]);
      }
      else if(strcasecmp("FORCEYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.forceyear[file_num]);
      }
      else if(strcasecmp("FORCEMONTH",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.forcemonth[file_num]);
      }
      else if(strcasecmp("FORCEDAY",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.forceday[file_num]);
      }
      else if(strcasecmp("FORCEHOUR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.forcehour[file_num]);
      }
      else if(strcasecmp("GRID_DECIMAL",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&options.GRID_DECIMAL);
      }
      else if(strcasecmp("WIND_H",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&global_param.wind_h);
      }
      else if(strcasecmp("MEASURE_H",optstr)==0) {
        sscanf(cmdstr,"%*s %lf",&global_param.measure_h);
      }
      else if(strcasecmp("ALMA_INPUT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.ALMA_INPUT=TRUE;
        else options.ALMA_INPUT = FALSE;
      }
      else if(strcasecmp("OUTPUT_FORCE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.OUTPUT_FORCE=TRUE;
        else options.OUTPUT_FORCE = FALSE;
      }

      /*************************************
       Define parameter files
      *************************************/

      else if(strcasecmp("SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->soil);
      }
      else if(strcasecmp("ARC_SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.ARC_SOIL=TRUE;
        else options.ARC_SOIL = FALSE;
      }
      else if(strcasecmp("SOIL_DIR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->soil_dir);
      }
      else if (strcasecmp("ARNO_PARAMS", optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) {
          nrerror("Please change \"ARNO_PARAMS  TRUE\" to \"BASEFLOW  NIJSSEN2001\" in your global parameter file.");
        }
        else {
          nrerror("Please change \"ARNO_PARAMS  FALSE\" to \"BASEFLOW  ARNO\" in your global parameter file.");
        }
      }
      else if (strcasecmp("NIJSSEN2001_BASEFLOW", optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) {
          nrerror("Please change \"NIJSSEN2001_BASEFLOW  TRUE\" to \"BASEFLOW  NIJSSEN2001\" in your global parameter file.");
        }
        else {
          nrerror("Please change \"NIJSSEN2001_BASEFLOW  FALSE\" to \"BASEFLOW  ARNO\" in your global parameter file.");
        }
      }
      else if (strcasecmp("BASEFLOW", optstr)==0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if(strcasecmp("NIJSSEN2001",flgstr)==0) options.BASEFLOW=NIJSSEN2001;
        else options.BASEFLOW = ARNO;
      }
      else if(strcasecmp("JULY_TAVG_SUPPLIED",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("FALSE",flgstr)==0) options.JULY_TAVG_SUPPLIED=FALSE;
	else options.JULY_TAVG_SUPPLIED=TRUE;
      }
      else if(strcasecmp("ORGANIC_FRACT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("FALSE",flgstr)==0) options.ORGANIC_FRACT=FALSE;
        else options.ORGANIC_FRACT=TRUE;
      }
      else if(strcasecmp("VEGLIB",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->veglib);
      }
      else if(strcasecmp("VEGPARAM",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->veg);
      }
      else if(strcasecmp("GLOBAL_LAI",optstr)==0) {
        fprintf(stderr, "WARNING: GLOBAL_LAI has been replaced by 2 new options: VEGPARAM_LAI (whether the vegparam file contains LAI values) and LAI_SRC (where to get LAI values).\n"); 
        fprintf(stderr, "\"GLOBAL_LAI  TRUE\" should now be: \"VEGPARAM_LAI  TRUE\" and \"LAI_SRC  LAI_FROM_VEGPARAM\".\n"); 
        fprintf(stderr, "\"GLOBAL_LAI  FALSE\" should now be: \"LAI_SRC  LAI_FROM_VEGLIB\".\n"); 
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) {
//          nrerror("Please replace \"GLOBAL_LAI  TRUE\" with the following in your global parameter file:\n\"VEGPARAM_LAI  TRUE\"\n\"LAI_SRC  LAI_FROM_VEGPARAM\"");
          options.VEGPARAM_LAI=TRUE;
          options.LAI_SRC=LAI_FROM_VEGPARAM;
        }
        else {
//          nrerror("Please replace \"GLOBAL_LAI  FALSE\" with the following in your global parameter file:\n\"LAI_SRC  LAI_FROM_VEGLIB\"");
          options.LAI_SRC=LAI_FROM_VEGLIB;
        }
      }
      else if(strcasecmp("VEGPARAM_LAI",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.VEGPARAM_LAI=TRUE;
        else options.VEGPARAM_LAI = FALSE;
      }
      else if(strcasecmp("LAI_SRC",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("LAI_FROM_VEGPARAM",flgstr)==0)
          options.LAI_SRC=LAI_FROM_VEGPARAM;
        else
          options.LAI_SRC=LAI_FROM_VEGLIB;
      }
      else if (strcasecmp("GLACIER_ID", optstr) == 0) {
        sscanf(cmdstr,"%*s %d", &options.GLACIER_ID);
      }
      else if(strcasecmp("ROOT_ZONES",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&options.ROOT_ZONES);
      }
      else if(strcasecmp("SNOW_BAND",optstr)==0) {
        sscanf(cmdstr,"%*s %d %s",&options.SNOW_BAND,names->snowband);
      }
      else if(strcasecmp("LAKES",optstr)==0) {
        sscanf(cmdstr, "%*s %s", flgstr);
        if (strcasecmp("FALSE", flgstr) == 0)
          options.LAKES = FALSE;
        else {
          options.LAKES = TRUE;
          strcpy(names->lakeparam, flgstr);
        }
      }
      else if(strcasecmp("LAKE_PROFILE",optstr)==0) {
        sscanf(cmdstr,"%*s %s", flgstr);
        if(strcasecmp("FALSE", flgstr) == 0) options.LAKE_PROFILE = FALSE;
        else {
	  options.LAKE_PROFILE = TRUE;
	}
      }

      /*************************************
       Define output files
      *************************************/
      else if(strcasecmp("RESULT_DIR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",names->result_dir);
      }
      else if (strcasecmp("NETCDF_OUTPUT_FILENAME", optstr) == 0) {
        sscanf(cmdstr, "%*s %s", names->netCDFOutputFileName);
      }
      else if(strcasecmp("OUT_STEP",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.out_dt);
      }
      else if(strcasecmp("SKIPYEAR",optstr)==0) {
        sscanf(cmdstr,"%*s %d",&global_param.skipyear);
      }
      else if(strcasecmp("COMPRESS",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.COMPRESS=TRUE;
        else options.COMPRESS = FALSE;
      }
      else if(strcasecmp("BINARY_OUTPUT",optstr)==0) {
        if (outputTypeSet) {
          throw VICException("ERROR: output format specified more than once. Check for multiples of BINARY_OUTPUT and OUTPUT_FORMAT in the global options file.\n");
        }
        outputTypeSet = true;
        fprintf(stderr, "Warning: the BINARY_OUTPUT option is now deprecated. Instead, use OUTPUT_FORMAT set to BINARY, ASCII, or NETCDF\n");
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.OUTPUT_FORMAT = OutputFormat::BINARY_FORMAT;
        else options.OUTPUT_FORMAT = OutputFormat::ASCII_FORMAT;
      }
      else if (strcasecmp("OUTPUT_FORMAT", optstr) == 0) {
        if (outputTypeSet) {
          throw VICException("ERROR: output format specified more than once. Check for multiples of BINARY_OUTPUT and OUTPUT_FORMAT in the global options file.\n");
        }
        outputTypeSet = true;
        sscanf(cmdstr, "%*s %s", flgstr);
        if (strcasecmp("BINARY", flgstr) == 0) {
          options.OUTPUT_FORMAT = OutputFormat::BINARY_FORMAT;
        } else if (strcasecmp("ASCII", flgstr) == 0) {
          options.OUTPUT_FORMAT = OutputFormat::ASCII_FORMAT;
        } else if (strcasecmp("NETCDF", flgstr) == 0) {
          options.OUTPUT_FORMAT = OutputFormat::NETCDF_FORMAT;
#if !NETCDF_OUTPUT_AVAILABLE
          throw VICException("If the NETCDF output is enabled, then VIC must be built with the NETCDF_OUTPUT_AVAILABLE define set to true for this support (In user_def.h)! Please recompile with NETCDF_OUTPUT_AVAILABLE TRUE or change the OUTPUT_FORMAT type (in the global input file) to BINARY or ASCII");
#endif
        } else {
          fprintf(stderr, "Warning, input for option OUTPUT_FORMAT was expecting either BINARY, ASCII, or NETCDF, but received: \"%s\"\n", optstr);
          fprintf(stderr, "OUTPUT_FORMAT will default to ASCII\n");
          options.OUTPUT_FORMAT = OutputFormat::ASCII_FORMAT;
        }
      }
      else if(strcasecmp("ALMA_OUTPUT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.ALMA_OUTPUT=TRUE;
        else options.ALMA_OUTPUT = FALSE;
      }
      else if(strcasecmp("MOISTFRACT",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.MOISTFRACT=TRUE;
        else options.MOISTFRACT = FALSE;
      }
      else if(strcasecmp("PRT_HEADER",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.PRT_HEADER=TRUE;
        else options.PRT_HEADER = FALSE;
      }
      else if(strcasecmp("PRT_SNOW_BAND",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) options.PRT_SNOW_BAND=TRUE;
        else options.PRT_SNOW_BAND = FALSE;
      }

      /*************************************
       Define output file contents
      *************************************/
      else if(strcasecmp("N_OUTFILES",optstr)==0) {
        ; // do nothing
      }
      else if(strcasecmp("OUTFILE",optstr)==0) {
        ; // do nothing
      }
      else if(strcasecmp("OUTVAR",optstr)==0) {
        ; // do nothing
      }

      /****************************************************
       Define netCDF output format global attributes files
       ****************************************************/
      else if (strcasecmp("NETCDF_ATTRIBUTE", optstr) == 0) {
        char attributeValue[MAXSTRING];
        sscanf(cmdstr, "%*s %s %[^\n]", flgstr, attributeValue);
        global_param.netCDFGlobalAttributes.push_back(std::pair<std::string, std::string>(flgstr, attributeValue));
      }

      /******************************
       Get Model Debugging Options
       ****************************/

#if LINK_DEBUG
      else if(strcasecmp("PRT_FLUX",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_FLUX=TRUE;
        else debug.PRT_FLUX = FALSE;
      }
      else if(strcasecmp("PRT_BALANCE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_BALANCE=TRUE;
        else debug.PRT_BALANCE = FALSE;
      }
      else if(strcasecmp("PRT_SOIL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_SOIL=TRUE;
        else debug.PRT_SOIL = FALSE;
      }
      else if(strcasecmp("PRT_VEGE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_VEGE=TRUE;
        else debug.PRT_VEGE = FALSE;
      }
      else if(strcasecmp("PRT_GLOBAL",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_GLOBAL=TRUE;
        else debug.PRT_GLOBAL = FALSE;
      }
      else if(strcasecmp("PRT_ATMOS",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_ATMOS=TRUE;
        else debug.PRT_ATMOS = FALSE;
      }
      else if(strcasecmp("PRT_SNOW",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_SNOW=TRUE;
        else debug.PRT_SNOW = FALSE;
      }
      else if(strcasecmp("PRT_MOIST",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_MOIST=TRUE;
        else debug.PRT_MOIST = FALSE;
      }
      else if(strcasecmp("PRT_LAKE",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_LAKE=TRUE;
        else debug.PRT_LAKE = FALSE;
      }
      else if(strcasecmp("PRT_TEMP",optstr)==0) {
        sscanf(cmdstr,"%*s %s",flgstr);
        if(strcasecmp("TRUE",flgstr)==0) debug.PRT_TEMP=TRUE;
        else debug.PRT_TEMP = FALSE;
      }
      else if(strcasecmp("DEBUG_DIR",optstr)==0) {
        sscanf(cmdstr,"%*s %s",debug.debug_dir);
      }
#endif

      /***********************************
        Unrecognized Global Parameter Flag
        ***********************************/
      else {
	fprintf(stderr,"WARNING: Unrecognized option in the global parameter file:\n\t%s is unknown - check your spelling\n", optstr);
      }
    }
    fgets(cmdstr,MAXSTRING,gp);
  }
  fclose(gp);

  finilizeForceTypes(this);

  /******************************************
    Check for undefined required parameters
  ******************************************/

  // Validate model time step
  if (IS_INVALID(global_param.dt))
    nrerror("Model time step has not been defined.  Make sure that the global file defines TIME_STEP.");
  else if (global_param.dt < 1) {
    sprintf(ErrStr,"The specified model time step (%d) < 1 hour.  Make sure that the global file defines a positive number of hours for TIME_STEP.",global_param.dt);
    nrerror(ErrStr);
  }

  // Validate the output step
  if (global_param.out_dt == 0 || IS_INVALID(global_param.out_dt)) {
    global_param.out_dt = global_param.dt;
  }
  else if (global_param.out_dt < global_param.dt || global_param.out_dt > 24 || (float)global_param.out_dt/(float)global_param.dt != (float)(global_param.out_dt/global_param.dt)){
    nrerror("Invalid output step specified.  Output step must be an integer multiple of the model time step; >= model time step and <= 24");
  }

  // Validate SNOW_STEP and set NR and NF
  if (global_param.dt < 24 && global_param.dt != options.SNOW_STEP)
    nrerror("If the model step is smaller than daily, the snow model should run at the same time step as the rest of the model.");
  if (global_param.dt % options.SNOW_STEP != 0 || options.SNOW_STEP > global_param.dt)
    nrerror("SNOW_STEP should be <= TIME_STEP and divide TIME_STEP evenly ");
  NF = global_param.dt/options.SNOW_STEP;
  if (NF == 1)
    NR = 0;
  else
    NR = NF;

  // Validate simulation start date
  if (IS_INVALID(global_param.startyear))
    nrerror("Simulation start year has not been defined.  Make sure that the global file defines STARTYEAR.");
  else if (global_param.startyear < 0) {
    sprintf(ErrStr,"The specified simulation start year (%d) < 0.  Make sure that the global file defines a positive integer for STARTYEAR.",global_param.startyear);
    nrerror(ErrStr);
  }
  if (IS_INVALID(global_param.startmonth))
    nrerror("Simulation start month has not been defined.  Make sure that the global file defines STARTMONTH.");
  else if (global_param.startmonth < 0) {
    sprintf(ErrStr,"The specified simulation start month (%d) < 0.  Make sure that the global file defines a positive integer for STARTMONTH.",global_param.startmonth);
    nrerror(ErrStr);
  }
  if (IS_INVALID(global_param.startday))
    nrerror("Simulation start day has not been defined.  Make sure that the global file defines STARTDAY.");
  else if (global_param.startday < 0) {
    sprintf(ErrStr,"The specified simulation start day (%d) < 0.  Make sure that the global file defines a positive integer for STARTDAY.",global_param.startday);
    nrerror(ErrStr);
  }
  if (IS_INVALID(global_param.starthour)) {
    if (global_param.dt == 24)
      global_param.starthour = 0;
    else
      nrerror("Simulation start hour has not been defined, yet model time step is less than 24 hours.  Make sure that the global file defines STARTHOUR.");
  }
  else if (global_param.starthour < 0) {
    sprintf(ErrStr,"The specified simulation start hour (%d) < 0.  Make sure that the global file defines a positive integer for STARTHOUR.",global_param.starthour);
    nrerror(ErrStr);
  }

  // Validate simulation end date and/or number of timesteps
  if (IS_INVALID(global_param.nrecs) && IS_INVALID(global_param.endyear) && IS_INVALID(global_param.endmonth) && IS_INVALID(global_param.endday))
    nrerror("The model global file MUST define EITHER the number of records to simulate (NRECS), or the year (ENDYEAR), month (ENDMONTH), and day (ENDDAY) of the last full simulation day");
  else if (IS_INVALID(global_param.nrecs)) {
    if (IS_INVALID(global_param.endyear))
      nrerror("Simulation end year has not been defined.  Make sure that the global file defines ENDYEAR.");
    else if (global_param.endyear < 0) {
      sprintf(ErrStr,"The specified simulation end year (%d) < 0.  Make sure that the global file defines a positive integer for ENDYEAR.",global_param.endyear);
      nrerror(ErrStr);
    }
    if (IS_INVALID(global_param.endmonth))
      nrerror("Simulation end month has not been defined.  Make sure that the global file defines ENDMONTH.");
    else if (global_param.endmonth < 0) {
      sprintf(ErrStr,"The specified simulation end month (%d) < 0.  Make sure that the global file defines a positive integer for ENDMONTH.",global_param.endmonth);
      nrerror(ErrStr);
    }
    if (IS_INVALID(global_param.endday))
      nrerror("Simulation end day has not been defined.  Make sure that the global file defines ENDDAY.");
    else if (global_param.endday < 0) {
      sprintf(ErrStr,"The specified simulation end day (%d) < 0.  Make sure that the global file defines a positive integer for ENDDAY.",global_param.endday);
      nrerror(ErrStr);
    }
    tmpstartdate = global_param.startyear*10000 + global_param.startmonth*100 + global_param.startday;
    tmpenddate = global_param.endyear*10000 + global_param.endmonth*100 + global_param.endday;
    if (tmpenddate < tmpstartdate) {
      sprintf(ErrStr,"The specified simulation end date (%04d-%02d-%02d) is EARLIER than the specified start date (%04d-%02d-%02d).",global_param.endyear,global_param.endmonth,global_param.endday,global_param.startyear,global_param.startmonth,global_param.startday);
      nrerror(ErrStr);
    }
  }
  else if (global_param.nrecs < 1) {
    sprintf(ErrStr,"The specified duration of simulation (%d) < 1 time step.  Make sure that the global file defines a positive integer for NRECS.",global_param.nrecs);
    nrerror(ErrStr);
  }

  // Validate forcing files and variables
  if ( strcmp ( names->f_path_pfx[0], "MISSING" ) == 0 )
    nrerror("No forcing file has been defined.  Make sure that the global file defines FORCING1.");
  for(i=0;i<2;i++) {
    if ( i == 0 || (i == 1 && IS_VALID(param_set.N_TYPES[i])) ) {
      if (IS_INVALID(param_set.N_TYPES[i])) {
        sprintf(ErrStr,"Need to specify the number forcing variables types in forcing file %d.", i);
        nrerror(ErrStr);
      }
      if (IS_INVALID(param_set.FORCE_FORMAT[i])) {
        sprintf(ErrStr,"Need to specify the INPUT_FORMAT (ASCII or BINARY) for forcing file %d.",i);
        nrerror(ErrStr);
      }
      if (IS_INVALID(param_set.FORCE_INDEX[i][param_set.N_TYPES[i]-1])) {
        sprintf(ErrStr,"Did not define enough forcing variables in forcing file %d.",i);
        nrerror(ErrStr);
      }
      if(IS_INVALID(param_set.FORCE_DT[i])) {
        sprintf(ErrStr,"Must define time steps (FORCE_DT <dt>) in control file for focing file %d.",file_num);
        nrerror(ErrStr);
      }
    }
  }
  if(IS_VALID(param_set.N_TYPES[1]) && IS_INVALID(global_param.forceyear[1])) {
    global_param.forceyear[1] = global_param.forceyear[0];
    global_param.forcemonth[1] = global_param.forcemonth[0];
    global_param.forceday[1] = global_param.forceday[0];
    global_param.forcehour[1] = global_param.forcehour[0];
    global_param.forceskip[1] = 0;
  }

  // Validate result directory
  if ( strcmp ( names->result_dir, "MISSING" ) == 0 )
    nrerror("No results directory has been defined.  Make sure that the global file defines the result directory on the line that begins with \"RESULT_DIR\".");

  // Validate soil parameter file information
  if ( strcmp ( names->soil, "MISSING" ) == 0 )
    nrerror("No soil parameter file has been defined.  Make sure that the global file defines the soil parameter file on the line that begins with \"SOIL\".");
  if (options.ARC_SOIL && strcmp ( names->soil_dir, "MISSING" ) == 0)
    nrerror("\"ARC_SOIL\" was specified as TRUE, but no soil parameter directory (\"SOIL_DIR\") has been defined.  Make sure that the global file defines the soil parameter directory on the line that begins with \"SOIL_DIR\".");

  /*******************************************************************************
    Validate parameters required for normal simulations but NOT for OUTPUT_FORCE
  *******************************************************************************/

  if (!options.OUTPUT_FORCE) {

    // Validate veg parameter information
    if ( strcmp ( names->veg, "MISSING" ) == 0 )
      nrerror("No vegetation parameter file has been defined.  Make sure that the global file defines the vegetation parameter file on the line that begins with \"VEGPARAM\".");
    if ( strcmp ( names->veglib, "MISSING" ) == 0 )
      nrerror("No vegetation library file has been defined.  Make sure that the global file defines the vegetation library file on the line that begins with \"VEGLIB\".");
    if(IS_INVALID(options.ROOT_ZONES) || options.ROOT_ZONES<0)
      nrerror("ROOT_ZONES must be defined to a positive integer greater than 0, in the global control file.");
    if (options.LAI_SRC == LAI_FROM_VEGPARAM && !options.VEGPARAM_LAI) {
      sprintf(ErrStr, "\"LAI_SRC\" was specified as \"LAI_FROM_VEGPARAM\", but \"VEGPARAM_LAI\" was set to \"FALSE\" in the global parameter file.  If you want VIC to read LAI values from the vegparam file, you MUST make sure the veg param file contains 1 line of 12 monthly LAI values for EACH veg tile in EACH grid cell, and you MUST specify \"VEGPARAM_LAI\" as \"TRUE\" in the global parameter file.  Alternatively, if you want VIC to read LAI values from the veg library file, set \"LAI_SRC\" ro \"LAI_FROM_VEGLIB\" in the global parameter file.  In either case, the setting of \"VEGPARAM_LAI\" must be consistent with the contents of the veg param file (i.e. whether or not it contains LAI values).");
      nrerror(ErrStr);
    }

    // Validate the elevation band file information
    if(options.SNOW_BAND > 1) {
      if ( strcmp ( names->snowband, "MISSING" ) == 0 ) {
        sprintf(ErrStr, "\"SNOW_BAND\" was specified with %d elevation bands, but no elevation band file has been defined.  Make sure that the global file defines the elevation band file on the line that begins with \"SNOW_BAND\" (after the number of bands).", options.SNOW_BAND);
        nrerror(ErrStr);
      }
      if(options.SNOW_BAND > MAX_BANDS) {
        sprintf(ErrStr,"Global file wants more snow bands (%d) than are defined by MAX_BANDS (%d).  Edit user_def.h and recompile.",options.SNOW_BAND,MAX_BANDS);
        nrerror(ErrStr);
      }
    }
    else if (options.SNOW_BAND <= 0) {
      sprintf(ErrStr,"Invalid number of elevation bands specified in global file (%d).  Number of bands must be >= 1.",options.SNOW_BAND);
      nrerror(ErrStr);
    }

    // Validate the input state file information
    if( options.INIT_STATE ) {
      if ( strcmp ( names->init_state, "MISSING" ) == 0 )
        nrerror("\"INIT_STATE\" was specified, but no input state file has been defined.  Make sure that the global file defines the inputstate file on the line that begins with \"INIT_STATE\".");
    }

    // Validate the output state file information
    if( options.SAVE_STATE ) {
      if ( strcmp ( names->statefile, "MISSING" ) == 0)
        nrerror("\"SAVE_STATE\" was specified, but no output state file has been defined.  Make sure that the global file defines the output state file on the line that begins with \"SAVE_STATE\".");
      if ( IS_INVALID(global_param.stateyear) || IS_INVALID(global_param.statemonth) || IS_INVALID(global_param.stateday) )  {
        sprintf(ErrStr,"Incomplete specification of the date to save state for state file (%s).\nSpecified date (yyyy-mm-dd): %04d-%02d-%02d\nMake sure STATEYEAR, STATEMONTH, and STATEDAY are set correctly in your global parameter file.\n", names->statefile, global_param.stateyear, global_param.statemonth, global_param.stateday);
        nrerror(ErrStr);
      }
      // Check for month, day in range
      lastvalidday = lastday[global_param.statemonth - 1];
      if ( global_param.statemonth == 2 ) {
        if ( (global_param.stateyear % 4) == 0 && ( (global_param.stateyear % 100) != 0 || (global_param.stateyear % 400) == 0 ) ){
          lastvalidday = 29;
        }
      }
      if ( global_param.stateday > lastvalidday || global_param.statemonth > 12 || global_param.statemonth < 1 || global_param.stateday > 31 || global_param.stateday < 1 ){
        sprintf(ErrStr,"Unusual specification of the date to save state for state file (%s).\nSpecified date (yyyy-mm-dd): %04d-%02d-%02d\nMake sure STATEYEAR, STATEMONTH, and STATEDAY are set correctly in your global parameter file.\n", names->statefile, global_param.stateyear, global_param.statemonth, global_param.stateday);
        nrerror(ErrStr);
      }
    }
    // Set the statename here to be able to compare with INIT_STATE name
    if( options.SAVE_STATE ) {
      sprintf(names->statefile,"%s_%04i%02i%02i", names->statefile,
              global_param.stateyear, global_param.statemonth, global_param.stateday);
    }
    if( options.INIT_STATE && options.SAVE_STATE && (strcmp( names->init_state, names->statefile ) == 0))  {
      sprintf(ErrStr,"The save state file (%s) has the same name as the initialize state file (%s).  The initialize state file will be destroyed when the save state file is opened.", names->statefile, names->init_state);
      nrerror(ErrStr);
    }

    // Validate soil parameter/simulation mode combinations
    if(options.QUICK_FLUX) {
      if(options.Nnode != 3) {
        fprintf(stderr,"WARNING: To run the model QUICK_FLUX=TRUE, you must define exactly 3 soil thermal nodes.  Currently Nnodes is set to %d.  Setting Nnodes to 3.\n",options.Nnode);
        options.Nnode = 3;
      }
      if(options.IMPLICIT || options.EXP_TRANS) {
        sprintf(ErrStr,"To run the model with QUICK_FLUX=TRUE, you cannot have IMPLICIT=TRUE or EXP_TRANS=TRUE.");
        nrerror(ErrStr);
      }
    }
    if(!options.QUICK_FLUX && !(options.FULL_ENERGY || options.FROZEN_SOIL)) {
      sprintf(ErrStr,"To run the model in water balance mode (both FULL_ENERGY and FROZEN_SOIL are FALSE) you MUST set QUICK_FLUX to TRUE (or leave QUICK_FLUX out of your global parameter file).");
      nrerror(ErrStr);
    }
    if((options.FULL_ENERGY || options.FROZEN_SOIL) && options.Nlayer<3) {
      sprintf(ErrStr,"You must define at least 3 soil moisture layers to run the model in FULL_ENERGY or FROZEN_SOIL modes.  Currently Nlayers is set to  %d.",options.Nlayer);
      nrerror(ErrStr);
    }
    if((!options.FULL_ENERGY && !options.FROZEN_SOIL) && options.Nlayer<1) {
      sprintf(ErrStr,"You must define at least 1 soil moisture layer to run the model.  Currently Nlayers is set to  %d.",options.Nlayer);
      nrerror(ErrStr);
    }
    if(options.IMPLICIT)  {
      if ( QUICK_FS ) 
        fprintf(stderr,"WARNING: IMPLICIT and QUICK_FS are both TRUE.\n\tThe QUICK_FS option is ignored when IMPLICIT=TRUE\n");
    }
    if( EXCESS_ICE ) {
      if ( !options.FULL_ENERGY )
        nrerror("set FULL_ENERGY = TRUE to run EXCESS_ICE option.");
      if ( !options.FROZEN_SOIL )
        nrerror("set FROZEN_SOIL = TRUE to run EXCESS_ICE option.");
      if ( options.QUICK_SOLVE ) {
        fprintf(stderr,"WARNING: QUICK_SOLVE and EXCESS_ICE are both TRUE.\n\tThis is an incompatible combination.  Setting QUICK_SOLVE to FALSE.\n");
        options.QUICK_SOLVE=FALSE;  
      }    
      if ( QUICK_FS ) 
        nrerror("QUICK_FS = TRUE and EXCESS_ICE = TRUE are incompatible options.");
    }
    if(options.Nlayer > MAX_LAYERS) {
      sprintf(ErrStr,"Global file wants more soil moisture layers (%d) than are defined by MAX_LAYERS (%d).  Edit user_def.h and recompile.",options.Nlayer,MAX_LAYERS);
      nrerror(ErrStr);
    }
    if(options.Nnode > MAX_NODES) {
      sprintf(ErrStr,"Global file wants more soil thermal nodes (%d) than are defined by MAX_NODES (%d).  Edit user_def.h and recompile.",options.Nnode,MAX_NODES);
      nrerror(ErrStr);
    }

    // Validate lake parameter information
    if (options.LAKES) {
      if (!options.FULL_ENERGY) {
        sprintf(ErrStr, "FULL_ENERGY must be TRUE if the lake model is to be run.");
        nrerror(ErrStr);
      }
      if ( strcmp ( names->lakeparam, "MISSING" ) == 0 )
        nrerror("\"LAKES\" was specified, but no lake parameter file has been defined.  Make sure that the global file defines the lake parameter file on the line that begins with \"LAKES\".");
      if (IS_INVALID(global_param.resolution) || global_param.resolution == 0) {
        sprintf(ErrStr, "The model grid cell resolution (RESOLUTION) must be defined in the global control file when the lake model is active.");
        nrerror(ErrStr);
      }
      if (global_param.resolution > 360 && !options.EQUAL_AREA) {
        sprintf(ErrStr, "For EQUAL_AREA=FALSE, the model grid cell resolution (RESOLUTION) must be set to the number of lat or lon degrees per grid cell.  This cannot exceed 360.");
        nrerror(ErrStr);
      }
      if (options.COMPUTE_TREELINE) {
        sprintf(ErrStr, "LAKES = TRUE and COMPUTE_TREELINE = TRUE are incompatible options.");
        nrerror(ErrStr);
      }
    }

    //validate and set Noutfiles based on set options
    options.Noutfiles = 2;  // This is the default for the OUTPUT_FORCE=FALSE case
    if (options.FROZEN_SOIL) {
      options.Noutfiles++;
    }
    if (options.PRT_SNOW_BAND) {
      options.Noutfiles++;
    }
    if (options.LAKES) {
      options.Noutfiles++;
    }

    /*********************************
    Output major options to stderr
    *********************************/
#if VERBOSE
    this->display_current_settings(DISP_ALL,names);
#else
    this->display_current_settings(DISP_VERSION,names);
#endif

#if VERBOSE
    fprintf(stderr,"Time Step = %d hour(s)\n",global_param.dt);
    fprintf(stderr,"Simulation start date = %02i/%02i/%04i\n",
            global_param.startday, global_param.startmonth, global_param.startyear);
    if (IS_VALID(global_param.nrecs) && global_param.nrecs > 0 )
      fprintf(stderr,"Number of Records = %d\n\n",global_param.nrecs);
    else 
      fprintf(stderr,"Simulation end date = %02i/%02i/%04i\n\n",
              global_param.endday, global_param.endmonth, global_param.endyear);
    fprintf(stderr,"Full Energy...................(%d)\n",options.FULL_ENERGY);
    fprintf(stderr,"Use Distributed Precipitation.(%d)\n",options.DIST_PRCP);
    if(options.DIST_PRCP)
      fprintf(stderr,"..Using Precipitation Exponent of %f\n",options.PREC_EXPT);
    fprintf(stderr,"Ground heat flux will be estimated ");
    if ( options.QUICK_FLUX ) 
      fprintf(stderr,"using Liang, Wood and Lettenmaier (1999).\n");
    else 
      fprintf(stderr,"using Cherkauer and Lettenmaier (1999).\n");
    fprintf(stderr,"Use Frozen Soil Model.........(%d)\n",options.FROZEN_SOIL);
    if( options.IMPLICIT ) 
      fprintf(stderr,".... Using the implicit solution for the soil heat equation.\n");
    else
      fprintf(stderr,".... Using the explicit solution for the soil heat equation.\n");
    if( options.EXP_TRANS )
      fprintf(stderr,".... Thermal nodes are exponentially distributed with depth.\n");
    else
      fprintf(stderr,".... Thermal nodes are linearly distributed with depth (except top two nodes).\n");
    if( EXCESS_ICE )
      fprintf(stderr,".... Excess ground ice is being considered.\n\t\tTherefore, ground ice (as a volumetric fraction) must be initialized for each\n\t\t   soil layer in the soil file.\n\t\tCAUTION: When excess ice melts, subsidence occurs.\n\t\t  Therefore, soil layer depths, damping depth, thermal node depths,\n\t\t     bulk densities, porosities, and other properties are now dynamic!\n\t\t  EXERCISE EXTREME CAUTION IN INTERPRETING MODEL OUTPUT.\n\t\t  It is recommended to add OUT_SOIL_DEPTH to your list of output variables.\n");
    if ( QUICK_FS ){
      fprintf(stderr,".... Using linearized UFWC curve with %d temperatures.\n", QUICK_FS_TEMPS);
    }
    fprintf(stderr,"Run Snow Model Using a Time Step of %d hours\n", 
            options.SNOW_STEP);
    fprintf(stderr,"Compress Output Files.........(%d)\n",options.COMPRESS);
    fprintf(stderr,"Correct Precipitation.........(%d)\n",options.CORRPREC);
    fprintf(stderr,"\n");
    fprintf(stderr,"Using %d Snow Bands\n",options.SNOW_BAND);
    fprintf(stderr,"Using %d Root Zones\n",options.ROOT_ZONES);
    if ( options.SAVE_STATE )
      fprintf(stderr,"Model state will be saved on = %02i/%02i/%04i\n\n",
              global_param.stateday, global_param.statemonth, global_param.stateyear);
    if ( options.OUTPUT_FORMAT == OutputFormat::BINARY_FORMAT ) {
      fprintf(stderr,"Model output is in standard BINARY format.\n");
    } else if ( options.OUTPUT_FORMAT == OutputFormat::ASCII_FORMAT){
      fprintf(stderr,"Model output is in standard ASCII format.\n");
    } else if (options.OUTPUT_FORMAT == OutputFormat::NETCDF_FORMAT) {
      fprintf(stderr, "Model output is in NETCDF format.\n");
    } else {
      fprintf(stderr, "Warning: model output is undefined (should be either BINARY, ASCII, or NETCDF).\n");
    }
    if ( LINK_DEBUG ) 
      fprintf(stderr,"Debugging code has been included in the executable.\n");
    else 
      fprintf(stderr,"Debugging code has not been compiled.\n");
#endif
  }

}
