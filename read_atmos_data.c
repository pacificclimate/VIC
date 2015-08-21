#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vicNl.h"
#include <netcdf.h>
#define __STDC_LIMIT_MACROS 1
#include <stdint.h>

static char vcid[] = "$Id$";

void read_atmos_data(FILE                 *infile,
                     int                   ncid,
                     int                   file_num,
                     int                   forceskip,
                     double              **forcing_data,
                     soil_con_struct      *soil_con,
                     const ProgramState   *state)
/**********************************************************************
  read_atmos_data
  
  This routine reads in atmospheric data values from a binary/ascii file.

  BINARY
  Binary data are always specified as unsigned or signed ints, and a
  multiplier is used to convert to float.

  atmos variable: type:            model units:
  
  precipitation   unsigned short   mm per file_dt
  temperature       signed short   C
  wind              signed short   m/s

  swap bytes code from Kernighan, Brian W. and Rob Pike, "The practice of
  programming", Addison-Wesley, Reading, Massachusetts, 1999, 267 pp,
  page 206.   		
  
  ASCII
  ASCII data should have the same units as given in the table above.

  
  Supported Input Field Combinations, options in parenthesis optional:
  
  Daily met data and daily model timestep - prec, tmax, tmin, (wind)
  Daily met data and subdaily model timestep - prec, tmax, tmin, (wind)
  
  If the code is modified check;
  * for BINARY, the number of fields is correct
  * get_global flags are implemented
  
  Modifications:
  01/10/00 Modified to read a generic Binary, or ASCII column 
           data file and read the contents into the provided  
           data arrays.                                    KAC
  ??-???-?? Replaced NF with global_param.dt in condition checking
	    whether forcing file contains enough records to cover
	    the time range of the simulation.	(documented by TJB)
  2006-08-23 Changed order of fread/fwrite statements from ...1, sizeof...
             to ...sizeof, 1,... 				GCT
  2006-Sep-23 Fixed RCS ID string.				TJB
  2007-Jan-15 Added PRT_HEADER option; now binary forcing files
	      might have headers, and these need to be skipped.	TJB
  2011-Nov-04 Fixed warning message dealing with insufficient
	      records.						TJB

  **********************************************************************/
{
  int rec;
  int skip_recs;
  int i;
  int endian;
  int fields;
  int Nfields;
  int day = 0;
  unsigned short ustmp;
  signed short stmp;
  char str[MAXSTRING + 1];
  char ErrStr[MAXSTRING + 1];
  unsigned short Identifier[4];
  int Nbytes;

  // get number of forcing variable types
  Nfields = state->param_set.N_TYPES[file_num];

  /** locate starting record **/
  /* if ascii then the following refers to the number of lines to skip,
   if binary the following needs multiplying by the number of input fields */
  skip_recs = (int) ((float) (state->global_param.dt * forceskip)) / (float) state->param_set.FORCE_DT[file_num];
  if ((((state->global_param.dt < 24
      && (state->param_set.FORCE_DT[file_num] * forceskip) % state->global_param.dt) > 0))
      || (state->global_param.dt == 24
          && (state->global_param.dt % state->param_set.FORCE_DT[file_num] > 0)))
    nrerror("Currently unable to handle a model starting date that does not correspond to a line in the forcing file.");

  /** Error checking - Model can be run at any time step using daily forcing
   data, but if sub-daily data is used, the model must be run at the
   same time step as the data.  That way aggregation and disaggragation
   techniques are left to the user. **/
  if (state->param_set.FORCE_DT[file_num] < 24
      && state->global_param.dt != state->param_set.FORCE_DT[file_num]) {
    sprintf(ErrStr,
        "When forcing the model with sub-daily data, the model must be run at the same time step as the forcing data.  Currently the model time step is %i hours, while forcing file %i has a time step of %i hours.",
        state->global_param.dt, file_num, state->param_set.FORCE_DT[file_num]);
    nrerror(ErrStr);
  }

  if (state->param_set.FORCE_FORMAT[file_num] == NETCDF) {
    /*****************************
     *  Read NetCDF Forcing Data  *
     *****************************/

    /* Assumptions, for now:
     * -dims are all (time, lat, lon)
     * -time, lat, lon are named accordingly
     * -vartype is double for dimvars, ... (?)
     * -need to implement proper FP comparison of dims
     */

    /* TODO some of these should be parameters once this is complete... */
    const int nforcesteps = state->global_param.nrecs * state->global_param.dt
        / state->param_set.FORCE_DT[file_num]; /* number of forcing timesteps to be loaded */

    int ndims; /* scratch for nc_inq_varndims */
    nc_type vartype; /* scratch for nc_inq_vartype */
    int ncerr;

    /* dims *//* rename ndays because it may be misleading */
    size_t ndays, nlats, nlons;
    int timevarid, latvarid, lonvarid, timetype, lattype, lontype, timedimid,
        latdimid, londimid;
    double *time, *lats, *lons; /* actual dimvar data *//* not using time, yet; just assuming... */
    size_t dimvarstart = 0, dimvarcount;

    /* vars */
    int varids[Nfields];
    int vardimids[3];

    /* hyperslab bounds and permutation */
    size_t starts[3] = { (size_t) skip_recs, SIZE_MAX, SIZE_MAX }; /* FIXME this may need to be multiplied by FORCE_DT or something?  not sure of the semantics here but for now we have 24h forcings so this shouldn't matter */
    size_t counts[3] = { (size_t) nforcesteps, 1, 1 }; /* TODO make this support multiple cells, and allow permutation to fit dim reordering... */
    int dimlengths[3] = { 1, 1, nforcesteps }; /* FIXME THIS SHOULD BE A PARAMETER */
    const ptrdiff_t perm[3] =
        { 1, dimlengths[1] * dimlengths[2], dimlengths[2] };

    /* handle dimvars */
    assert(nc_inq_varid(ncid, "time", &timevarid) == NC_NOERR);
    assert(nc_inq_vartype(ncid, timevarid, &vartype) == NC_NOERR);
    assert(vartype == NC_DOUBLE || vartype == NC_FLOAT); // supports legacy double and new disaggregated float types
    assert(nc_inq_varndims(ncid, timevarid, &ndims) == NC_NOERR);
    assert(ndims == 1);
    assert(nc_inq_vardimid(ncid, timevarid, &timedimid) == NC_NOERR);
    assert(nc_inq_dimlen(ncid, timedimid, &ndays) == NC_NOERR);
    /* use this later to determine whether we have enough data */

    assert(nc_inq_varid(ncid, "lat", &latvarid) == NC_NOERR);
    assert(nc_inq_vartype(ncid, latvarid, &vartype) == NC_NOERR);
    assert(vartype == NC_DOUBLE || vartype == NC_FLOAT); // supports legacy double and new disaggregated float types
    assert(nc_inq_varndims(ncid, latvarid, &ndims) == NC_NOERR);
    assert(ndims == 1);
    assert(nc_inq_vardimid(ncid, latvarid, &latdimid) == NC_NOERR);
    assert(nc_inq_dimlen(ncid, latdimid, &nlats) == NC_NOERR);

    assert(nc_inq_varid(ncid, "lon", &lonvarid) == NC_NOERR);
    assert(nc_inq_vartype(ncid, lonvarid, &vartype) == NC_NOERR);
    assert(vartype == NC_DOUBLE || vartype == NC_FLOAT); // supports legacy double and new disaggregated float types
    assert(nc_inq_varndims(ncid, lonvarid, &ndims) == NC_NOERR);
    assert(ndims == 1);
    assert(nc_inq_vardimid(ncid, lonvarid, &londimid) == NC_NOERR);
    assert(nc_inq_dimlen(ncid, londimid, &nlons) == NC_NOERR);


    /* get lat/lon *//* FIXME also get time and scan for it */
    lats = (double *) malloc(nlats * sizeof(*lats));
    assert(lats != NULL);
    dimvarcount = nlats;
    nc_get_vara_double(ncid, latvarid, &dimvarstart, &dimvarcount, lats);
    lons = (double *) malloc(nlons * sizeof(*lons));
    assert(lons != NULL);
    dimvarcount = nlons;
    nc_get_vara_double(ncid, lonvarid, &dimvarstart, &dimvarcount, lons);

    /* figure out slices */
    for (unsigned int latidx = 0; latidx < nlats; latidx++) {
      if (lats[latidx] == soil_con->lat) { /* FIXME, proper FP comparison would be really good!! */
        starts[1] = latidx;
        break;
      }
      assert(latidx != (nlats - 1));
      /* die if we make it to the last iteration without finding a match FIXME better spew */
    }
    for (unsigned int lonidx = 0; lonidx < nlons; lonidx++) {
      if (lons[lonidx] == soil_con->lng) {
        starts[2] = lonidx;
        break;
      }
      assert(lonidx != (nlons - 1));
    }

    /* handle vars */
    for (int varidx = 0; varidx < Nfields; ++varidx) {
      int attlen;
      float scale_factor = NAN, inverse_scale_factor = NAN; /* leave uninitialized in case of has_inverse_scale_factor */
      int has_inverse_scale_factor = 0;

      /* Get varid + check type */
      // Get variable id
      std::string variableKey = std::string(state->param_set.TYPE[state->param_set.FORCE_INDEX[file_num][varidx]].varname);
    	if (state->forcing_mapping.find(variableKey) == state->forcing_mapping.end()) {
    	        throw VICException("Error: could not find forcing variable in forcing_mapping: " + variableKey);
    	}
    	std::string varName = state->forcing_mapping.at(variableKey);
      fprintf(stderr,
          "Reading NetCDF forcing variable %s (NetCDF input variable name %s), slice [%d..%d,%d..%d,%d..%d] ... ",
          state->param_set.TYPE[state->param_set.FORCE_INDEX[file_num][varidx]].varname,
					varName.c_str(),
					(int) starts[0],
					(int) (starts[0] + counts[0] - 1), (int) starts[1],
          (int) (starts[1] + counts[1] - 1), (int) starts[2],
          (int) (starts[2] + counts[2] - 1));

      assert(nc_inq_varid(ncid, varName.c_str(), &varids[varidx]) == NC_NOERR);

      // Get variable type
      assert(nc_inq_vartype(ncid, varids[varidx], &vartype) == NC_NOERR);

      /* Check ndims + dim order */
      assert(nc_inq_varndims(ncid, varids[varidx], &ndims) == NC_NOERR);
      assert(ndims == 3);
      assert(nc_inq_vardimid(ncid, varids[varidx], vardimids) == NC_NOERR);
      assert((vardimids[0] == timedimid) && (vardimids[1] == latdimid) && (vardimids[2] == londimid));

      /* Get and convert data */
      switch (vartype) {
        case NC_SHORT: { // Legacy VIC short integer type input with scaling factors, for backward compatibility
          /* TODO check for relevant return code values instead of just NC_NOERR for cases where value might just not be present, although require at least one of scale_factor and inverse_scale_factor for integer-packed data */
          if (nc_get_att_float(ncid, varids[varidx], "inverse_scale_factor",
              &inverse_scale_factor) == NC_NOERR)	//TODO: move outside of switch
            has_inverse_scale_factor = 1;
          else
            assert(nc_get_att_float(ncid, varids[varidx], "scale_factor", &scale_factor) == NC_NOERR);
          short int *data = (short int *) malloc(
              (state->global_param.nrecs * state->global_param.dt) * sizeof(short));

          if ((ncerr = nc_get_varm_short(ncid, varids[varidx], starts, counts,
              NULL, perm, data)) != NC_NOERR) {
            fprintf(stderr, "Error reading NetCDF variable data: %s\n",
                nc_strerror(ncerr));
            exit(1);
          }
          fprintf(stderr, "done\n");
          /* FIXME (and below) handle missing value (probably by bailing, because this data should be contiguous) */
          if (has_inverse_scale_factor)
            /* Implemented for numerically-identical operation to classic VIC input */
            for (int rec = 0; rec < nforcesteps; rec++)
              forcing_data[state->param_set.FORCE_INDEX[file_num][varidx]][rec] = (double) data[rec]
                  / inverse_scale_factor;
          else
            for (int rec = 0; rec < nforcesteps; rec++)
              forcing_data[state->param_set.FORCE_INDEX[file_num][varidx]][rec] = (double) data[rec]
                  * scale_factor;
          free(data);
          break;
        }
        case NC_USHORT: {  // Legacy VIC unsigned integer type input with scaling factors, for backward compatibility
          if (nc_get_att_float(ncid, varids[varidx], "inverse_scale_factor",
              &inverse_scale_factor) == NC_NOERR)
            has_inverse_scale_factor = 1;
          else
            assert(
                nc_get_att_float(ncid, varids[varidx], "scale_factor", &scale_factor) == NC_NOERR);
          unsigned short int *data = (unsigned short int *) malloc(
              (state->global_param.nrecs * state->global_param.dt) * sizeof(short));
          if ((ncerr = nc_get_varm_ushort(ncid, varids[varidx], starts, counts,
              NULL, perm, data)) != NC_NOERR) {
            fprintf(stderr, "Error reading NetCDF forcing variable data: %s\n",
                nc_strerror(ncerr));
            exit(1);
          }
          if (has_inverse_scale_factor)
            /* Implemented for numerically-identical operation to classic VIC input */
            for (int rec = 0; rec < nforcesteps; rec++)
              forcing_data[state->param_set.FORCE_INDEX[file_num][varidx]][rec] = (double) data[rec] / inverse_scale_factor;
          else
            for (int rec = 0; rec < nforcesteps; rec++)
              forcing_data[state->param_set.FORCE_INDEX[file_num][varidx]][rec] = (double) data[rec] * scale_factor;
          free(data);
          break;
        }
        case NC_FLOAT: { // Supports new disaggregated forcing input types
        	float *data = (float *) malloc(
        	              (state->global_param.nrecs * state->global_param.dt) * sizeof(double));

					if ((ncerr = nc_get_varm_float(ncid, varids[varidx], starts, counts, NULL, perm, data)) != NC_NOERR) {
						fprintf(stderr, "Error reading NetCDF forcing variable data: %s\n",
								nc_strerror(ncerr));
						exit(1);
					}
          fprintf(stderr, "done\n");
        	for (int rec = 0; rec < nforcesteps; rec++)
        	  forcing_data[state->param_set.FORCE_INDEX[file_num][varidx]][rec] = (double) data[rec];
					free(data);
					break;
        }
        case NC_DOUBLE: { // Why not handle doubles too?
        	double *data = (double *) malloc(
        	              (state->global_param.nrecs * state->global_param.dt) * sizeof(double));

					if ((ncerr = nc_get_varm_double(ncid, varids[varidx], starts, counts, NULL, perm, data)) != NC_NOERR) {
						fprintf(stderr, "Error reading NetCDF forcing variable data: %s\n",
								nc_strerror(ncerr));
						exit(1);
					}
          fprintf(stderr, "done\n");
        	for (int rec = 0; rec < nforcesteps; rec++)
        	  forcing_data[state->param_set.FORCE_INDEX[file_num][varidx]][rec] = data[rec];
					free(data);
					break;
        }
        default:
        	fprintf(stderr,
        	          "Error reading NetCDF forcing variable %s. Type not supported.",
        	          varids[varidx], state->param_set.TYPE[state->param_set.FORCE_INDEX[file_num][varidx]].varname);
          assert(0);
          break;
      }

      rec = nforcesteps;
    }

    /* NC_MAX_NAME, NC_MAX_VAR_DIMS, nc_inq_varid(ncid, "name", &varid); (for expected dim) 
     *    nc_inq_var(ncid, varid, );  Don't need this until we're actually figuring out dims dynamically
     * nc_inq_dim(ncid, dimid, name_str, &dimlength); */

    /* TODO:
     * -position dim support?
     * -better matching of lats + lons
     * -support for other possible synonymous variable names (check that only ONE matches!)
     * -close files (inconsequential for now)
     */
  }

  /***************************
   Read BINARY Forcing Data
   ***************************/

  else if (state->param_set.FORCE_FORMAT[file_num] == BINARY) {

    /** test whether the machine is little-endian or big-endian **/
    i = 1;
    if (*(char *) &i == 1)
      endian = LITTLE;
    else
      endian = BIG;

    // Check for presence of a header, & skip over it if appropriate.
    // A VIC header will start with 4 instances of the identifier,
    // followed by number of bytes in the header (Nbytes).
    // Nbytes is assumed to be the byte offset at which the data records start.
    fseek(infile, 0, SEEK_SET);
    if (feof(infile))
      nrerror("No data in the forcing file.  Model stopping...");
    for (i = 0; i < 4; i++) {
      fread(&ustmp, sizeof(unsigned short), 1, infile);
      if (endian != state->param_set.FORCE_ENDIAN[file_num]) {
        ustmp = ((ustmp & 0xFF) << 8) | ((ustmp >> 8) & 0xFF);
      }
      Identifier[i] = ustmp;
    }
    if (Identifier[0] != 0xFFFF || Identifier[1] != 0xFFFF
        || Identifier[2] != 0xFFFF || Identifier[3] != 0xFFFF) {
      Nbytes = 0;
    } else {
      fread(&ustmp, sizeof(unsigned short), 1, infile);
      if (endian != state->param_set.FORCE_ENDIAN[file_num]) {
        ustmp = ((ustmp & 0xFF) << 8) | ((ustmp >> 8) & 0xFF);
      }
      Nbytes = (int) ustmp;
    }
    fseek(infile, Nbytes, SEEK_SET);

    /** if forcing file starts before the model simulation, 
     skip over its starting records **/
    fseek(infile, skip_recs * Nfields * sizeof(short), SEEK_CUR);
    if (feof(infile))
      nrerror(
          "No data for the specified time period in the forcing file.  Model stopping...");

    /** Read BINARY forcing data **/
    rec = 0;

    while (!feof(infile)
        && (rec * state->param_set.FORCE_DT[file_num]
            < state->global_param.nrecs * state->global_param.dt)) {

      for (i = 0; i < Nfields; i++) {
        if (state->param_set.TYPE[state->param_set.FORCE_INDEX[file_num][i]].SIGNED) {
          fread(&stmp, sizeof(short int), 1, infile);
          if (endian != state->param_set.FORCE_ENDIAN[file_num]) {
            stmp = ((stmp & 0xFF) << 8) | ((stmp >> 8) & 0xFF);
          }
          forcing_data[state->param_set.FORCE_INDEX[file_num][i]][rec] = (double) stmp / state->param_set.TYPE[state->param_set.FORCE_INDEX[file_num][i]].multiplier;
        } else {
          fread(&ustmp, sizeof(unsigned short int), 1, infile);
          if (endian != state->param_set.FORCE_ENDIAN[file_num]) {
            ustmp = ((ustmp & 0xFF) << 8) | ((ustmp >> 8) & 0xFF);
          }
          forcing_data[state->param_set.FORCE_INDEX[file_num][i]][rec] = (double) ustmp / state->param_set.TYPE[state->param_set.FORCE_INDEX[file_num][i]].multiplier;
        }
      }

      rec++;

    }
  }

  /**************************
   Read ASCII Forcing Data
   **************************/

  else {

    // No need to skip over a header here, since ascii file headers are skipped
    // in open_file().  However, if we wanted to read information from the header,
    // we'd want to do it here, after rewinding to the beginning of the file (or
    // moving the code that deals with headers from open_file() to this function
    // and to any other functions that read the files, so that those functions could
    // also read the headers if necessary).

    /* skip to the beginning of the required met data */
    for (i = 0; i < skip_recs; i++) {
      if (fgets(str, MAXSTRING, infile) == NULL)
        nrerror(
            "No data for the specified time period in the forcing file.  Model stopping...");
    }

    /* read forcing data */
    rec = 0;

    while (!feof(infile)
        && (rec * state->param_set.FORCE_DT[file_num]
            < state->global_param.nrecs * state->global_param.dt)) {
      for (i = 0; i < Nfields; i++)
        fscanf(infile, "%lf", &forcing_data[state->param_set.FORCE_INDEX[file_num][i]][rec]);
      fgets(str, MAXSTRING, infile);
      rec++;
    }
  }

  if (rec * state->param_set.FORCE_DT[file_num] < state->global_param.nrecs * state->global_param.dt) {
    sprintf(ErrStr, "Not enough records in forcing file %i (%i * %i = %i) to run the number of records defined in the global file (%i * %i = %i).  Check forcing file time step, and global file",
        file_num + 1, rec, state->param_set.FORCE_DT[file_num],
        rec * state->param_set.FORCE_DT[file_num], state->global_param.nrecs, state->global_param.dt,
        state->global_param.nrecs * state->global_param.dt);
    nrerror(ErrStr);
  }

}
