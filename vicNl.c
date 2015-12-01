#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vicNl.h"
#include "global.h"
#include "StateIOContext.h"
#include <assert.h>
#include <omp.h>
#include <unistd.h>
#include <sstream>
#include <vector>

#include <chrono>
#include <ctime>

#include "WriteOutputAscii.h"
#include "WriteOutputBinary.h"
#include "WriteOutputNetCDF.h"
#include "WriteOutputAllCells.h"

static char vcid[] = "$Id: vicNl.c,v 5.14.2.19 2011/01/05 22:35:53 vicadmin Exp $";

void readSoilData(std::vector<cell_info_struct>& cell_data_structs,
    filep_struct filep, filenames_struct filenames,
    dmy_struct* dmy, ProgramState& state);

void runModel(std::vector<cell_info_struct>& cell_data_structs,
    filep_struct filep, filenames_struct filenames,
    out_data_file_struct* out_data_files_template, out_data_struct* out_data_list,
    dmy_struct* dmy, const ProgramState* state);

int initializeCell(cell_info_struct& cell,
    filep_struct filep, dmy_struct* dmy, filenames_struct filenames,
    const ProgramState* state);

int main(int argc, char *argv[])
/**********************************************************************
	vicNl.c		Dag Lohmann		January 1996

  This program controls file I/O and variable initialization as well as
  being the primary driver for the model.

  For details about variables, input files and subroutines check:
	http://ce.washington.edu/~hydro/Lettenmaier/Models/VIC/VIC_home.html

  UNITS: unless otherwise marked:
         all water balance components are in mm
	 all energy balance components are in mks
	 depths, and lengths are in m

  modifications:
  1997-98 Model was updated from simple 2 layer water balance to 
          an extension of the full energy and water balance 3 layer
	  model.                                                  KAC
  02-27-01 added controls for lake model                          KAC
  11-18-02 Updated storage of lake water for water balance 
           calculations.                                          LCB
  03-12-03 Modifed to add AboveTreeLine to soil_con_struct so that
           the model can make use of the computed treeline.     KAC
  04-10-03 Modified to initialize storm parameters using the state
           file.                                                KAC
  04-10-03 Modified to start the model by skipping records until the
           state file date is found.  This replaces the previous method
           of modifying the global file start date, which can change 
           the interpolation of atmospheric forcing data.        KAC
  04-15-03 Modified to store wet and dry fractions when intializing 
           water balance storage.  This accounts for changes in model
           state initialization, which now stores wet and dry fractions
           rather than just averagedvalues.                      KAC
  29-Oct-03 Modified the version display banner to print the version
	    string defined in global.h.					TJB
  01-Nov-04 Updated arglist for make_dist_prcp(), as part of fix for
	    QUICK_FLUX state file compatibility.			TJB
  02-Nov-04 Updated arglist for read_lakeparam(), as part of fix for
	    lake fraction readjustment.					TJB
  2005-Apr-13 OUTPUT_FORCE option now calls close_files().		TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
              out_data, out_data_files, and save_data structures.	TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      This included merging builtnames into filenames.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2006-Nov-07 Changed statefile to init_state in call to
	      check_state_file().					TJB
  2007-Jan-15 Added PRT_HEADER option; added call to
	      write_header().						TJB
  2007-Apr-04 Added option to continue run after a cell fails. 		GCT/KAC
  2007-Apr-21 Added calls to free_dmy(), free_out_data_files(),
	      free_out_data(), and free_veglib().  Added closing of
	      all parameter files.					TJB
  2007-Aug-21 Return ErrorFlag from initialize_model_state.		JCA
  2007-Sep-14 Excluded calls to free_veglib() and closing of parameter
	      files other than the soil param file for the case
	      when OUTPUT_FORCE=TRUE.					TJB
  2007-Nov-06 Moved computation of cell_area from read_lakeparam() to
	      read_soilparam() and read_soilparam_arc().		TJB
  2008-May-05 Added prcp fraction (mu) to initial water storage
	      computation.  This solves water balance errors for the
	      case where DIST_PRCP is TRUE.				TJB
  2009-Jan-16 Added soil_con.avgJulyAirTemp to argument list of
	      initialize_atmos().					TJB
  2009-Jun-09 Modified to use extension of veg_lib structure to contain
	      bare soil information.					TJB
  2009-Jul-07 Added soil_con.BandElev[] to read_snowband() arg list.	TJB
  2009-Jul-31 Replaced references to N+1st veg tile with references
	      to index of lake/wetland tile.				TJB
  2009-Sep-28 Replaced initial water/energy storage computations and
	      calls to calc_water_balance_error/calc_energy_balance_error
	      with an initial call to put_data.  Modified the call to
	      read_snowband().						TJB
  2009-Dec-11 Removed save_data structure from argument list of 
	      initialize_model_state().					TJB
  2010-Mar-31 Added cell_area to initialize_atmos().			TJB
  2010-Apr-28 Removed individual soil_con variables from argument list
	      of initialize_atmos() and replaced with *soil_con.	TJB
  2010-Nov-10 Added closing of state files.				TJB
  2011-Jan-04 Made read_soilparam_arc() a sub-function of
	      read_soilparam().						TJB
**********************************************************************/
{
  /** Read Model Options **/
  ProgramState state;
  state.initialize_global();

  filenames_struct filenames;
  cmd_proc(argc, argv, filenames.global, &state);


#if VERBOSE
  state.display_current_settings(DISP_VERSION, (filenames_struct*) NULL);
#endif
  /* Build input forcing variable name mappings */
  state.build_forcing_variable_mapping();
  /* Build default output variable name mappings */
  state.build_output_variable_mapping();
  /** Read Global Control File **/
  state.init_global_param(&filenames, filenames.global);
  /** Set up output data structures **/
  out_data_struct *out_data_list = create_output_list(&state);
  out_data_file_struct *out_data_files = set_output_defaults(out_data_list, &state);
  parse_output_info(filenames.global, out_data_files, out_data_list, &state);

  /** Check and Open Files **/
  filep_struct filep = get_files(&filenames, &state);

  if (!state.options.OUTPUT_FORCE) {
#if LINK_DEBUG
  state.open_debug();
#endif
  /** Read Vegetation Library File **/
  state.veg_lib = read_veglib(filep.veglib, &state.num_veg_types, state.options.LAI_SRC);
  }

  /** Make Date Data Structure **/
  dmy_struct* dmy = make_dmy(&state.global_param, &state);

  /** Initialize state **/
  std::vector<cell_info_struct> cell_data_structs; // Stores physical parameters for each grid cell
  readSoilData(cell_data_structs, filep, filenames, dmy, state);
  state.initGrid(cell_data_structs); // Calculate the grid cell parameters. This is used for NetCDF outputs.
  initializeNetCDFOutput(&filenames, out_data_files, out_data_list, &state);

  // Initialize state input/output if necessary.
  if (!state.options.OUTPUT_FORCE) {
    if (state.options.INIT_STATE)
      check_state_file(filenames.init_state, &state);
    /** open state file if model state is to be saved **/
    if (state.options.SAVE_STATE && strcmp(filenames.statefile, "NONE") != 0) {
      StateIOContext context(filenames.statefile, StateIO::Writer, &state);
      context.stream->initializeOutput();
    }
  }

  runModel(cell_data_structs, filep, filenames, out_data_files, out_data_list, dmy, &state);

  /** cleanup **/
  free_dmy(&dmy);
  delete [] out_data_files;
  free_out_data(&out_data_list);
  if (!state.options.OUTPUT_FORCE) {
    free_veglib(&state.veg_lib);
  }
  fclose(filep.soilparam);
  if (!state.options.OUTPUT_FORCE) {
    fclose(filep.vegparam);
    fclose(filep.veglib);
    if (state.options.SNOW_BAND > 1)
      fclose(filep.snowband);
    if (state.options.LAKES)
      fclose(filep.lakeparam);
  }

#if VERBOSE
  fprintf(stderr, "VIC exiting.\n");
#endif /* VERBOSE */

  return EXIT_SUCCESS;
} /* End Main Program */

//This method produces a warning based on the number of cells, and the specified RAM limit of the computer
//This warning is relevant only when running in image mode (where all cells should fit in memory).
//If the model is just being run sequentially, then only a single cell is in memory at a time so
//the amount of RAM available shouldn't be an issue (although it may still take a while).
void sanityCheckNumberOfCells(const int nCells, const ProgramState* state) {
  double GigsOfRam = state->options.MAX_MEMORY;
  const double approxBytesPerCell = 96000; //excluding the atmos forcing data
  double estimatedGigsOfRamUsed = approxBytesPerCell * nCells / (1024 * 1024 * 1024);
  if (GigsOfRam == 0.0) {
    fprintf(stderr, "Unlimited memory assumed.\n");
    return;
  }
  fprintf(stderr, "\nRAM limitation: %f Gb, estimated amount required: %f Gb, for %d cells\n", GigsOfRam, estimatedGigsOfRamUsed, nCells);
  if (estimatedGigsOfRamUsed > GigsOfRam) {
    fprintf(stderr, "Only continue if you know what you are doing, or you are not running in image mode.\n");
    fprintf(stderr, "Otherwise, consider running again using fewer cells.\nContinue anyways? [y/n] ");
    char c = getchar();
    if (c != 'y' && c != 'Y') {
      exit(0);
    }
  }
}

void readSoilData(std::vector<cell_info_struct>& cell_data_structs,
    filep_struct filep, filenames_struct filenames,
    dmy_struct* dmy, ProgramState& state) {

  /*****************************************
   * Read soil for all "active" grid cells *
   *****************************************/
  char done_reading_soil_file = FALSE, is_valid_soil_cell;
  int nallocatedcells = 10; /* arbitrary */
  int currentCellNumber = 0;

  double *lat = NULL;
  double *lng = NULL;
  int    *cellnum = NULL;

  soil_con_struct temp_soil_con;
  while (!done_reading_soil_file) {
    if (state.options.ARC_SOIL) {
      assert(0); // presently unsupported; maybe support vector functionality later?
      int  Ncells = 0;  // This will be initialized in read_soilparam_arc
      temp_soil_con = read_soilparam_arc(filep.soilparam, filenames.soil_dir, &Ncells, &is_valid_soil_cell, currentCellNumber, lat, lng, cellnum, &state);
      currentCellNumber++;
      if (currentCellNumber == Ncells)
        done_reading_soil_file = TRUE;
    } else {
      temp_soil_con = read_soilparam(filep.soilparam, filenames.soil_dir, &is_valid_soil_cell, &done_reading_soil_file, &state);
    }
    if (is_valid_soil_cell) {
      cell_info_struct currentCell;
      currentCell.soil_con = temp_soil_con;

      if (state.options.OUTPUT_FORMAT == OutputFormat::BINARY_FORMAT) {
      	currentCell.outputFormat = new WriteOutputBinary(&state);
  #if NETCDF_OUTPUT_AVAILABLE
      }
      else if (state.options.OUTPUT_FORMAT == OutputFormat::NETCDF_FORMAT) {
	#if PARALLEL
      	currentCell.outputFormat = new WriteOutputAllCells(&state);
	#else
      	currentCell.outputFormat = new WriteOutputNetCDF(&state);
  #endif
  #endif
      } else {
      	currentCell.outputFormat = new WriteOutputAscii(&state);
      }
      cell_data_structs.push_back(currentCell); // add an element to cell_data_structs vector
    }
  }
  sanityCheckNumberOfCells(cell_data_structs.size(), &state);
}

int initializeCell(cell_info_struct& cell,
    filep_struct filep, dmy_struct* dmy, filenames_struct filenames,
    const ProgramState* state) {

  const int Ndist = state->options.DIST_PRCP ? 2 : 1;

  #if LINK_DEBUG
      if (state->debug.PRT_SOIL)
        write_soilparam(&cell.soil_con, state);
  #endif

  #if QUICK_FS
      /** Allocate Unfrozen Water Content Table **/
      if(state->options.FROZEN_SOIL) {
        for(int i=0;i<MAX_LAYERS;i++) {
          cell.soil_con.ufwc_table_layer[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));
          for(int j=0;j<QUICK_FS_TEMPS+1;j++)
          cell.soil_con.ufwc_table_layer[i][j] = (double *)malloc(2*sizeof(double));
        }
        for(int i=0;i<MAX_NODES;i++) {
          cell.soil_con.ufwc_table_node[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));

          for(int j=0;j<QUICK_FS_TEMPS+1;j++)
          cell.soil_con.ufwc_table_node[i][j] = (double *)malloc(2*sizeof(double));
        }
      }
  #endif /* QUICK_FS */


  if (!state->options.OUTPUT_FORCE) {
    make_in_files(&filep, &filenames, &cell.soil_con, state);
    /** Read Grid Cell Vegetation Parameters **/
    read_vegparam(filep.vegparam, cell, state);
    calc_root_fractions(cell.prcp.hruList, &cell.soil_con, state);
#if LINK_DEBUG
    if (state->debug.PRT_VEGE) {
      write_vegparam(cell, state);
    }
#endif /* LINK_DEBUG*/
    if (state->options.LAKES) {
      cell.lake_con = read_lakeparam(filep.lakeparam, cell.soil_con, cell.prcp.hruList, state);
    }
  }
  else if (state->options.OUTPUT_FORCE) {
    make_in_files(&filep, &filenames, &cell.soil_con, state);
  }
  if (!state->options.OUTPUT_FORCE) {
      /** Read Elevation Band Data if Used **/
      read_snowband(filep.snowband, &cell.soil_con, state->options.SNOW_BAND);

  }
      /**************************************************
       Initialize Meteorological Forcing Values That
       Have not Been Specifically Set
       **************************************************/
  #if VERBOSE
      fprintf(stderr, "Initializing Forcing Data for cell at %4.5f %4.5f\n", cell.soil_con.lat, cell.soil_con.lng);
  #endif /* VERBOSE */

  /** allocate memory for the atmos_data_struct **/
  cell.atmos = alloc_atmos(state->global_param.nrecs, state->NR);
  initialize_atmos(cell.atmos, dmy, filep.forcing, filep.forcing_ncid, &cell.soil_con, state);

#if LINK_DEBUG
  if (state->debug.PRT_ATMOS)
    write_atmosdata(cell.atmos, state->global_param.nrecs, state);
#endif
  cell.writeDebug.initialize(cell.prcp.hruList.size(), state);
  /**************************************************
   Initialize Energy Balance and Snow Variables
   **************************************************/
  #if VERBOSE
      fprintf(stderr, "Model State Initialization\n"); //TODO: add information about which cell
  #endif /* VERBOSE */
  int ErrorFlag = initialize_model_state(&cell, dmy[0], filep, Ndist, filenames.init_state, state);

  if (ErrorFlag == ERROR) {
    if (state->options.CONTINUEONERROR == TRUE) {
      // Handle grid cell solution error
      fprintf(stderr,
          "Error initializing the model state (energy balance, water balance, and snow components) for cell %d (method initialize_model_state).  Cell has been marked as invalid and will be skipped for remainder of model run.\n",
          cell.soil_con.gridcel);
      return ERROR;
    } else {
      // Else exit program on cell solution error as in previous versions
      sprintf(cell.ErrStr,
          "Error initializing cell %d (method initialize_model_state).  Check your inputs before rerunning the simulation.  Exiting.\n",
          cell.soil_con.gridcel);
      vicerror(cell.ErrStr);
    }
  }
      return 0;
}

void printThreadInformation() {
  // Obtain thread number
  int tid = omp_get_thread_num();
  // Only master thread does this
  if (tid == 0) {
    printf("Thread %d getting environment info...\n", tid);
    // Get environment information
    int procs = omp_get_num_procs();
    int nthreads = omp_get_num_threads();
    int maxt = omp_get_max_threads();
    int inpar = omp_in_parallel();
    int dynamic = omp_get_dynamic();
    int nested = omp_get_nested();
    // Print environment information
    printf("Number of processors = %d\n", procs);
    printf("Number of threads = %d\n", nthreads);
    printf("Max threads = %d\n", maxt);
    printf("In parallel? = %d\n", inpar);
    printf("Dynamic threads enabled? = %d\n", dynamic);
    printf("Nested parallelism supported? = %d\n", nested);
  }
}

/************************************
 Run Model for all Active Grid Cells
 ************************************/
void runModel(std::vector<cell_info_struct>& cell_data_structs,
    filep_struct filep, filenames_struct filenames,
    out_data_file_struct* out_data_files_template, out_data_struct* out_data_list,
    dmy_struct* dmy, const ProgramState* state) {

	// Create vector for holding output data from one time iteration for all cells.
	std::vector<out_data_struct*> current_output_data;

	// Single WriteOutputAllCells object outputwriter takes care of writing one or all cells' data at a given time step
  WriteOutputAllCells *outputwriter = new WriteOutputAllCells(state);
  outputwriter->openFile();

#if PARALLEL_AVAILABLE
  unsigned int num_threads_allowed = std::max(atoi(std::getenv("OMP_NUM_THREADS")), 1);
	std::chrono::time_point<std::chrono::system_clock> start, end;
  if (num_threads_allowed > 1){
  	if (state->options.OUTPUT_FORCE) {
  		start = std::chrono::system_clock::now();
  	}
  }
#endif

  // Initializations
  for (unsigned int cellidx = 0; cellidx < cell_data_structs.size(); cellidx++) {

  	int initError = 0;
    initError = initializeCell(cell_data_structs[cellidx], filep, dmy, filenames, state);
    if (initError == ERROR) {
  	  cell_data_structs[cellidx].isValid = FALSE;
    }

    // Copy the format of the out_data_files_template and allocate in cell_data_structs[cellidx].outputFormat->dataFiles
    copy_data_file_format(out_data_files_template, cell_data_structs[cellidx].outputFormat->dataFiles, state);

    // Copy the format of the out_data_list (which is specific to this model run) and allocate space for this cell's element of current_output_data
    copy_output_data(current_output_data, out_data_list, state);

    /* Create output filename(s), and open (if already created, just open for appending).
        ASCII/binary output format will make two files per grid cell; NetCDF will make one file to rule them all */
    make_out_files(&filep, &filenames, &cell_data_structs[cellidx].soil_con, cell_data_structs[cellidx].outputFormat, state);

    // Write output file headers at initialization (does nothing in the NetCDF output format case)
    if (state->options.PRT_HEADER) {
      // Use out_data_list as a template for constructing the header of the output file (in ASCII mode)
//      out_data_struct out_data_template = copy_output_data(out_data_list, state);

    	// NOTE: commented out next 3 lines and replaced with vector-based implementation, as per vector implementation of copy_output_data.  NEED TO TEST THIS IN ASCII MODE
//      out_data_struct* out_data_template;
//      copy_output_data(out_data_template, out_data_list, state);
//      cell_data_structs[cellidx].outputFormat->write_header(out_data_template, dmy, state);
//
       std::vector<out_data_struct*> out_data_template;
       copy_output_data(out_data_template, out_data_list, state);
       cell_data_structs[cellidx].outputFormat->write_header(out_data_template[0], dmy, state);
    }
  } // for - grid cell loop

  /* If OUTPUT_FORCE is set to TRUE in the global parameters file then the full disaggregated
  forcing data array is written to file(s), and the full model run is skipped. */
  if (state->options.OUTPUT_FORCE) {
  	fprintf(stderr, "Writing forcing file...\n");
  	for (int rec = 0; rec < state->global_param.nrecs; rec++) {
#if PARALLEL_AVAILABLE
#pragma omp parallel for
#endif
      for (unsigned int cellidx = 0; cellidx < cell_data_structs.size(); cellidx++) {
    	write_forcing_file(&cell_data_structs[cellidx], rec, cell_data_structs[cellidx].outputFormat, current_output_data[cellidx], state, dmy);
      }
      outputwriter->write_data(current_output_data, out_data_files_template, &dmy[rec], state->global_param.out_dt, state);
			// Reset the aggdata for all variables
			for (int var_idx=0; var_idx<N_OUTVAR_TYPES; var_idx++) {
				for (unsigned int cell_idx = 0; cell_idx < cell_data_structs.size(); cell_idx++) {
					for (int elem=0; elem<out_data_list[var_idx].nelem; elem++) {
						current_output_data[cell_idx][var_idx].aggdata[elem] = 0;
					}
				}
			}
  	}
  }

  if (!state->options.OUTPUT_FORCE) {
#if VERBOSE
        fprintf(stderr, "Running Model...\n");
#endif /* VERBOSE */
#if PARALLEL_AVAILABLE
        if (num_threads_allowed > 1){
        	start = std::chrono::system_clock::now();
        }
#endif
  }

  /********************************************************
     Run Model for all Grid Cells, one Time Step at a time
  ********************************************************/
  for (int rec = 0; rec < state->global_param.nrecs; rec++) {

  	// Disaggregated meteorological forcings for entire time range is written in one shot for each grid cell
  	if (state->options.OUTPUT_FORCE) break;

#if PARALLEL_AVAILABLE
#pragma omp parallel for
#endif
    for (unsigned int cellidx = 0; cellidx < cell_data_structs.size(); cellidx++) {

//    	printThreadInformation();

    	// If this cell has been deemed invalid due to an error in an earlier time step, we don't process it.
    	if (cell_data_structs[cellidx].isValid == FALSE) continue;

    //TODO: These error files should not be global like this
    /** Update Error Handling Structure **/
    //state->Error.filep = filep;
    //state->Error.out_data_files = out_data_files;

      // Initialize storage terms on first time step
      if (rec == 0) {
        // Initialize the storage terms in the water and energy balances
        int putDataError = put_data(&cell_data_structs[cellidx], cell_data_structs[cellidx].outputFormat, current_output_data[cellidx], &dmy[0],
                		-state->global_param.nrecs, state);

        // Skip the rest of this cell if there is an error here.
        if (putDataError == ERROR) {
        	cell_data_structs[cellidx].isValid = FALSE;
          if (state->options.CONTINUEONERROR == TRUE) {
            fprintf(stderr, "Error initializing storage terms for cell %d (method put_data).  Cell has been marked as invalid and will be skipped for remainder of model run.\n", cell_data_structs[cellidx].soil_con.gridcel);
            continue;
          }
          else {
            sprintf(cell_data_structs[cellidx].ErrStr, "Error initializing storage terms for cell %d (method put_data).  Exiting.\n", cell_data_structs[cellidx].soil_con.gridcel);
            vicerror(cell_data_structs[cellidx].ErrStr);
          }
        }
      }

      int distPrecError = dist_prec(&cell_data_structs[cellidx], dmy, &filep, cell_data_structs[cellidx].outputFormat, current_output_data[cellidx], rec, FALSE, state);

      if (distPrecError == ERROR) {
      	cell_data_structs[cellidx].isValid = FALSE;
        if (state->options.CONTINUEONERROR == TRUE) {
          // Handle grid cell solution error
          fprintf(stderr,
              "Error processing cell %d (method dist_prec) at record (time step) %d.  Cell has been marked as invalid and will be skipped for remainder of model run.  An incomplete output file has been generated, check your inputs before re-running the simulation.\n",
              cell_data_structs[cellidx].soil_con.gridcel, rec);
        } else {
          // Else exit program on cell solution error as in previous versions
          sprintf(cell_data_structs[cellidx].ErrStr,
              "Error processing cell %d (method dist_prec) at record (time step) %d so the simulation has ended. Check your inputs before re-running the simulation.\n",
              cell_data_structs[cellidx].soil_con.gridcel, rec);
          vicerror(cell_data_structs[cellidx].ErrStr);
        }
      }

      // FIXME: should accumulateGlacierMassBalance have error checking?
      if (cell_data_structs[cellidx].isValid)
        accumulateGlacierMassBalance(&(cell_data_structs[cellidx].gmbEquation), dmy, rec, &(cell_data_structs[cellidx].prcp), &(cell_data_structs[cellidx].soil_con), state);

      /************************************
       Save model state at assigned date
       (after the final time step of the assigned date)
       ************************************/
      if (state->options.SAVE_STATE == TRUE
          && (dmy[rec].year == state->global_param.stateyear
          && dmy[rec].month == state->global_param.statemonth
          && dmy[rec].day == state->global_param.stateday
          && (rec + 1 == state->global_param.nrecs
          || dmy[rec + 1].day != state->global_param.stateday)))
      {
        write_model_state(&cell_data_structs[cellidx], filenames.statefile, state);
      }

#if QUICK_FS
      if(options.FROZEN_SOIL) {
        for(int i=0;i<MAX_LAYERS;i++) {
          for(int j=0;j<6;j++)
          free((char *)cell_data_structs[cellidx].soil_con.ufwc_table_layer[i][j]);
          free((char *)cell_data_structs[cellidx].soil_con.ufwc_table_layer[i]);
        }
        for(int i=0;i<MAX_NODES;i++) {
          for(int j=0;j<6;j++)
          free((char *)cell_data_structs[cellidx].soil_con.ufwc_table_node[i][j]);
          free((char *)cell_data_structs[cellidx].soil_con.ufwc_table_node[i]);
        }
      }
#endif /* QUICK_FS */
    } // for - grid cell loop

    // Write output data for all cells on this time iteration to file
    if(rec >= state->global_param.skipyear) {
    	outputwriter->write_data(current_output_data, out_data_files_template, &dmy[rec], state->global_param.out_dt, state);

      // Reset the aggdata for all variables (even those not necessarily being written, as some variables' aggdata values are derived from other variables)
    	for (unsigned int var_idx=0; var_idx<N_OUTVAR_TYPES; var_idx++) {
    		for (unsigned int cell_idx = 0; cell_idx < current_output_data.size(); cell_idx++) {
    			for (unsigned int elem=0; elem<out_data_list[var_idx].nelem; elem++) {
    				current_output_data[cell_idx][var_idx].aggdata[elem] = 0;
    			}
    		  // Reset the step count
    			cell_data_structs[cell_idx].fallBackStats.step_count = 0;
    		}
    	}
    }
//  	fprintf(stderr,".");
  } // for - time loop

	delete outputwriter;

#if PARALLEL_AVAILABLE
  //    	outputwriter->closeFile();

	if (num_threads_allowed > 1){
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		std::time_t end_time = std::chrono::system_clock::to_time_t(end);
		std::cout << "elapsed time in parallel loop: " << elapsed_seconds.count() << " seconds\n";
	}
#endif

  // Close NetCDF forcing file
  if(state->param_set.FORCE_FORMAT[0] == NETCDF)
    close_files(&filep, &filenames, state->options.COMPRESS, state);

  // Free up cell_data_structs
  for (unsigned int cellidx = 0; cellidx < cell_data_structs.size(); cellidx++) {
    cell_data_structs[cellidx].writeDebug.cleanup(cell_data_structs[cellidx].prcp.hruList.size(), state);
    free_atmos(state->global_param.nrecs, &cell_data_structs[cellidx].atmos);
    free_vegcon(cell_data_structs[cellidx]);
    free(cell_data_structs[cellidx].soil_con.AreaFract);
    free(cell_data_structs[cellidx].soil_con.BandElev);
    free(cell_data_structs[cellidx].soil_con.Tfactor);
    free(cell_data_structs[cellidx].soil_con.Pfactor);
    free(cell_data_structs[cellidx].soil_con.AboveTreeLine);
    delete cell_data_structs[cellidx].outputFormat;
  }

}  // runModel
