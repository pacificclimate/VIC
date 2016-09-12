#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"
#include <string.h>
#include <string>
#include <sstream>

static char vcid[] = "$Id$";

void ttrim( char *string );

HRU initHRU(veg_con_struct& veg, const ProgramState* state) {
  HRU hru;
  hru.energy.frozen = FALSE;
  hru.isGlacier = (veg.vegClass == state->options.GLACIER_ID);
  hru.isArtificialBareSoil = false;
  hru.veg_con = veg;
  hru.mu = 1;
  return hru;
}

void ProgramState::update_max_num_HRUs(int numHRUs){
	max_num_HRUs = numHRUs;
}

// For any type T, get the next value from the provided stream.
template<class T> void getValueFromStream(std::stringstream& stream, T& value) {
  if (!(stream >> value)) {
    throw VICException("Error in read_vegparam, not enough inputs on line containing " + stream.str());
  }
}

int getVegIndex(int vegClass, const ProgramState* state) {
// Identify current vegetation class
  int veg_index = INVALID_INT;
  for (int j = 0; j < state->num_veg_types; j++) {
    if (vegClass == state->veg_lib[j].veg_class) {
      veg_index = j;
      break;
    }
  }
  if (IS_INVALID(veg_index)) {
    std::stringstream ss;
    ss << "The vegetation class id " << vegClass
        << " defined for above-treeline is not defined in the vegetation library file.";
    throw VICException(ss.str());
  }
  return veg_index;
}

// MDF: changed return type to int so we can return numHRUs for main program to calculate the max # of HRUs across all cells,
// used to allocate space for in the state file
int read_vegparam(FILE *vegparam,
                   cell_info_struct& cell,
                   const ProgramState* state)

//void read_vegparam(FILE *vegparam,
//                   cell_info_struct& cell,
//                   const ProgramState* state)
/**********************************************************************
  read_vegparam.c    Keith Cherkauer and Dag Lohmann       1997

  This routine reads in vegetation parameters for the current grid cell.
  It also relates each vegetation class in the cell to the appropriate
  parameters in the vegetation library.

  Modifications:
  09-24-98  Modified to read root zone distribution information so
           that soil layer root fractions can be computed for new 
	   soil layer depths - see calc_root_fractions.c           KAC
  07-15-99 Modified to read LAI values from a new line in the vegetation
           parameter file.  Added specifically to work with the new
	   global LAI files.
  11-18-02 Added code to read in blowing snow parameters.          LCB
  03-27-03 Modified code to update Wdmax based on LAI values read in
           for the current grid cell.  If LAI is not obtained from this
           function, then the values cacluated in read_veglib.c are
           left unchanged.						DP & KAC
  2006-Nov-07 Allocates MaxVeg+1 veg tiles.				TJB
  2007-May-11 Changed some 'fscanf' statements to 'fgets' and 'sscanf' 
	      to count rootzone and BLOWING fields. Also tests for
	      fetch < 1.						GCT
  2007-Oct-31 Added missing brackets in if(options.GLOBAL_LAI) block.	TJB
  2008-Oct-23 Added blocks to free vegarr[].				LCB via TJB
  2009-Jan-16 Added logic for COMPUTE_TREELINE option.			TJB
  2009-Jun-09 Modified to use extension of veg_lib structure to contain
	      bare soil information.					TJB
  2009-Jun-17 Modified to understand both tabs and spaces as delimiters.TJB
  2009-Jun-17 Fixed incorrect placement of free vegarr[] for case of
	      GLOBAL_LAI==FALSE.					TJB
  2009-Jul-26 Allocate extra veg tile for COMPUTE_TREELINE case.	TJB
  2009-Jul-31 Removed extra veg tile for lake/wetland case.		TJB
  2009-Sep-14 Made error messages clearer.				TJB
  2009-Oct-01 Added error message for case of LAI==0 and overstory==1.	TJB
  2010-Apr-28 Replaced GLOBAL_LAI with VEGPARAM_LAI and LAI_SRC.	TJB
**********************************************************************/
{
  int             vegcel, numHRUs, skip;
  int             NoOverstory;
  char            str[500];
  char            ErrStr[MAXSTRING];
  char            line[MAXSTRING];
  char            tmpline[MAXSTRING];
  const char      delimiters[] = " \t";
  char            *token;
  size_t	  length;

  if(state->options.VEGPARAM_LAI) skip=2;
  else skip=1;

  NoOverstory = 0;

#if !NO_REWIND
  rewind(vegparam);
#endif  

  while ( ( fscanf(vegparam, "%d %d", &vegcel, &numHRUs) == 2 ) && vegcel != cell.soil_con.gridcel ){
    if (numHRUs < 0) {
      sprintf(ErrStr,"ERROR number of vegetation tiles (%i) given for cell %i is < 0.\n",numHRUs,vegcel);
      nrerror(ErrStr);
    }
    for (int i = 0; i <= numHRUs * skip; i++){
      if ( fgets(str, 500, vegparam) == NULL ){
        sprintf(ErrStr,"ERROR unexpected EOF for cell %i while reading root zones and LAI\n",vegcel);
        nrerror(ErrStr);
      }
    }
  }
  fgets(str, 500, vegparam); // read newline at end of veg class line to advance to next line
  if (vegcel != cell.soil_con.gridcel) {
    fprintf(stderr, "Error in vegetation file.  Grid cell %d not found\n", cell.soil_con.gridcel);
    exit(99);
  }

  cell.Cv_sum = 0.0;

  for (int i = 0; i < numHRUs; i++) {
    veg_con_struct curVeg;
    curVeg.zone_depth = (float*)calloc(state->options.ROOT_ZONES,sizeof(float));
    curVeg.zone_fract = (float*)calloc(state->options.ROOT_ZONES,sizeof(float));

    // Read the root zones line
    if ( fgets( line, MAXSTRING, vegparam ) == NULL ){
      sprintf(ErrStr,"ERROR unexpected EOF for cell %i while reading HRU number %d\n",vegcel,i);
      nrerror(ErrStr);
    }
    strcpy(tmpline, line);
    ttrim( tmpline );
    std::stringstream stream(tmpline);

    getValueFromStream(stream, curVeg.vegClass);
    getValueFromStream(stream, curVeg.Cv);

    float depth_sum = 0;
    float sum = 0.;
    for (int zone = 0; zone < state->options.ROOT_ZONES; zone++) {
      getValueFromStream(stream, curVeg.zone_depth[zone]);
      getValueFromStream(stream, curVeg.zone_fract[zone]);
      depth_sum += curVeg.zone_depth[zone];
      sum += curVeg.zone_fract[zone];
    }
    if(depth_sum <= 0) {
      throw VICException("Root zone depths must sum to a value greater than 0.");
    }
    if (sum != 1.) {
      fprintf(stderr, "WARNING: Root zone fractions sum to more than 1 ( = %f), normalizing fractions.  If the sum is large, check that your vegetation parameter file is in the form - <zone 1 depth> <zone 1 fract> <zone 2 depth> <zone 2 fract> ...\n", sum);
      for (int j = 0; j < state->options.ROOT_ZONES; j++) {
        curVeg.zone_fract[j] /= sum;
      }
    }

    if (state->options.BLOWING) {
      getValueFromStream(stream, curVeg.sigma_slope);
      getValueFromStream(stream, curVeg.lag_one);
      getValueFromStream(stream, curVeg.fetch);
      if( curVeg.sigma_slope <= 0. || curVeg.lag_one <= 0.) {
        sprintf(str,"Deviation of terrain slope must be greater than 0.");
        nrerror(str);
      }
      if( curVeg.fetch < 1.0  ) {
        sprintf(str,"ERROR - BLOWING parameter fetch should be >> 1 but cell %i has fetch = %.2f\n", cell.soil_con.gridcel, curVeg.fetch );
        nrerror(str);
      }
    }

    int curBandIndex = INVALID_INT;
    getValueFromStream(stream,curBandIndex);

    int NfieldsMax = 2 + 2 * state->options.ROOT_ZONES;  /* Number of expected fields this line */
    if( state->options.BLOWING ){
      NfieldsMax += 3;
    }

    curVeg.LAKE = 0;

    curVeg.vegIndex = getVegIndex(curVeg.vegClass, state);

    cell.Cv_sum += curVeg.Cv;

    if ( state->options.VEGPARAM_LAI && state->options.LAI_SRC == LAI_FROM_VEGPARAM) {
      // Read the LAI line
      if (fgets( line, MAXSTRING, vegparam ) == NULL) {
        sprintf(ErrStr,"ERROR unexpected EOF for cell %i while reading LAI for HRU number %d\n",vegcel,i);
        nrerror(ErrStr);
      }
      strcpy(tmpline, line);
      ttrim( tmpline );
      std::stringstream laiStream(tmpline);

      for (int j = 0; j < 12; j++) {
        //TODO: it is wrong to change program state here, should LAI be moved to the HRU class?
        try {
          getValueFromStream(laiStream, state->veg_lib[curVeg.vegIndex].LAI[j]);
        } catch (std::exception& e) {
          fprintf(stderr, "Error reading LAI values at gridcell %d for month %d\n", cell.soil_con.gridcel, j);
          throw;
        }
        if (state->veg_lib[curVeg.vegIndex].overstory && state->veg_lib[curVeg.vegIndex].LAI[j] == 0) {
          std::stringstream ss("Error: cell ");
          ss << cell.soil_con.gridcel << ", veg type " << curVeg.vegClass << " the specified veg class is listed as an overstory class in the veg LIBRARY, but the LAI given in the veg PARAM FILE for this tile for month " << j << " is 0.";
          throw VICException(ss.str());
        }
        state->veg_lib[curVeg.vegIndex].Wdmax[j] = LAI_WATER_FACTOR * state->veg_lib[curVeg.vegIndex].LAI[j];
      }
    }

    // Determine if cell contains non-overstory vegetation
    if (state->options.COMPUTE_TREELINE && !state->veg_lib[curVeg.vegIndex].overstory )
      NoOverstory++;

    //Create the HRU and add it to the vector.
    HRU e = initHRU(curVeg, state);
    e.bandIndex = curBandIndex;
    cell.prcp.hruList.push_back(e);
  } // end of loop

  // Determine if we have bare soil
  if(cell.Cv_sum>1.0){
    fprintf(stderr,"WARNING: Cv exceeds 1.0 at grid cell %d, fractions being adjusted to equal 1\n", cell.soil_con.gridcel);
    for(std::vector<HRU>::iterator hru = cell.prcp.hruList.begin(); hru != cell.prcp.hruList.end(); ++hru) {
      hru->veg_con.Cv = hru->veg_con.Cv / cell.Cv_sum;
    }
    cell.Cv_sum = 1.;
  }
  else if(cell.Cv_sum>0.99 && cell.Cv_sum<1.0){
    fprintf(stderr,"WARNING: Cv > 0.99 and Cv < 1.0 at grid cell %d, model assuming that bare soil is not to be run - fractions being adjusted to equal 1\n", cell.soil_con.gridcel);
    for(std::vector<HRU>::iterator hru = cell.prcp.hruList.begin(); hru != cell.prcp.hruList.end(); ++hru) {
      hru->veg_con.Cv = hru->veg_con.Cv / cell.Cv_sum;
    }
    cell.Cv_sum = 1.;
  }

  // Handle veg above the treeline
  if ( state->options.SNOW_BAND > 1 && state->options.COMPUTE_TREELINE
       && ( !NoOverstory && cell.Cv_sum == 1. ) ) {

    // All vegetation in the current cell is defined with overstory.
    // Add default non-overstory vegetation so that snow bands above treeline
    // can be sucessfully simulated.

    if ( state->options.AboveTreelineVeg < 0 ) {

      // Above treeline snowband should be treated as bare soil
      for(std::vector<HRU>::iterator hru = cell.prcp.hruList.begin(); hru != cell.prcp.hruList.end(); ++hru) {
        hru->veg_con.Cv -= ( 0.001 / (float)numHRUs );
      }
      cell.Cv_sum -= 0.001;

    }
    else {

      // Above treeline snowband should use the defined vegetation
      // add vegetation to typenum
      // check that veg type exists in library and does not have overstory
      veg_con_struct treeVeg;

      if(numHRUs > 0) {

        for(std::vector<HRU>::iterator hru = cell.prcp.hruList.begin(); hru != cell.prcp.hruList.end(); ++hru) {
          hru->veg_con.Cv -= ( 0.001 / (float)numHRUs );
        }


        treeVeg.Cv         = 0.001;
        treeVeg.vegClass  = state->options.AboveTreelineVeg;
        treeVeg.zone_depth = (float*)calloc( state->options.ROOT_ZONES,sizeof(float));
        treeVeg.zone_fract = (float*)calloc( state->options.ROOT_ZONES,sizeof(float));

        // Since root zones are not defined they are copied from the last
        // vegetation type.
        veg_con_struct* lastVeg = &cell.prcp.hruList[cell.prcp.hruList.size() - 1].veg_con;
        for (int j = 0; j < state->options.ROOT_ZONES; j++ ) {
          treeVeg.zone_depth[j] = lastVeg->zone_depth[j];
          treeVeg.zone_fract[j] = lastVeg->zone_fract[j];
        }
      }

      treeVeg.vegIndex = getVegIndex(treeVeg.vegClass, state);

      HRU treeHru = initHRU(treeVeg, state);
      treeHru.bandIndex = state->options.SNOW_BAND - 1; //TODO: currently set to the top snow band. Should a new HRU be added to all bands instead?
      cell.prcp.hruList.push_back(treeHru);

      if (state->veg_lib[treeVeg.vegIndex].overstory ) {
        std::stringstream ss;
        ss << "Vegetation class " << state->veg_lib[treeVeg.vegIndex].veg_class << "is defined to have overstory, so it cannot be used as the default vegetation type for above canopy snow bands.";
        throw VICException(ss.str());
      }
    }
  }

  // Bare soil tile
  if (cell.Cv_sum < 1.) {
    double CvPerBand = (1.0 - cell.Cv_sum) / (double)state->options.SNOW_BAND;
    // A bare soil HRU is added to each elevation.
    for (int band = 0; band < state->options.SNOW_BAND; band++) {
      veg_con_struct bareSoilVeg;
      bareSoilVeg.vegClass = state->num_veg_types; // Create a veg_class ID for bare soil, which is not mentioned in the veg library
      bareSoilVeg.vegIndex = state->num_veg_types;
      bareSoilVeg.Cv = CvPerBand;
      // Don't allocate any root-zone-related arrays
      if(state->options.BLOWING) {
        if (numHRUs > 0) {
          bareSoilVeg.sigma_slope = cell.prcp.hruList[0].veg_con.sigma_slope;
          bareSoilVeg.lag_one = cell.prcp.hruList[0].veg_con.lag_one;
          bareSoilVeg.fetch = cell.prcp.hruList[0].veg_con.fetch;
        }
        else {
          bareSoilVeg.sigma_slope = 0.005;
          bareSoilVeg.lag_one = 0.95;
          bareSoilVeg.fetch = 2000;
        }
      }
      HRU bareSoil = initHRU(bareSoilVeg, state);
      bareSoil.bandIndex = band;
      bareSoil.isArtificialBareSoil = true;
      cell.prcp.hruList.push_back(bareSoil);
    }
    cell.Cv_sum = 1;
  }
// MDF: added return of numHRUs for state file dimension allocation
  return numHRUs;
} 


/* trim trailing newlines */

#define END '\0'
#define NEW '\n'

void ttrim( char *c ) 
{
  while( (*c++ != END) );
    --c;
  for( ; *--c == NEW; *c = END );

}

