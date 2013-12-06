#include "VegConditions.h"
#include "vicNl.h"

VegConditions::VegConditions() : snowFree(INVALID), canopyIfOverstory(INVALID), snowCovered(INVALID), glacierSurface(INVALID) {
}

double& VegConditions::operator[](VegSurfType type) {
  switch (type) {
    case SNOW_FREE_CASE:
      return snowFree;
    case CANOPY_IF_OVERSTORY_CASE:
      return canopyIfOverstory;
    case SNOW_COVERED_CASE:
      return snowCovered;
    case GLACIER_SURFACE_CASE:
      return glacierSurface;
    case NUM_VEGETATION_CONDITIONS:
      throw VICException("Error: NUM_VEGETATION CONDITIONS is not an actual variable! (VegConditions::operator[]");
      break;
  }
  throw VICException("Error: undefined Vegetation type passed to VegConditions::operator[]");
}



