#ifndef VEGCONDITIONS_H_
#define VEGCONDITIONS_H_

struct VegConditions {
  VegConditions();
  double snowFree;            // Previously index 0
  double canopyIfOverstory;   // Previously index 1
  double snowCovered;         // Previously index 2
  double glacierSurface;

  enum VegetationConditions { //TODO: rename shorter
    SNOW_FREE_CASE,
    CANOPY_IF_OVERSTORY_CASE,
    SNOW_COVERED_CASE,
    GLACIER_SURFACE_CASE,
    NUM_VEGETATION_CONDITIONS // End marker to automatically count the number of veg conditions.
  };
  // This operator overload allows array like indexing for this class for the VegetationConditions enum type.
  double& operator[](VegetationConditions type);
};


#endif /* VEGCONDITIONS_H_ */
