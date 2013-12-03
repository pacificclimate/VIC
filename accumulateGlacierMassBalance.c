
#include "GlacierMassBalanceResult.h"
#include "vicNl.h"

void resetAccumulationValues(std::vector<HRU>* hruList) {
  for (std::vector<HRU>::iterator hru = hruList->begin(); hru != hruList->end(); ++hru) {
    if (hru->isGlacier) {
      hru->glacier.cum_mass_balance = 0;
    }
  }
}

void accumulateGlacierMassBalance(const dmy_struct* dmy, int rec, dist_prcp_struct* prcp, const ProgramState* state) {

  if (IS_INVALID(state->global_param.glacierAccumStartYear) || IS_INVALID(state->global_param.glacierAccumStartMonth)
      || IS_INVALID(state->global_param.glacierAccumStartDay) || IS_INVALID(state->global_param.glacierAccumInterval)) {
    return; // If these have not been set in the global file, then don't bother with accumulation.
  }

  if (rec + 1 >= state->global_param.nrecs) {
    return; // Reached the end of the model simulation, don't do anything to avoid dmy[rec + 1] being out of bounds.
  }

  const dmy_struct nextDate = dmy[rec + 1];
  // If the next time step is the start of the new accumulation interval then output results (pass to glacier model)
  // and reinitialize all the accumulation values.
  if ((abs(nextDate.year - state->global_param.glacierAccumStartYear) % state->global_param.glacierAccumInterval == 0)
      && nextDate.month == state->global_param.glacierAccumStartMonth && nextDate.day == state->global_param.glacierAccumStartDay) {

    GlacierMassBalanceResult result(prcp->hruList, dmy[rec]);
    result.printForDebug();
    //TODO: pass the result object to the glacier dynamics module somehow.

    resetAccumulationValues(&(prcp->hruList));

  } else {  // Not the final time step in the specified interval, accumulate mass balance.
    // Initialize on the first time step.
    if (rec == 0) {
      resetAccumulationValues(&prcp->hruList);
    }

    // Accumulate mass balance for each glacier hru.
    for (std::vector<HRU>::iterator hru = prcp->hruList.begin(); hru != prcp->hruList.end(); ++hru) {
      if (hru->isGlacier) {
        if (IS_VALID(hru->glacier.mass_balance)) {
          hru->glacier.cum_mass_balance += hru->glacier.mass_balance;
        }
      }
    }
  }

}
