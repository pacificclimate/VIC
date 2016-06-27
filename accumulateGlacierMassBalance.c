
#include "GlacierMassBalanceResult.h"
#include "vicNl.h"

void resetAccumulationValues(std::vector<HRU>* hruList) {
	for (std::vector<HRU>::iterator hru = hruList->begin(); hru != hruList->end(); ++hru) {
		if (hru->isGlacier) {
			hru->glacier.cum_mass_balance = 0;
		}
	}
}

void accumulateGlacierMassBalance(GraphingEquation* gmbEquation, const dmy_struct* dmy, int rec, dist_prcp_struct* prcp, const soil_con_struct* soil, ProgramState* state) {

	if (IS_INVALID(state->global_param.glacierAccumStartYear) || IS_INVALID(state->global_param.glacierAccumStartMonth)
			|| IS_INVALID(state->global_param.glacierAccumStartDay) || IS_INVALID(state->global_param.glacierAccumInterval)) {
		return; // If these have not been set in the global file, then don't bother with accumulation.
	}

	if (rec == state->global_param.nrecs) {
		return; // Reached the end of the model simulation, don't do anything to avoid dmy[rec + 1] being out of bounds.
	}

	// Initialize on the first time step.
	if (rec == 0) {
		resetAccumulationValues(&prcp->hruList);
	}

	// Check if we have reached the glacier accumulation start point of the simulation
	if ( (dmy[rec].year == state->global_param.glacierAccumStartYear)
			&& (dmy[rec].month == state->global_param.glacierAccumStartMonth)
			&& (dmy[rec].day == state->global_param.glacierAccumStartDay) ) {

		state->glacier_accum_started = true;
	}

	// If we have reached the glacier accumulation start point of the simulation, accumulate mass balance for each glacier hru.
	if (state->glacier_accum_started) {
		for (std::vector<HRU>::iterator hru = prcp->hruList.begin(); hru != prcp->hruList.end(); ++hru) {
			if (hru->isGlacier) {
				if (IS_VALID(hru->glacier.mass_balance)) {
					hru->glacier.cum_mass_balance += hru->glacier.mass_balance;
				}
			}
		}
	}

	// nextDate is exactly one day ahead of the current time step
	const dmy_struct nextDate = dmy[rec + 1];

	// If the next time step is the start of the next accumulation interval then output results (to pass to glacier model)
	// and reinitialize all the accumulation values.
	if ( (nextDate.year > state->global_param.glacierAccumStartYear)
			&& (abs(nextDate.year - state->global_param.glacierAccumStartYear) % state->global_param.glacierAccumInterval == 0)
			&& nextDate.month == state->global_param.glacierAccumStartMonth
			&& nextDate.day == state->global_param.glacierAccumStartDay
			&& (((state->global_param.dt <= 12) && (dmy[rec].hour == 23)) || (state->global_param.dt == 24 )) ) {

		GlacierMassBalanceResult result(prcp->hruList, soil, dmy[rec]);
#if LINK_DEBUG
		fprintf(stderr, "accumulateGlacierMassBalance for cell at %4.5f %4.5f:\n", soil->lat, soil->lng );
		result.printForDebug();
#endif /* LINK_DEBUG */
		resetAccumulationValues(&(prcp->hruList));
		*gmbEquation = result.equation; // update GMB polynomial for this cell
	}
}
