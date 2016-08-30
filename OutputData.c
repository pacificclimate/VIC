/*
 * OutputData.c
 *
 *  Created on: Aug 30, 2016
 *      Author: mfischer
 */

#include "OutputData.h"

OutputData::OutputData() {
	data = NULL;
	aggdata = NULL;
}

OutputData::~OutputData() {
	delete [] data;
	delete [] aggdata;
}


