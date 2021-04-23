#ifndef TSP_H_  

#define TSP_H_

#include "tsp_data.h"
#include "log.h"
#include "error.h"
#include "tsp_utility.h"
#include "cpx_models.h"
#include "chrono.h"

#define OPT_TL_PENALIZATION_MULT 10

typedef enum
{
	OPT_OK,
	OPT_TL_EXPIRED
} opt_result;

// optimization functions
int TSPopt(instance* inst, char getsol);


#endif   /* TSP_H_ */ 