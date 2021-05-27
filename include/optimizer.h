#ifndef OPTIMIZER_H_  

#define OPTIMIZER_H_

#include "cpx_models.h"
#include "heuristics.h"

#define OPT_TL_PENALIZATION_MULT 10

typedef enum
{
	OPT_OK,
	OPT_TL_EXPIRED,
	OPT_HEUR_OK
} OptResult;

// optimization functions
OptResult TSPopt(Instance* inst);


#endif   /* OPTIMIZER_H_ */ 