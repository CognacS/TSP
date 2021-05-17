#ifndef OPTIMIZER_DATA_H_

#define OPTIMIZER_DATA_H_

#include "tsp_data.h"

typedef enum
{
	OPT_OK,
	OPT_TL_EXPIRED,
	OPT_HEUR_OK

} opt_result;

typedef struct
{
	CPXENVptr env;
	CPXLPptr lp;
} CplexData;

typedef struct
{
	instance* inst;
	CplexData* cpx;
	LinkedList* perflog;
} OptData;

inline int using_cplex(OptData* optdata) { return optdata->cpx != NULL; }


#endif