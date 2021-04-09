#ifndef MIP_UTILITY_H_  

#define MIP_UTILITY_H_

#include <cplex.h>

#include "chrono.h"
#include "log.h"
#include "error.h"
#include "tsp_data.h"

#define CUT_STATIC			-1
#define CUT_LAZY			-2
#define	CUT_USER			-3
#define CUT_CALLBACK_LOCAL	-4
#define CUT_CALLBACK_REJECT	-5

/**
* Unified driver to add constraints to the LP.
* */
void mip_add_cut(void* env, void* cbdata,
	int wherefrom, int nnz, double rhs, char sense,
	int* index, double* value, char* cname, int purgeable);

void mip_timelimit(CPXENVptr env, double timelimit, instance* inst);

#endif