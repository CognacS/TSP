#ifndef CPX_MODELS_H_  

#define CPX_MODELS_H_

#include "random.h"
#include <cut.h>

#include "mip_utility.h"
#include "graphs.h"
#include "graph_plot.h"

#define CC_EPSILON 0.1
#define CC_CUTOFF 2.0 - CC_EPSILON


// ************************** USEFUL STRUCTURES FOR PACKING **************************
typedef struct
{
	instance* inst;
	void* args;
	int (*sep_procedure)(void*, void*, instance*, void*, double*, int, int, int);
	int (*rej_procedure)(void*, void*, instance*, void*, double*, int, int, int);
} callback_instance;

typedef struct
{
	instance* inst;
	CPXCALLBACKCONTEXTptr context;
	double* value;
	int* index;
	int* labels;
} concorde_instance;


// ****************************** SYMMETRIC TSP MODELS *******************************
// ******** CALLBACKS *********
static int CPXPUBLIC sec_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle);

// *********** CUTS ***********
// add SE constraints on subtour from an available infeasable solution
// conform to all types of cuts (e.g. static, lazy, callback, etc...)
#define CUT_FLAGS_MINCUT_ENABLED 0b00000001
int add_sec_on_subtours(void* env, void* cbdata, instance* inst, void* args, double* xstar, int wherefrom, int purgeable, int flags);
int CC_add_sec_on_subtours(void* env, void* cbdata, instance* inst, void* args, double* xstar, int wherefrom, int purgeable, int flags);

// callback function to add a user cut in CPlex from Concorde
int doit_fn_concorde2cplex(double cutval, int cutcount, int* cut, void* args);

// ********** MODELS **********
void build_model_base_undirected(instance* inst, CPXENVptr env, CPXLPptr lp);

// ******** SOLUTIONS *********
void solve_benders(instance* inst, CPXENVptr env, CPXLPptr lp);
void solve_callback(instance* inst, CPXENVptr env, CPXLPptr lp);
void solve_symmetric_tsp(instance* inst, CPXENVptr env, CPXLPptr lp);
// ***********************************************************************************


// ****************************** ASYMMETRIC TSP MODELS ******************************
// *********** CUTS ***********
// add SE constraints on pairs of nodes
void add_sec2_asymmetric(instance* inst, CPXENVptr env, CPXLPptr lp);

// ********** MODELS **********
void build_model_base_directed(instance* inst, CPXENVptr env, CPXLPptr lp);
void build_model_mtz(instance* inst, CPXENVptr env, CPXLPptr lp);
void build_model_gg(instance* inst, CPXENVptr env, CPXLPptr lp);

// ******** SOLUTIONS *********
void solve_asymmetric_tsp(instance* inst, CPXENVptr env, CPXLPptr lp);
// ***********************************************************************************

#endif