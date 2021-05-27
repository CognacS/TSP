#ifndef CPX_MODELS_H_  

#define CPX_MODELS_H_

#include "plot.h"
#include "cpx_utility.h"
#include "opt_utility.h"

#define CC_EPSILON 0.1
#define CC_CUTOFF 2.0 - CC_EPSILON


// ************************** USEFUL STRUCTURES FOR PACKING **************************
typedef struct
{
	Instance* inst;
	void* args;
	int (*sep_procedure)(void*, void*, Instance*, void*, double*, int, int, int);
	int (*rej_procedure)(void*, void*, Instance*, void*, double*, int, int, int);
	double (*prob_function)(int, double);
	double prob_decay;
} CallbackInstance;

typedef struct
{
	Instance* inst;
	CPXCALLBACKCONTEXTptr context;
	double* value;
	int* index;
	int* labels;
	int local;
} ConcordeInstance;

// ****************************** SYMMETRIC TSP MODELS *******************************
// ******** CALLBACKS *********
// callback separation function
int CPXPUBLIC sec_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle);

// *********** CUTS ***********
// add SE constraints on subtour from an available infeasable solution
// conform to all types of cuts (e.g. static, lazy, callback, etc...)
#define CUT_FLAGS_MINCUT_ENABLED 0b00000001
#define CUT_LO_COEFF_DECAY 0.05
#define CUT_HI_COEFF_DECAY 0.3
inline double which_decay(ModelType mt) {
	return (model_variant(mt) & MODEL_VAR_HI_DEC) ?
		CUT_HI_COEFF_DECAY : CUT_LO_COEFF_DECAY;
}
int add_sec_on_subtours(void* env, void* lp, Instance* inst, void* args, double* xstar, int purgeable, int flags, int local);
int CC_add_sec_on_subtours(void* env, void* lp, Instance* inst, void* args, double* xstar, int purgeable, int flags, int local);

// callback function to add a user cut in CPlex from Concorde
int doit_fn_concorde2cplex(double cutval, int cutcount, int* cut, void* args);


// ********** MODELS **********
void build_model_base_undirected(OptData* optdata);

// ******** SOLUTIONS *********
void solve_benders(OptData* optdata);
void solve_callback(OptData* optdata);
void solve_symmetric_tsp(OptData* optdata);
// ***********************************************************************************

// ****************************** ASYMMETRIC TSP MODELS ******************************
// *********** CUTS ***********
// add SE constraints on pairs of nodes
void add_sec2_asymmetric(OptData* optdata);

// ********** MODELS **********
void build_model_base_directed(OptData* optdata);
void build_model_mtz(OptData* optdata);
void build_model_gg(OptData* optdata);

// ******** SOLUTIONS *********
void solve_asymmetric_tsp(OptData* optdata);
// ***********************************************************************************

#endif