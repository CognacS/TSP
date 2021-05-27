#ifndef CPX_UTILITY_H_  

#define CPX_UTILITY_H_

#include <cut.h>
#include "tsp.h"

typedef enum
{
	CUT_STATIC = -1,
	CUT_LAZY = -2,
	CUT_USER = -3,
	CUT_CALLBACK_LOCAL = -4,
	CUT_CALLBACK_REJECT = -5
} CutType;


int upos(int i, int nnodes);
int ypos(int i, int j, int nnodes);

// ************************** SEPARATION PROB FUNCTIONS ******************************
// functions for computing the separation probability
inline double hyp_decay_prob(int depth, double decay) { return 1.0 / (decay * depth + 1); }
inline double exp_decay_prob(int depth, double decay) { return exp(-decay * depth); }
inline double fixed_prob(int depth, double decay) { return decay + 0.2; }
inline double cutoff_depth(int depth, double decay) { return depth < (decay * 100 + 1); }

// ***************************** MODEL TYPE FUNCTIONS ********************************
/**
* takes the model type and returns the type of cut/constraint to give to cpx_add_cut
* in other words it matches the variant of mt with the bit mask for lazy constraints
*/
inline int variant2constr(ModelType mt) { return CUT_STATIC - (model_variant(mt) & MODEL_VAR_LZ); }
/**
* takes the model type and returns whether SEC's are needed
* in other words it matches the variant of mt with the bit mask for SEC
*/
inline int need_sec(ModelType mt) { return model_variant(mt) & MODEL_VAR_SEC; }

/**
* takes the model type and returns whether to use Concorde for computing separation
*/
inline int use_cc_on_sep(ModelType mt) { return !(model_variant(mt) & MODEL_VAR_NOSEP); }


// ************************* INTERFACE FUNCTIONS TO CPLEX ****************************
/**
* Unified driver to add constraints/cuts to the LP.
* */
void cpx_add_cut(void* env, void* lp,
	int nnz, double rhs, char sense,
	int* index, double* value, int purgeable, char* name, int local);

// setup cplex parameters from instance
void cpx_setup_cplex(OptData* optdata);
// closes cplex problem and environment
void cpx_close_cplex(OptData* optdata);

int cpx_solution_available(OptData* optdata);

void cpx_warmstart(OptData* optdata, double* xstar);

void cpx_extract_sol_obj_lb(OptData* optdata, char* zone_name);
void cpx_extract_sol_obj(OptData* optdata, Solution* sol, char* zone_name);

void cpx_timelimit(OptData* optdata, double timelimit);


#endif