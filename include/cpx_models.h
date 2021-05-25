#ifndef CPX_MODELS_H_  

#define CPX_MODELS_H_

#include "random.h"
#include <cut.h>

#include "optimizer_data.h"
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
	double (*prob_function)(int, double);
	double prob_decay;
} callback_instance;

typedef struct
{
	instance* inst;
	CPXCALLBACKCONTEXTptr context;
	double* value;
	int* index;
	int* labels;
	int local;
} concorde_instance;


// ****************************** SYMMETRIC TSP MODELS *******************************
// ******** CALLBACKS *********
int CPXPUBLIC sec_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle);

// *********** CUTS ***********
// add SE constraints on subtour from an available infeasable solution
// conform to all types of cuts (e.g. static, lazy, callback, etc...)
#define CUT_FLAGS_MINCUT_ENABLED 0b00000001
#define CUT_LO_COEFF_DECAY 0.05
#define CUT_HI_COEFF_DECAY 0.3
inline double which_decay(modeltype mt) {
	return (model_variant(mt) & MODEL_VAR_HI_DEC) ?
			CUT_HI_COEFF_DECAY : CUT_LO_COEFF_DECAY;
}
int add_sec_on_subtours(void* env, void* lp, instance* inst, void* args, double* xstar, int purgeable, int flags, int local);
int CC_add_sec_on_subtours(void* env, void* lp, instance* inst, void* args, double* xstar, int purgeable, int flags, int local);

// callback function to add a user cut in CPlex from Concorde
int doit_fn_concorde2cplex(double cutval, int cutcount, int* cut, void* args);

// 2-opt move
double move_2opt(int* succ, graph* g, char allow_unimproving);
double remove_crossings(int* succ, graph* g);

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