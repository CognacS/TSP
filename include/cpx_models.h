#ifndef CPX_MODELS_H_  

#define CPX_MODELS_H_

#include "mip_utility.h"
#include "graphs.h"
#include "graph_plot.h"


// ****************************** SYMMETRIC TSP MODELS *******************************
// *********** CUTS ***********
// add SE constraints on subtour from an available infeasable solution
// conform to all types of cuts (e.g. static, lazy, callback, etc...)
int add_sec_on_subtours(void* env, void* cbdata, instance* inst, double* xstar, int wherefrom, int purgeable);

// ********** MODELS **********
void build_model_base_undirected(instance* inst, CPXENVptr env, CPXLPptr lp);

// ******** SOLUTIONS *********
void solve_benders(instance* inst, CPXENVptr env, CPXLPptr lp);
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