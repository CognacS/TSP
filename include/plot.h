#ifndef PLOT_H_  

#define PLOT_H_

#include "solution.h"
#include "datastructs.h"
#include <gnuplot_c.h>

h_GPC_Plot* setup_tsp_gnuplot();

void plot_tsp_xstar_undirected(Graph* g, double* xstar);
void plot_tsp_xstar_directed(Graph* g, double* xstar);

void plot_tsp_succ_undirected(Graph* g, int* succ);
void plot_tsp_succ_directed(Graph* g, int* succ);

void plot_tsp_solution_undirected(Graph* g, Solution* sol);
void plot_tsp_solution_directed(Graph* g, Solution* sol);

void plot_tsp_hardfixing_undirected(Graph* g, int* succ, char* fixed, int* bridges);

void plot_heuristic_perflog(LinkedList* ll);

#endif