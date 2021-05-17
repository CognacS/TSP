#ifndef GRAPH_PLOT_H_  

#define GRAPH_PLOT_H_

#include "tsp_data.h"
#include "tsp_utility.h"

#include <gnuplot_c.h>

h_GPC_Plot* setup_tsp_gnuplot();

void plot_tsp_xstar_undirected(graph* g, double* xstar);
void plot_tsp_xstar_directed(graph* g, double* xstar);

void plot_tsp_succ_undirected(graph* g, int* succ);
void plot_tsp_succ_directed(graph* g, int* succ);

void plot_tsp_solution_undirected(graph* g, Solution* sol);
void plot_tsp_solution_directed(graph* g, Solution* sol);

void plot_tsp_hardfixing_undirected(graph* g, int* succ, char* fixed, int* bridges);

void plot_heuristic_perflog(LinkedList* ll);

#endif