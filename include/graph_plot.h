#ifndef GRAPH_PLOT_H_  

#define GRAPH_PLOT_H_

#include "tsp.h"
#include "tsp_utility.h"

#include <gnuplot_c.h>

void plot_tsp_solution_undirected(graph* g, double* xstar);
void plot_tsp_solution_directed(graph* g, double* xstar);

#endif