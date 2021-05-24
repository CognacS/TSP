#ifndef GRAPHS_H_  

#define GRAPHS_H_

#include "tsp_data.h"
#include "tsp_utility.h"

//void find_conncomps_kruskal(graph* g, const double* xstar, int* succ, int* comp, int* ncomp);

/**
* Find connected components on graph with edges xstar, and return the number
* of connected components
*/
int find_conncomps_dfs(graph* g, const double* xstar, int* succ, int* comp, int* ncomp);

char counterclockwise(graph* g, int p, int q, int r);
void build_convex_hull(graph* g, Solution* sol);
void furthest_nodes_sol(graph* g, Solution* sol);

#endif