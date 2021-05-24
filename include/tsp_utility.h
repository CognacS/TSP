#ifndef TSP_UTILITY_H_  

#define TSP_UTILITY_H_

#define _USE_MATH_DEFINES
#include <math.h>

#include "tsp_data.h"
#include "log.h"
#include "error.h"

#define TYPE_GRAPH 0
#define TYPE_PARAM 1
#define TYPE_GLOB  2
#define TYPE_MODEL 3
#define TYPE_INST  4

// constants
#define RRR 6378.388	// used in geo distance

// ********* coordinates transformation *********
// master function
void coord_transform(graph* g);
// partial functions
void geo_transform(graph* g);
// utility functions
double degrees_to_radians(double coord);
// get coordinates
double ggetx(graph* g, int idx);
double ggety(graph* g, int idx);

// ********* distance computation functions *********
// master function
double dist(int i, int j, graph* g);
// partial functions
double euc_dist(int i, int j, graph* g);
double att_dist(int i, int j, graph* g);
double geo_dist(int i, int j, graph* g);
int dist_to_int(double distance);
// array functions
void compute_dists(double* dists, graph* g);
double delta_cost(int a, int b, int c, graph* g);

// variable index in CPLEX functions
int xpos(int i, int j, int nnodes);
int xxpos(int i, int j, int nnodes);
int upos(int i, int nnodes);
int ypos(int i, int j, int nnodes);

// variable inference of 0/1 value functions
int is_one(double num);
int is_zero_strict(double num);
int is_one_strict(double num);

// FUNCTIONS FOR HANDLING SOLUTIONS
// count how many edges are activated in xstar
size_t count_active_edges(size_t size, double* xstar);
// extract active edges and return a list of pairs (i, j) with indices (2k, 2k+1), k = 0,1,...,nnodes
void extract_active_edges(int nnodes, double* xstar, int* ecount, int** elist);
// produce integer distances from a list of edges
void compute_idx(int nnodes, int ecount, int* elist, int** idxlist);
// produce integer distances from a list of edges
void compute_idx_dst(graph* g, int ecount, int* elist, int** idxlist, unsigned int** dstlist);
// radix sort: sort edges by distance. Warning: arguments are changed!!!
void radix_sort(int ecount, int** idxlist, unsigned int** dstlist);

// conversion methods for a feasable solution
int xstar2succ(double* xstar, int* succ, int nnodes);	// O(n^2) complexity
int succ2xstar(int* succ, double* xstar, int nnodes);	// O(n) complexity
int succ2chromo(int* succ, int* chromo, int nnodes);	// O(n) complexity
int chromo2succ(int* chromo, int* succ, int nnodes);	// O(n) complexity
int convert_solution(Solution* sol, solformat format);

// handle solution
void empty_solution(Solution* sol, solformat format, int nnodes);
void clear_solution(Solution* sol);

/**
* add the edge (i,j) into the solution with the following complexity:
* -	xstar	O(1)
* -	succ	O(1)
* -	chromo	O(n)	inefficient!
*/
void add_edge_solution(Solution* sol, int i, int j, double cost);
char rem_edge_solution(Solution* sol, int i, int j, double cost);
void shallow_copy_solution(Solution* dst, Solution* src);

typedef struct
{
	Solution* sol;
	int start;
	int curr;
	int info;
} SolutionIterator;

void initialize_sol_iterator(SolutionIterator* iter, Solution* sol, int start_node);
char next_node_in_solution(SolutionIterator* iter);

int compute_list_unvisited_nodes(Solution* sol, int* nodelist);

// log datastructure function
void log_datastruct(void* object, int type, int runlvl, int loglvl);

// print solutions
void print_directed_sol(graph* g, double* xstar);
void print_undirected_sol(graph* g, double* xstar);
void print_xstar(size_t size, double* xstar);

#define FORMAT2FORMAT(sol, f1, f2, result)\
if ((f2) == SOLFORMAT_XSTAR)\
{\
	arr_malloc_s(sol->xstar, (int)(sol->nnodes*(sol->nnodes-1)/2.0), double);\
	if		((f1) == SOLFORMAT_SUCC)	result = succ2xstar(sol->succ, sol->xstar, sol->nnodes);\
	else if ((f1) == SOLFORMAT_CHROMO)	print_error(ERR_GENERIC_NOT_IMPL, "chromo2xstar");\
}\
else if ((f2) == SOLFORMAT_SUCC)\
{\
	arr_malloc_s(sol->succ, sol->nnodes, int);\
	if		((f1) == SOLFORMAT_XSTAR)	result = xstar2succ(sol->xstar, sol->succ, sol->nnodes);\
	else if ((f1) == SOLFORMAT_CHROMO)	result = chromo2succ(sol->chromo, sol->succ, sol->nnodes);\
}\
else if ((f2) == SOLFORMAT_CHROMO)\
{\
	arr_malloc_s(sol->chromo, sol->nnodes, int);\
	if		((f1) == SOLFORMAT_XSTAR)	print_error(ERR_GENERIC_NOT_IMPL, "xstar2chromo");\
	else if ((f1) == SOLFORMAT_SUCC)	result = succ2chromo(sol->succ, sol->chromo, sol->nnodes);\
}

#define REMOVE_FORMAT(sol, f)\
if		(f == SOLFORMAT_XSTAR)	free_s(sol->xstar);\
else if (f == SOLFORMAT_SUCC)	free_s(sol->succ);\
else if (f == SOLFORMAT_CHROMO) free_s(sol->chromo);

#endif
