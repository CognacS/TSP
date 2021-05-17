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

// ********* distance computation functions *********
// master function
double dist(int i, int j, graph* g);
// partial functions
double euc_dist(int i, int j, graph* g);
double att_dist(int i, int j, graph* g);
double geo_dist(int i, int j, graph* g);
int dist_to_int(double distance);

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
int convert_solution(Solution* sol, solformat format);

// handle solution
void empty_solution(Solution* sol, solformat format, int nnodes);

// log datastructure function
void log_datastruct(void* object, int type, int runlvl, int loglvl);

// print solutions
void print_directed_sol(graph* g, double* xstar);
void print_undirected_sol(graph* g, double* xstar);
void print_xstar(size_t size, double* xstar);


#endif
