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

// log datastructure function
void log_datastruct(void* object, int type, int runlvl, int loglvl);

// print solutions
void print_directed_sol(graph* g, double* xstar);
void print_undirected_sol(graph* g, double* xstar);


#endif
