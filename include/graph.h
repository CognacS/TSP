#ifndef GRAPH_H_  

#define GRAPH_H_

#define _USE_MATH_DEFINES
#include <math.h>

#include "global.h"
#include "error.h"

/* ***********************************************************************************
*						GRAPH STRUCTURE DEFINITION
*********************************************************************************** */
#define DEF_INTEGER_COSTS 0

/**
* Type of distance to be computed on the graph
*/
typedef enum
{
	ATT,
	EUC_2D,
	GEO
} DistanceType;

/**
* Graph data structure containing nnodes nodes, with coordinates xcoord and ycoord.
* Distances are calculated with the eclidian distance, which may be floating point or
* integer
*/
typedef struct {

	int nnodes;					// number of nodes of the graph

	double* xcoord;				// x coordinates of points
	double* ycoord;				// y coordinates of points

	double* tr_xcoord;			// transformed x coordinates of points (useful for some distypes)
	double* tr_ycoord;			// transformed y coordinates of points (useful for some distypes)

	// distance matrix containing distance for each edge (i,j).
	// Note that dist(i,j) = dist(j,i), also dist(i,i) is undefined
	double** distance_matrix;

	char integer_costs;			// flag for integer costs (rounded distances)
	DistanceType distance_type;	// type of distance to be computed on the nodes

} Graph;


/* ***********************************************************************************
*									GRAPH FUNCTIONS
*********************************************************************************** */

// ********************************* GRAPH HANDLING **********************************
void print_graph(Graph* g);
void default_graph(Graph* g);
void free_graph(Graph* g);

// *********************************** COORDINATES ***********************************
// constants
#define RRR 6378.388	// used in geo distance

// ********* coordinates transformation *********
// master function
void coord_transform(Graph* g);
// partial functions
void geo_transform(Graph* g);
// utility functions
double degrees_to_radians(double coord);
// get coordinates (when distance type may be ambiguous)
double ggetx(Graph* g, int node);
double ggety(Graph* g, int node);

// ************************************* DISTANCE ************************************

// master function
double dist(Graph* g, int i, int j);
double calc_dist(Graph* g, int i, int j);
// partial functions
double euc_dist(Graph* g, int i, int j);
double att_dist(Graph* g, int i, int j);
double geo_dist(Graph* g, int i, int j);
int dist_to_int(double distance);
// setup funtions
void compute_distance_matrix(Graph* g);
// distance-related functions
double delta_cost(Graph* g, int a, int b, int c);

// ************************************ INDEXING *************************************
// variable index in CPLEX functions
int xpos(int i, int j, int nnodes);
int xxpos(int i, int j, int nnodes);

// ********************************** GENERAL PURPOSE ********************************
// count all edges (not self edges) in a group of nodes
inline int complete_graph_edges(int nnodes) { return (int)(nnodes * (nnodes - 1) / 2.0); }
char counterclockwise(Graph* g, int p, int q, int r);

// ************************************* ALGORITHMS **********************************
/**
* Find connected components on graph with edges xstar, and return the number
* of connected components
*/
int find_conncomps_dfs(Graph* g, const double* xstar, int* succ, int* comp, int* ncomp);



#endif