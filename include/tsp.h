#ifndef TSP_H_  

#define TSP_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>
#include <time.h>

#include <cplex.h>


/* ***********************************************************************************
*						GENERAL CONSTANTS
*********************************************************************************** */

//hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-9		// 1e-9		// very small numerical tolerance 
#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ
                                 
//data structures

/* ***********************************************************************************
*						INSTANCE GRAPH SECTION
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
} distancetype;

/**
* Graph data structure containing nnodes nodes, with coordinates xcoord and ycoord.
* Distances are calculated with the eclidian distance, which may be floating point or
* integer
*/
typedef struct {

	int nnodes;					// number of nodes of the graph
	
	double* xcoord;				// x coordinates of points
	double* ycoord;				// y coordinates of points

	double* tr_xcoord;			// transformed x coordinates of points (only used in some cases)
	double* tr_ycoord;			// transformed y coordinates of points (only used in some cases)

	int integer_costs;			// flag for integer costs (rounded distances)
	distancetype distance_type;	// type of distance to be computed on the nodes

} graph;

/* ***********************************************************************************
*						INSTANCE PARAMETERS SECTION
*********************************************************************************** */
#define DEF_MODEL_TYPE 0
#define DEF_RANDOMSEED 123456
#define DEF_TIMELIMIT CPX_INFBOUND
#define DEF_INPUT_FILE ""
#define DEF_BATCH_FILE ""
#define DEF_MAX_NODES -1
#define DEF_CUTOFF 0

/**
* Type xy of model to be solved, where:
* x = model specification
* y = variant of the specification
*/
typedef enum
{
	MTZ_S = 10,		// MTZ static constraints
	MTZ_L = 11,		// MTZ lazy constraints
	MTZ_STE = 12,	// MTZ subtour elimination constraints
	GG = 20			// GG
} modeltype;

inline int model_variant(modeltype mt) { return mt % 10; }

/**
* Parameters of the instance that defines the optimization procedure and the input data
*/
typedef struct {

	modeltype model_type;					// criteria for method selection based on the model type
	int randomseed;							// seed for RNG of CPlex
	double timelimit;						// overall time limit, in sec.s
	char input_file[100];		  			// input file
	char batch_file[100];		  			// batch file for instance runs automatization
	int max_nodes; 							// max n. of branching nodes in the final run (-1 unlimited)
	double cutoff; 							// cutoff (upper bound) for master

} params;

/* ***********************************************************************************
*						INSTANCE GLOBAL DATA SECTION
*********************************************************************************** */

/**
* Define the global data to be shared during optimization
*/
typedef struct {

	double z_best;							// best solution values

} global_data;

/* ***********************************************************************************
*						INSTANCE MODEL SECTION
*********************************************************************************** */
/**
* Define the model variables used during optimization
*/
typedef struct {

	int plc_holder;

} model;

/* ***********************************************************************************
*						INSTANCE DEFINITION AND FUNCTIONS
*********************************************************************************** */

/**
* Define the instance as a graph, its parameters, its global data and the model
*/
typedef struct {   
	
	graph inst_graph;
	params inst_params;
	global_data inst_global_data;
	model inst_model;

} instance;

#include "log.h"
#include "error.h"
#include "tsp_utility.h"
#include "cpx_models.h"

// to string functions for printing purpuses
void tostring_graph			(char* buffer, graph* inst_graph);
void tostring_params		(char* buffer, params* inst_p);
void tostring_global_data	(char* buffer, global_data* inst_global);
void tostring_model			(char* buffer, model* inst_model);
void tostring_instance		(char* buffer, instance* inst);

// fill an instance with default values
void fill_inst_default(instance* inst);

// free functions for destruction
void free_graph		(graph* inst_graph);
void free_instance	(instance* inst);

// optimization functions
int TSPopt(instance* inst, double** xstar);

//inline
inline int imax(int i1, int i2) { return ( i1 > i2 ) ? i1 : i2; } 
inline double dmin(double d1, double d2) { return ( d1 < d2 ) ? d1 : d2; } 
inline double dmax(double d1, double d2) { return ( d1 > d2 ) ? d1 : d2; }

#endif   /* TSP_H_ */ 
