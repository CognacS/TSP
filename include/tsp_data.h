#ifndef TSP_DATA_H_  

#define TSP_DATA_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#define DEF_CUTOFF 1e-3

// tree of implemented models
#define MODEL_TSP_ASYMM 1
#define	MODEL_AS_MTZ 1
#define	MODEL_AS_GG	2

#define MODEL_TSP_SYMM 2
#define	MODEL_SY_BEND 1
#define MODEL_SY_CLBCK 2

// masks for variants
#define MODEL_VAR_STAT		0b00000000
#define MODEL_VAR_LZ		0b00000001
#define MODEL_VAR_NO_SEC	0b00000000
#define MODEL_VAR_SEC		0b00000010

/**
* Type xyz of model to be solved, where:
* x = TSP problem variant
* y = model specification (archetype) for the proposed TSP problem
* z = variant of the specification
*/
typedef enum
{
	MTZ_ST =		110,	// MTZ static constraints
	MTZ_LZ =		111,	// MTZ lazy constraints
	MTZ_ST_SEC =	112,	// MTZ static constraints + SEC
	MTZ_LZ_SEC =	113,	// MTZ lazy constraints + SEC

	GG_ST =			120,	// GG  static constraints
	GG_LZ =			121,	// GG  lazy constraints
	GG_ST_SEC =		122,	// GG  static constraints + SEC
	GG_LZ_SEC =		123,	// GG  lazy constraints + SEC

	BENDERS =		210,	// Benders' method
	CALLBACK =		220		// Callback's method
} modeltype;


inline int model_tsptype(modeltype mt) { return (int)(mt / 100) % 10; }
inline int model_archetype(modeltype mt) { return (int)(mt / 10) % 10; }
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
#define DEF_XSTAR NULL

/**
* Define the global data to be shared during optimization
*/
typedef struct {

	double tstart;		// time from start of the optimization procedure
	double texec;		// time of execution updated during optimization

	double* xstar;		// best solution available
	double zbest;		// obj value of the best integer solution
	double lbbest;		// lower bound of the best solution

} global_data;

/* ***********************************************************************************
*						INSTANCE MODEL SECTION
*********************************************************************************** */
#define DEF_NCOLS 0

/**
* Define the model variables used during optimization
*/
typedef struct {

	int ncols;

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

// functions for printing purposes
void print_graph		(graph* inst_graph);
void print_params		(params* inst_p);
void print_global_data	(global_data* inst_global);
void print_model		(model* inst_model);
void print_instance		(instance* inst);

// fill an instance with default values
void fill_inst_default(instance* inst);

// free functions for destruction
void free_graph			(graph* inst_graph);
void free_global_data	(global_data* inst_global);
void free_instance		(instance* inst);

#endif   /* TSP_DATA_H_ */ 
