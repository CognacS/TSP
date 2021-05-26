#ifndef TSP_DATA_H_  

#define TSP_DATA_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>
#include "random.h"
#include "tokens.h"


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

#define MODEL_HEURISTICS 3
// heuristics methods are implemented through the heuristic code
// due to the complexity in possible combinations

// masks for variants in asymmetric TSP
#define MODEL_VAR_STAT		0b00000000	// static constraints
#define MODEL_VAR_LZ		0b00000001	// lazy constraints
#define MODEL_VAR_NO_SEC	0b00000000	// SEC on pairs disabled
#define MODEL_VAR_SEC		0b00000010	// SEC on pairs enabled

// masks for variants in symmetric TSP
#define MODEL_VAR_LO_DEC	0b00000000	// use low decay coefficient
#define MODEL_VAR_HI_DEC	0b00000001	// use high decay coefficient
#define MODEL_VAR_SEP_HYP	0b00000000	// use hyperbolic decay
#define MODEL_VAR_SEP_EXP	0b00000010	// Concorde separation with exp decay
#define MODEL_VAR_SEP_FIX	0b00000100	// use low decay coefficient
#define MODEL_VAR_SEP_CUT	0b00000110	// use high decay coefficient
#define MODEL_VAR_NOSEP		0b00001000	// use Concorde procedure for rejects
#define MODEL_VAR_SEP_MSK	0b00000110	// get separation prob method bits

/**
* Type xyz of model to be solved, where:
* x = TSP problem variant
* y = model specification (archetype) for the proposed TSP problem
* z = variant of the specification
*/
typedef enum
{
	// COMPACT MODELS FOR SOLVING EXACT ASYMMETRIC TSP
	MTZ_ST =			110,	// MTZ static constraints
	MTZ_LZ =			111,	// MTZ lazy constraints
	MTZ_ST_SEC =		112,	// MTZ static constraints + SEC
	MTZ_LZ_SEC =		113,	// MTZ lazy constraints + SEC

	GG_ST =				120,	// GG  static constraints
	GG_LZ =				121,	// GG  lazy constraints
	GG_ST_SEC =			122,	// GG  static constraints + SEC
	GG_LZ_SEC =			123,	// GG  lazy constraints + SEC

	// METHODS FOR SOLVING EXACT SYMMETRIC TSP
	BENDERS =			210,	// Benders' method

	CLBK_HYP_LO =		220,	// Callback's method with Separation: Concorde + hyperbolic low decay
	CLBK_HYP_HI =		221,	// Callback's method with Separation: Concorde + hyperbolic high decay
	CLBK_EXP_LO =		222,	// Callback's method with Separation: Concorde + exponential low decay
	CLBK_EXP_HI =		223,	// Callback's method with Separation: Concorde + exponential high decay
	CLBK_FIX_LO =		224,	// Callback's method with Separation: Concorde + fixed low prob
	CLBK_FIX_HI =		225,	// Callback's method with Separation: Concorde + fixed high prob
	CLBK_CUT_LO =		226,	// Callback's method with Separation: Concorde + cutoff low depth
	CLBK_CUT_HI =		227,	// Callback's method with Separation: Concorde + cutoff high depth
	CLBK_NOSEP	=		228,	// Callback's method without Separation

	// HEURISTICS
	HEURISTICS	=		300

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
	char heuristic_code[200];				// code containing heuristic method information

} params;

/* ***********************************************************************************
*						INSTANCE GLOBAL DATA SECTION
*********************************************************************************** */
#define DEF_XSTAR NULL
#define DEF_SUCC NULL

#define NUM_SOLFORMATS 3
typedef enum
{
	SOLFORMAT_SUCC		= 0b00000001,
	SOLFORMAT_XSTAR		= 0b00000010,
	SOLFORMAT_CHROMO	= 0b00000100
} solformat;

#define SOLFORMAT_NOFORMAT 0b00000000
#define SOLFORMAT_ALLFORMAT (SOLFORMAT_SUCC | SOLFORMAT_XSTAR | SOLFORMAT_CHROMO)

/**
* Solution package for compatibility purpose
*/
typedef struct
{
	double* xstar;		// xstar format compatible with cplex methods
	int* succ;			// succ format compatible with heuristic methods
	int* chromo;		// chromosome format compatible with genetic algorithms

	int handle_node;	// node which is assured to belong to the tour
	int nnodes;			// size of the solution
	double cost;		// cost of solution
	solformat format;	// describes which format(s) is(are) currently correct
	char optimal;		// flag to signal global (or local) optimum

} Solution;

/**
* Global data to be shared during optimization
*/
typedef struct {

	double tstart;		// time from start of the optimization procedure
	double texec;		// time of execution updated during optimization

	double* xstar;		// best solution available in xstar format
	double zbest;		// obj value of the best integer solution
	double lbbest;		// lower bound of the best solution

	double perf_measure; // measure to compare performances in different methods

} global_data;

/* ***********************************************************************************
*						INSTANCE MODEL SECTION
*********************************************************************************** */
#define DEF_NCOLS 0

/**
* Define the model variables used during optimization
*/
typedef struct {

	int ncols;				// number of columns in the MIP, more broadly the number of variables

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

#include "error.h"

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
void free_solution		(Solution* sol);
void free_instance		(instance* inst);

/* ***********************************************************************************
*						AUXILIARY STRUCTURES AND FUNCTIONS
*********************************************************************************** */

// ************* solformats *************
int solformat2num(solformat format);
solformat solformat_collapse(solformat format);

// ************* linked cell *************
typedef struct Cell
{
	double x;
	double y;
	struct Cell* next;
} Cell;

// ************* linked list *************
typedef struct
{
	Cell* start;
	Cell* end;
} LinkedList;

LinkedList* newLinkedList();
int LL_is_empty(LinkedList* list);
void LL_add_value(LinkedList* list, double x, double y);
void LL_free(LinkedList* list);

// ************* ordered indexed array *************
#define MAX_INDEX_SIZE 3
typedef struct
{
	int arr[MAX_INDEX_SIZE];
} Index;

typedef struct
{
	Index index;
	double value;
} IndexedValue;

// Ordered Indexed array has elements ordered from bigger to smaller, i.e.: oiarr[0] has the biggest value
inline char OIA_eligible(IndexedValue* oiarr, double newvalue) { return oiarr[0].value > newvalue; }
char OIA_insert(IndexedValue* oiarr, int size, IndexedValue newelem);
void OIA_clear(IndexedValue* oiarr, int size);
IndexedValue OIA_pack1(int index, double value);
IndexedValue OIA_pack3(int a, int b, int c, double value);
inline IndexedValue OIA_best(IndexedValue* oiarr, int size) { return oiarr[size-1]; }
IndexedValue OIA_choose(IndexedValue* oiarr, int size, char includebest);

void IV_deepcopy(IndexedValue* dst, IndexedValue* src);


#endif   /* TSP_DATA_H_ */ 
