#ifndef TSP_H_  

#define TSP_H_

#include <cplex.h>
#include "solution.h"
#include "random.h"
#include "tokens.h"
#include "chrono.h"
#include "datastructs.h"

/* ***********************************************************************************
*						GENERAL CONSTANTS
*********************************************************************************** */

#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ
                                 
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

} ModelType;


inline int model_tsptype(ModelType mt) { return (int)(mt / 100) % 10; }
inline int model_archetype(ModelType mt) { return (int)(mt / 10) % 10; }
inline int model_variant(ModelType mt) { return mt % 10; }

/**
* Parameters of the instance that defines the optimization procedure and the input data
*/
typedef struct {

	ModelType model_type;					// criteria for method selection based on the model type
	int randomseed;							// seed for RNG of CPlex
	double timelimit;						// overall time limit, in sec.s
	char input_file[100];		  			// input file
	char batch_file[100];		  			// batch file for instance runs automatization
	int max_nodes; 							// max n. of branching nodes in the final run (-1 unlimited)
	double cutoff; 							// cutoff (upper bound) for master
	char heuristic_code[200];				// code containing heuristic method information

} Params;

/* ***********************************************************************************
*						INSTANCE GLOBAL DATA SECTION
*********************************************************************************** */

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

} GlobalData;

/* ***********************************************************************************
*						INSTANCE MODEL SECTION
*********************************************************************************** */
#define DEF_NCOLS 0

/**
* Define the model variables used during optimization
*/
typedef struct {

	int ncols;				// number of columns in the MIP, more broadly the number of variables

} Model;

/* ***********************************************************************************
*						INSTANCE DEFINITION AND FUNCTIONS
*********************************************************************************** */

/**
* Define the instance as a graph, its parameters, its global data and the model
*/
typedef struct {   
	
	Graph inst_graph;
	Params inst_params;
	GlobalData inst_global_data;
	Model inst_model;

} Instance;

// ******************************* INSTANCE HANDLING *********************************
void print_params		(Params* inst_p);
void print_global_data	(GlobalData* inst_global);
void print_model		(Model* inst_model);
void print_instance		(Instance* inst);

// fill an instance with default values
void fill_inst_default(Instance* inst);

// free functions for destruction

void free_global_data	(GlobalData* inst_global);
void free_instance		(Instance* inst);

typedef enum
{
	TYPE_GRAPH = 0,
	TYPE_PARAM = 1,
	TYPE_GLOB  = 2,
	TYPE_MODEL = 3,
	TYPE_INST  = 4
} DataType;

// log datastructure function
void log_datastruct(void* object, DataType type, int runlvl, int loglvl);

// *************************** EXECUTION TIME HANDLING *******************************
double residual_time(Instance* inst);
int time_limit_expired(Instance* inst);

/* ***********************************************************************************
*						OPTIMIZATION DATA PACKAGE
*********************************************************************************** */
typedef struct
{
	CPXENVptr env;
	CPXLPptr lp;
} CplexData;

typedef struct
{
	Instance* inst;
	CplexData* cpx;
	LinkedList* perflog;
} OptData;

inline int using_cplex(OptData* optdata) { return optdata->cpx != NULL; }

#endif   /* TSP_ */