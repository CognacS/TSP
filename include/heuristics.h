#ifndef HEURISTICS_H_  

#define HEURISTICS_H_

#include "cpx_models.h"
#include "chrono.h"

// *********************************** HEURISTIC CODE ********************************
#define HEUR_SECTIONS 3
#define HEUR_NUMPARAM_MAX 10

typedef enum
{
	HEUR_BACKBONE	= 0,
	HEUR_CONSTRUCT	= 1,
	HEUR_REFINE		= 2
} HeuristicSection;

typedef struct
{
	char methods[HEUR_SECTIONS];					// heuristic method chars for each section
	char* params[HEUR_SECTIONS][HEUR_NUMPARAM_MAX];	// method-specific parameters for each section

} HeuristicCode;

void decompose_heuristic_code(HeuristicCode* decomp, char* heuristic_code);

// ****************************** HEURISTIC DATASTRUCTURE ****************************
// prefer succ format
#define HEUR_PREFERRED_FORMAT SOLFORMAT_SUCC

typedef struct Heuristic
{
	void (*backbone_heur) (OptData*, Solution*, struct Heuristic*, void*);
	void* backbone_data;
	SolFormat backbone_sol_format;

	void (*construct_heur) (OptData*, Solution*, void*, double);
	void* construct_data;
	SolFormat construct_sol_format;
	double construct_timelimit;

	void (*refine_heur) (OptData*, Solution*, void*, double);
	void* refine_data;
	SolFormat refine_sol_format;
	double refine_timelimit;

	char requires_cplex;
	char to_timelimit;

} Heuristic;

// ****************************** BACKBONES DATASTRUCTURES ***************************
typedef struct
{
	int max_k;
} VNSData;

// ************************** CONSTRUCTIVE HEUR DATASTRUCTURES ***********************
typedef struct
{
	double p_dev;		// probability of deviating from the greedy choice
	int candpool_size;	// size of the candidate pool
} GraspData;

typedef struct
{
	GraspData grasp;
} GreedyAlgoData;

typedef struct
{
	GraspData grasp;
	char variant_startingchoice;	// variant on the starting edges choice
} ExtraMileageData;

// *************************** REFINEMENT HEUR DATASTRUCTURES ************************
#define HEUR_HARDFIX_VALUES 5
typedef struct
{
	char variant;
	double values[HEUR_HARDFIX_VALUES];
	int rounds_nonimprove;
	int max_rounds;
	int progress;
} HardfixingData;

#define HEUR_LOCALBRANCH_VALUES 4
typedef struct
{
	char variant;
	int values[HEUR_LOCALBRANCH_VALUES];
	int rounds_nonimprove;
	int max_rounds;
	int progress;
} LocalbranchingData;

// ************************** SOLUTION HANDLING DATASTRUCTURES ***********************
typedef struct
{
	// best solution found until now
	Solution* best_sol;

	// current solution in use
	Solution* curr_sol;

} SolutionPool;

// ******************************** AUXILIARY FUNCTIONS ******************************
void decode_heuristic(OptData* optdata, Heuristic* heur);
void resolve_solformats(Heuristic* heur);
char extramileage_move(Solution* sol, Graph* g, SetOfNodes* ext_nodes, GraspData* grasp);
void kopt_kick(Solution* sol, Graph* g, int k);

// ********************************* BACKBONE METHODS ********************************
void solve_heuristically(OptData* optdata);
/**
* solve problem by constructing a starting solution, and refine it iteratively until
* the time limit is hit
*/
void backbone_iter_local_search(OptData* optdata, Solution* sol, Heuristic* heur, void* data);
void backbone_var_neighborhood_search(OptData* optdata, Solution* sol, Heuristic* heur, void* data);


// ***************************** CONSTRUCTIVE HEURISTICS *****************************
// solution construction procedures which generate a starting point from which to
// improve using refinement techniques
void construct_tspsolver(OptData* optdata, Solution* sol, void* data, double timelim);
void construct_greedy(OptData* optdata, Solution* sol, void* data, double timelim);
void construct_extramileage(OptData* optdata, Solution* sol, void* data, double timelim);

// ***************************** REFINEMENT HEURISTICS *****************************
// solution refinement procedures which expect a feasable solution and always return a
// better solution if it is possible
#define IMPROVEMENT_RATIO 0.995

void refine_hardfixing(OptData* optdata, Solution* sol, void* data, double timelim);
void refine_localbranching(OptData* optdata, Solution* sol, void* data, double timelim);
void refine_2opt(OptData* optdata, Solution* sol, void* data, double timelim);


#endif


