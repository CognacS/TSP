#ifndef HEURISTICS_H_  

#define HEURISTICS_H_

#include "cpx_models.h"
#include "chrono.h"

// ****************************** HEURISTIC DATASTRUCTURE ****************************
// prefer succ format
#define HEUR_PREFERRED_FORMAT SOLFORMAT_SUCC

typedef struct Heuristic
{
	void (*backbone_heur) (OptData*, Solution*, struct Heuristic*, void*);
	void* backbone_data;
	solformat backbone_sol_format;

	void (*construct_heur) (OptData*, Solution*, void*, double);
	void* construct_data;
	solformat construct_sol_format;

	void (*refine_heur) (OptData*, Solution*, void*, double);
	void* refine_data;
	solformat refine_sol_format;

	char requires_cplex;

} Heuristic;

// ****************************** BACKBONES DATASTRUCTURES ***************************
typedef struct
{
	double construct_timelimit;
	double refine_timelimit;
} IterativeData;

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
	int progress;
} HardfixingData;

#define HEUR_LOCALBRANCH_VALUES 4
typedef struct
{
	char variant;
	int values[HEUR_LOCALBRANCH_VALUES];
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

// ********************************* BACKBONE METHODS ********************************
void solve_heuristically(OptData* optdata);
/**
* solve problem by constructing a starting solution, and refine it iteratively until
* the time limit is hit
*/
void backbone_iter_local_search(OptData* optdata, Solution* sol, Heuristic* heur, void* data);

// ***************************** CONSTRUCTIVE HEURISTICS *****************************
// solution construction procedures which generate a starting point from which to
// improve using refinement techniques
void construct_tspsolver(OptData* optdata, Solution* sol, void* data, double timelim);
void construct_greedy(OptData* optdata, Solution* sol, void* data, double timelim);
void construct_extramileage(OptData* optdata, Solution* sol, void* data, double timelim);

// ***************************** REFINEMENT HEURISTICS *****************************
// solution refinement procedures which expect a feasable solution and always return a
// better solution if it is possible
#define IMPROVEMENT_RATIO 0.90

void refine_hardfixing(OptData* optdata, Solution* sol, void* data, double timelim);
void refine_localbranching(OptData* optdata, Solution* sol, void* data, double timelim);

#endif


