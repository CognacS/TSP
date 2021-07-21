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

/**
* Decompose heuristic strings using the following structure:
* B(...,...)C(...,...)R(...,...)
* ##### B = BACKBONE HEURISTICS #####
* - iterative local search:			i(construct_timelim,refine_timelim)
* - variable neighborhood search:	v(//,//,maxkick_k)
* - tabu search:					t(//,//,min_tenure,max_tenure,period,intermediate_levels)
* - genetic algorithm:				g(//,//,crossover_variant,pool_size,elite_ratio,mutation_prob)
*	> crossover variants = {aex: a, pmx: p, mix: m}
* ##### C = CONSTRUCTIVE HEURISTICS #####
* - cplex tsp solver:				s(,)
* - greedy construction:			g(0,any)
* - greedy construction (GRASP):	g(p_deviation,candidate_pool_size)
* - extra mileage ins.:				e(0,any,variant)
* - extra mileage ins. (GRASP):		e(p_deviation,candidate_pool_size,variant)
*	> extramil variants = {hull: h, random edge: r, furthest nodes: f}
* ##### R = REFINING HEURISTICS #####
* - hard fixing:					h(variant,max_rounds_noimprovement,values[5])
*	> variants = {scheduled[v0,v1,v2,v3,v4]: s, fixed[v0,...]: f, uniform[v0,v1,...]: u}
* - local branching:				l(max_rounds_noimprovement,values[4])
*	> only scheduled
* - mix of hard fix and local br.:	m(max_rounds_noimprovement,values[6])
*	> first 3 values are for scheduled hard fixing, last 3 for local branching
* - 2-OPT iterative refinement:		2(,)						
*/
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

typedef struct
{
	int min_tenure;
	int max_tenure;
	int period;
	int interm_levels;
} TabuSearchData;

typedef struct
{
	char crossover_variant;
	int pool_size;
	double elite_ratio;
	double mutation_p;

} GeneticAlgData;

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
	int num_values;
	int rounds_nonimprove;
	int max_rounds;
	int progress;
	double values[HEUR_HARDFIX_VALUES];
} HardfixingData;

#define HEUR_LOCALBRANCH_VALUES 4
typedef struct
{
	char variant;
	int num_values;
	int rounds_nonimprove;
	int max_rounds;
	int progress;
	double values[HEUR_LOCALBRANCH_VALUES];
} LocalbranchingData;

#define HEUR_MIXHARDSOFT_VALUES 6
typedef struct
{
	char variant;
	int num_values;
	int rounds_nonimprove;
	int max_rounds;
	int progress;
	double values[HEUR_MIXHARDSOFT_VALUES];
} MixHardSoftData;

typedef struct
{
	TabuList* tabu;
	char allow_worsening;
} Opt2MoveData;

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
void backbone_tabu_search(OptData* optdata, Solution* sol, Heuristic* heur, void* data);
void backbone_genetic_algorithm(OptData* optdata, Solution* sol, Heuristic* heur, void* data);

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
void refine_mixhardsoft(OptData* optdata, Solution* sol, void* data, double timelim);
void refine_2opt(OptData* optdata, Solution* sol, void* data, double timelim);


#endif


