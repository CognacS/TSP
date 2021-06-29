#ifndef SOLUTION_H_  

#define SOLUTION_H_

#include "graph.h"
#include "datastructs.h"

/* ***********************************************************************************
*						SOLUTION STRUCTURE DEFINITION
*********************************************************************************** */
#define NUM_SOLFORMATS 3
typedef enum
{
	SOLFORMAT_SUCC = 0b00000001,
	SOLFORMAT_XSTAR = 0b00000010,
	SOLFORMAT_CHROMO = 0b00000100
} SolFormat;

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
	SolFormat format;	// describes which format(s) is(are) currently correct
	char optimal;		// flag to signal global (or local) optimum

} Solution;


/* ***********************************************************************************
*									SOLUTION FUNCTIONS
*********************************************************************************** */

// ******************************* SOLUTION HANDLING *********************************
void print_solution(Solution* sol);
void create_solution(Solution* sol, SolFormat format, int nnodes);
void erase_solution(Solution* sol);
void shallow_copy_solution(Solution* dst, Solution* src);
void deep_copy_solution(Solution* dst, Solution* src);
void free_solution(Solution* sol);

// print solutions
void print_directed_sol(Graph* g, double* xstar);
void print_undirected_sol(Graph* g, double* xstar);
void print_xstar(int size, double* xstar);

#define DEEP_COPY_ARR(dst, src, size)\
do{\
for(int i = 0; i < (size); i++)\
	dst[i] = src[i];\
} while(0)

// **************************** SOLUTION FORMAT FUNCTIONS ****************************
int solformat2num(SolFormat format);
SolFormat solformat_collapse(SolFormat format);

// ********************************** COMPUTE COSTS **********************************
double compute_cost(int* chromo, Graph* g);

// *********************************** CONVERSION ************************************
// conversion methods for a feasable solution

int xstar2succ(double* xstar, int* succ, int nnodes);	// O(n^2) complexity
int succ2xstar(int* succ, double* xstar, int nnodes);	// O(n) complexity
int succ2chromo(int* succ, int* chromo, int nnodes);	// O(n) complexity
int chromo2succ(int* chromo, int* succ, int nnodes);	// O(n) complexity
int chromo2index(int* chromo, int* index, int nnodes);	// O(n) complexity

int convert_solution(Solution* sol, SolFormat format);

#define FORMAT2FORMAT(sol, f1, f2, result)\
if ((f2) == SOLFORMAT_XSTAR)\
{\
	arr_malloc_s(sol->xstar, (int)(sol->nnodes*(sol->nnodes-1)/2.0), double);\
	if		((f1) == SOLFORMAT_SUCC)	result = succ2xstar(sol->succ, sol->xstar, sol->nnodes);\
	else if ((f1) == SOLFORMAT_CHROMO)	print_error(ERR_GENERIC_NOT_IMPL, "chromo2xstar");\
}\
else if ((f2) == SOLFORMAT_SUCC)\
{\
	arr_malloc_s(sol->succ, sol->nnodes, int);\
	if		((f1) == SOLFORMAT_XSTAR)	result = xstar2succ(sol->xstar, sol->succ, sol->nnodes);\
	else if ((f1) == SOLFORMAT_CHROMO)	result = chromo2succ(sol->chromo, sol->succ, sol->nnodes);\
}\
else if ((f2) == SOLFORMAT_CHROMO)\
{\
	arr_malloc_s(sol->chromo, sol->nnodes, int);\
	if		((f1) == SOLFORMAT_XSTAR)	print_error(ERR_GENERIC_NOT_IMPL, "xstar2chromo");\
	else if ((f1) == SOLFORMAT_SUCC)	result = succ2chromo(sol->succ, sol->chromo, sol->nnodes);\
}

#define REMOVE_FORMAT(sol, f)\
if		(f == SOLFORMAT_XSTAR)	free_s(sol->xstar);\
else if (f == SOLFORMAT_SUCC)	free_s(sol->succ);\
else if (f == SOLFORMAT_CHROMO) free_s(sol->chromo);


// ****************************** SOLUTION MANIPULATION ******************************
/**
* add the edge (i,j) into the solution with the following complexity:
* -	xstar	O(1)
* -	succ	O(1)
* -	chromo	O(n)	inefficient!
*/
void add_edge_solution(Solution* sol, int i, int j, double cost);
char rem_edge_solution(Solution* sol, int i, int j, double cost);

// ******************************* FEASABILITY CHECK *********************************
// check functions
char check_succ(int* succ, int nnodes);
char check_chromo(int* chromo, int nnodes);

// ********************************** GENERAL PURPOSE ********************************
int count_active_edges(int size, double* xstar);
SetOfNodes* compute_set_unvisited_nodes(Solution* sol);

// *********************************** ALGORITHMS ************************************
void build_convex_hull(Graph* g, Solution* sol);
void furthest_nodes_sol(Graph* g, Solution* sol);


/* ***********************************************************************************
*						SOLUTION ITERATOR STRUCTURE DEFINITION
*********************************************************************************** */
typedef struct
{
	Solution* sol;
	int start;		// starting node
	int curr;		// current node in iterator
	int info;		// additional info based on the type of solution
} SolutionIterator;

/* ***********************************************************************************
*							SOLUTION ITERATOR FUNCTIONS
*********************************************************************************** */

void initialize_sol_iterator(SolutionIterator* iter, Solution* sol, int start_node);
inline int curr_node_sol_iterator(SolutionIterator* iter) { return iter->curr; }
char next_node_in_solution(SolutionIterator* iter);

#endif