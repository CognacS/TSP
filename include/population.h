#ifndef POPULATION_H_

#define POPULATION_H_

#include "solution.h"

/* ***********************************************************************************
*						POPULATION STRUCTURE DEFINITION
*********************************************************************************** */

typedef struct
{
	int* chromo;
	double cost;

} Specimen;

typedef struct
{
	Specimen* pool;		// pool of solutions currently in the population
	Specimen* champion;	// best solution in the solutions pool

	int pool_size;		// size of the population
	int elite_size;		// quantity of elites in the population
	int nnodes;

} Population;

/* ***********************************************************************************
*								POPULATION FUNCTIONS
*********************************************************************************** */

// ******************************* POPULATION HANDLING *******************************
Population* POP_new(int pool_size, int elite_size, int nnodes);
Population* POP_new_like(Population* pop);
void POP_free(Population* pop);

// ******************************** CHAMPION HANDLING ********************************
inline Specimen* POP_getchampion(Population* pop) { return pop->champion; }

// ***************************** POPULATION'S STATISTICS *****************************
double POP_fitness_avg(Population* pop);
double POP_fitness_best(Population* pop);

// ********************************* SPECIMEN CHOICE *********************************
Specimen* POP_spec_choose_rand(Population* pop); // choose from the entire specimen pool
Specimen* POP_spec_choose_good(Population* pop); // assuming first division comprises best specimen
Specimen* POP_spec_tournament(Population* pop); // choose 3 individuals and pick best

// ************************************ OPERATIONS ***********************************
Specimen* POP_choose(Population* pop);

#endif