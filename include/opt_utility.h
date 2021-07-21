#ifndef OPT_UTILITY_H_  

#define OPT_UTILITY_H_

#include "tsp.h"
#include "population.h" 

// 2-opt move
void exec_2opt(int* succ, int a, int c);
void exec_2opt_c(int* chromo, int i, int j);
double move_2opt(int* succ, Graph* g, TabuList* tabu, char allow_worsening);
double move_2opt_c(int* chromo, Graph* g);
double remove_crossings(int* succ, Graph* g);
double remove_crossings_c(int* succ, Graph* g, int howmuch);

void randflip_2opt(int* succ, int nnodes);
void mutation_rms(int* chromo, int nnodes);
void mutation_swap(int* chromo, int nnodes);

void crossover_aex(int* parent1, int* parent2, int* child, int nnodes);
void crossover_pmx(int* parent1, int* parent2, int* child, int nnodes);
void crossover_mix(int* parent1, int* parent2, int* child, int nnodes);

// 
void quicksort_elites(Population* pop, int start, int end);
void merge_pops(Population* pop_dst, Population* pop_other);
void populate(Population* pop, Graph* g);
/**
* perform the whole mating process in the population, from mating, to mutation, to 2-OPT
* refinement, and finally change of generation.
* fitness_emphasis is a number in [0,1] which accounts for how many offsprings should
* come from the best specimen, where 1 means only mating between the "good" division are
* allowed and 0 means that mating is totally random in the population
*/
void new_generation(Population* pop, void (*crossover)(int*, int*, int*, int),
	Graph* g, double mutation_p, double timelim);

#endif