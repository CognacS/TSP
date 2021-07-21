#include "../include/population.h"

/* **************************************************************************************************
*							POPULATION HANDLING
************************************************************************************************** */

Population* POP_new(int pool_size, int elite_size, int nnodes)
{
	// allocate structure
	Population* pop = NULL; malloc_s(pop, Population);

	// set sizes
	pop->pool_size = pool_size;
	pop->elite_size = elite_size;
	// allocate solution pool as NULL pointers
	calloc_s(pop->pool, pool_size, Specimen);
	// set champion
	pop->champion = NULL;
	// set nnodes
	pop->nnodes = nnodes;

	return pop;
}

Population* POP_new_like(Population* pop)
{
	return POP_new(pop->pool_size, pop->elite_size, pop->nnodes);
}

void POP_free(Population* pop)
{
	// free all solutions
	for (int i = 0; i < pop->pool_size; i++)
	{
		free(pop->pool[i].chromo);
	}

	// free pool
	free(pop->pool);
	// free population
	free(pop);
}

/* **************************************************************************************************
*							CHAMPION HANDLING
************************************************************************************************** */


/* **************************************************************************************************
*							POPULATION'S STATISTICS
************************************************************************************************** */

double POP_fitness_avg(Population* pop)
{
	double sum = 0;
	for (int i = 0; i < pop->pool_size; i++)
	{
		sum += pop->pool[i].cost;
	}
	return sum / pop->pool_size;
}

double POP_fitness_best(Population* pop)
{
	if (pop->champion != NULL)	return pop->champion->cost;
	else return INFINITY;
}

double POP_entropy(Population* pop)
{
	int nnodes = pop->nnodes;
	int edges = complete_graph_edges(nnodes);
	int* freqs = NULL; calloc_s(freqs, edges, int);
	
	for (int i = 0; i < pop->pool_size; i++)
	{
		int* att_chromo = pop->pool[i].chromo;
		for (int j = 0; j <nnodes-1; j++)
		{
			freqs[xpos(att_chromo[j], att_chromo[j + 1], nnodes)]++;
		}
		freqs[xpos(att_chromo[0], att_chromo[nnodes-1], nnodes)]++;
	}

	double entropy = 0;
	double prob;
	for (int i = 0; i < edges; i++)
	{
		prob = (double)freqs[i] / pop->pool_size;
		if (prob > 0) entropy -= prob * log2(prob);
	}
	free(freqs);
	return entropy;
}

/* **************************************************************************************************
*							SPECIMEN CHOICE
************************************************************************************************** */

Specimen* POP_spec_choose_rand(Population* pop)
{
	int idx = rand_int(pop->pool_size);
	return pop->pool + idx;
}

Specimen* POP_spec_choose_good(Population* pop)
{
	int idx = rand_int(pop->elite_size);
	return pop->pool + idx;
}

Specimen* POP_spec_tournament(Population* pop)
{
	Specimen* spec1 = POP_spec_choose_rand(pop);
	Specimen* spec2 = POP_spec_choose_rand(pop);
	Specimen* spec3 = POP_spec_choose_rand(pop);

	Specimen* spec12 = spec1->cost < spec2->cost ? spec1 : spec2;

	return spec12->cost < spec3->cost ? spec12 : spec3;
}

/* **************************************************************************************************
*							OPERATIONS
************************************************************************************************** */

Specimen* POP_choose(Population* pop)
{
	//return POP_spec_choose_rand(pop);
	return POP_spec_tournament(pop);
}

