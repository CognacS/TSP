#include "../include/opt_utility.h"

void exec_2opt(int* succ, int a, int c)
{
	// move to execute:
	// a->b	} => { a->c
	// c->d	} => { b->d
	// then reverse all edges from c to b (done from b to c)

	int b = succ[a], d = succ[c];

	// remember next of b
	int new_next = b;
	int new_curr = succ[b];
	int new_prev;

	succ[a] = c;
	succ[b] = d;

	// reverse path from b to c
	do
	{
		// swap direction
		new_prev = succ[new_curr];
		succ[new_curr] = new_next;

		// shift
		new_next = new_curr;
		new_curr = new_prev;

	} while (new_next != c);
}

void exec_2opt_c(int* chromo, int i, int j)
{
	// invert sequence
	for (int k = 0; k < ((j - i + 1) / 2); k++)
	{
		swap(chromo[i + k], chromo[j - k]);
	}
}

/* **************************************************************************************************
*						2-OPT MOVE: removes a crossing if it exists
************************************************************************************************** */
double move_2opt(int* succ, Graph* g, TabuList* tabu, char allow_worsening)
{
	// define edges ab and cd as those to replace
	// define edges ac and bd as those to insert
	double ab_dist, cd_dist, ac_dist, bd_dist;
	int a, b, c, d;
	int a_min = 0, b_min = 0, c_min = 0, d_min = 0;
	double delta_min = INFINITY;
	double delta_curr;

	a = 0;
	// for each node of the tour
	for (int i = 0; i < g->nnodes - 1; i++)
	{
		b = succ[a];
		// if node is not tabu (if using tabu)
		if (!using_tabu(tabu) || (using_tabu(tabu) && !(TABU_istabu(tabu, a) || TABU_istabu(tabu, b))))
		{
			ab_dist = dist(g, a, b);

			// c starts after b
			c = succ[b];
			// for each remaining edge in the tour
			// apart from (?, a), (a,b), (b, ?)
			// also avoid repeating a comparisons, (a,b);(c,d) == (c,d);(a,b)
			for (int j = 0; (j < (g->nnodes - 3)) && (c != 0); j++)
			{
				d = succ[c];
				if (!using_tabu(tabu) || (using_tabu(tabu) && !(TABU_istabu(tabu, c) || TABU_istabu(tabu, d))))
				{

					cd_dist = dist(g, c, d);
					ac_dist = dist(g, a, c);
					bd_dist = dist(g, b, d);

					delta_curr = (ac_dist + bd_dist) - (ab_dist + cd_dist);
					if (delta_curr < delta_min)
					{
						delta_min = delta_curr;
						a_min = a;
						b_min = b;
						c_min = c;
						d_min = d;
					}
				}
				c = d;
				d = succ[d];

			}
		}
		// get to the next edge
		a = b;
		b = succ[b];
	}

	// if this is an improving move (delta<0)
	// or if using tabu list
	if (delta_min < 0 || allow_worsening)
	{
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] 2-OPT move: (%d,%d);(%d,%d)=>(%d,%d);(%d,%d)",
			a_min, b_min, c_min, d_min, a_min, c_min, b_min, d_min);

		exec_2opt(succ, a_min, c_min);

		// if the move was worsening, make the reverse move tabu
		if (using_tabu(tabu) && delta_min >= 0)
		{
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Added to tabu list node %d", a_min);
			TABU_update(tabu, a_min);
		}

		return delta_min;
	}


	return 0;
}

double move_2opt_c(int* chromo, Graph* g)
{
	// define edges ab and cd as those to replace
	// define edges ac and bd as those to insert
	double ab_dist, cd_dist, ac_dist, bd_dist;
	int a, b, c, d;
	int i_min = 0, j_min = 0;
	double delta_min = INFINITY;
	double delta_curr;

	for (int i = 0; i < g->nnodes - 1; i++)
	{
		a = chromo[i];
		b = chromo[i + 1];
		ab_dist = dist(g, a, b);

		for (int j = i + 2; j < g->nnodes; j++)
		{
			c = chromo[j];
			d = chromo[(j + 1) % g->nnodes];
			cd_dist = dist(g, c, d);
			ac_dist = dist(g, a, c);
			bd_dist = dist(g, b, d);

			delta_curr = (ac_dist + bd_dist) - (ab_dist + cd_dist);
			if (delta_curr < delta_min)
			{
				delta_min = delta_curr;
				i_min = i+1;
				j_min = j;
			}
		}
	}

	if (delta_min < 0)
	{
		exec_2opt_c(chromo, i_min, j_min);
	}

	return delta_min;
}

double remove_crossings(int* succ, Graph* g)
{
	double delta_sum = 0;
	double delta;
	while ((delta = move_2opt(succ, g, NULL, 0)) < 0) delta_sum += delta;
	return delta_sum;
}

double remove_crossings_c(int* chromo, Graph* g, int howmuch)
{
	double delta_sum = 0;
	double delta;
	int times = 0;
	while ((delta = move_2opt_c(chromo, g)) < 0 && (times < howmuch))
	{
		delta_sum += delta; times++;
	}
	return delta_sum;
}

void randflip_2opt(int* succ, int nnodes)
{
	int node1 = rand_int(nnodes);
	int node2 = rand_int(nnodes);

	exec_2opt(succ, node1, node2);
}

void mutation_rms(int* chromo, int nnodes)
{
	// select indices (nodes in sequence)
	int idx1 = rand_int(nnodes);
	int idx2 = rand_int(nnodes);

	if (idx1 > idx2) swap(idx1, idx2);

	// invert sequence
	exec_2opt_c(chromo, idx1, idx2);
}

void mutation_swap(int* chromo, int nnodes)
{
	// select indices (nodes in sequence)
	int idx1 = rand_int(nnodes);
	int idx2 = rand_int(nnodes);

	swap(chromo[idx1], chromo[idx2]);

}

void crossover_aex(int* parent1, int* parent2, int* child, int nnodes)
{
	SetOfNodes* missing_nodes = SETN_new_allnodes(nnodes);
	int* succ1;	arr_malloc_s(succ1, nnodes, int);
	int* succ2;	arr_malloc_s(succ2, nnodes, int);

	chromo2succ(parent1, succ1, nnodes);
	chromo2succ(parent2, succ2, nnodes);

	int* succs[] = { succ1, succ2 };

	int curr = parent1[0];
	char turn = 0;
	for (int i = 0; i < nnodes; i++)
	{
		// if the node is already in the child, extract a random missing node
		if (!SETN_exists(missing_nodes, curr))
		{
			curr = SETN_rand_node(missing_nodes);
		}
		child[i] = curr;
		SETN_remove(missing_nodes, curr);
		curr = succs[turn][curr];
		turn = (turn + 1) % 2;
	}
	SETN_free(missing_nodes);
	free(succ1); free(succ2);
}

void crossover_pmx(int* parent1, int* parent2, int* child, int nnodes)
{
	char* visited = NULL; calloc_s(visited, nnodes, char);

	// produce crossover points
	int px1 = rand_int(nnodes);
	int px2 = rand_int(nnodes);
	if (px1 > px2) swap(px1, px2);

	int* index2;	arr_malloc_s(index2, nnodes, int);
	chromo2index(parent2, index2, nnodes);

	for (int i = 0; i < nnodes; i++) child[i] = -1;

	// copy slice
	for (int i = px1; i <= px2; i++)
	{
		child[i] = parent1[i];
		visited[child[i]] = 1;
	}

	// map same slice in parent2 to child using indices from parent1
	int idx;
	for (int i = px1; i <= px2; i++)
	{
		if (!visited[parent2[i]])
		{
			idx = i;
			while (child[idx] != -1)
			{
				idx = index2[parent1[idx]];
			}
			child[idx] = parent2[i];
		}
	}

	for (int i = 0; i < nnodes; i++)
	{
		if (child[i] == -1) child[i] = parent2[i];
	}

	free(visited); free(index2);
}

void crossover_mix(int* parent1, int* parent2, int* child, int nnodes)
{
	// 50:50 chances
	if (rand_int(2) == 1)
		crossover_aex(parent1, parent2, child, nnodes);
	else
		crossover_pmx(parent1, parent2, child, nnodes);
}

void quicksort_elites(Population* pop, int start, int end)
{
	int i, j;
	int pivot = rand_int_range(start, end);
	double pivot_cost = pop->pool[pivot].cost;
	Specimen t;
	if (start < end)
	{
		i = start - 1;
		j = end + 1;

		while (1)
		{
			do { i++; } while (pop->pool[i].cost < pivot_cost);
			do { j--; } while (pop->pool[j].cost > pivot_cost);

			if (i >= j) break;
			else swap_ext(pop->pool[i], pop->pool[j], t);
		}

		quicksort_elites(pop, start, j);
		quicksort_elites(pop, j + 1, end);
	}
}

void merge_pops(Population* pop_dst, Population* pop_other)
{
	Specimen t;
	int elite_size = pop_dst->elite_size;

	// do a full reorder of the offsprings
	quicksort_elites(pop_other, 0, pop_other->pool_size-1);
	pop_other->champion = pop_other->pool;

	// replace best elements of other population with worst elements
	// of the destination (on the elite side)
	int qty_replace;
	for (qty_replace = 0; qty_replace < pop_other->nnodes; qty_replace++)
	{
		if (pop_other->pool[qty_replace].cost >= pop_dst->pool[elite_size - qty_replace - 1].cost) break;
	}
	int dst_start = elite_size - qty_replace;
	for (int i = 0; i < qty_replace; i++)
	{	// insert better elements in destination and worse elements in other
		swap_ext(pop_other->pool[i], pop_dst->pool[dst_start + i], t);
	}

	// final sort of elites
	quicksort_elites(pop_dst, 0, elite_size - 1);
	pop_dst->champion = pop_dst->pool;

	for (int i = elite_size; i < pop_dst->pool_size; i++)
	{
		// 50/50 chance of retaining or inserting new element on the random side of pop
		if (rand_int(3) == 0) {
			int pos = pop_other->elite_size + rand_int(pop_other->pool_size - pop_other->elite_size);
			swap_ext(pop_other->pool[i], pop_dst->pool[pos], t);
		}
	}

}

void populate(Population* pop, Graph* g)
{
	int nnodes = pop->nnodes;
	int* remaining_nodes = NULL;	arr_malloc_s(remaining_nodes, nnodes, int);

	for (int i = 0; i < pop->pool_size; i++)
	{
		// refill remaining nodes
		for (int i = 0; i < nnodes; i++) remaining_nodes[i] = i;

		// allocate specimen
		Specimen* spec = pop->pool + i;
		arr_malloc_s(spec->chromo, nnodes, int);

		// generate specimen
		int idx;
		for (int i = 0; i < nnodes; i++)
		{
			// select node and place
			idx = rand_int(nnodes - i);
			spec->chromo[i] = remaining_nodes[idx];
			// place selected at the bottom
			remaining_nodes[idx] = remaining_nodes[nnodes - i - 1];
		}

		// check correctness
		log_line_ext(VERBOSITY, LOGLVL_PEDANTIC, "%d - CORRECTNESS AFTER GENERATION: %d", i, check_chromo(spec->chromo, pop->nnodes));

		// compute cost of specimen
		spec->cost = compute_cost(spec->chromo, g);


		// nominate champion if needed
		if (spec->cost < POP_fitness_best(pop))
		{
			pop->champion = spec;
		}

		log_line_ext(VERBOSITY, LOGLVL_EXTINFO, "[INFO+] %d/%d done", i+1, pop->pool_size);
	}

	// check that the best specimen from the other population is
	// on the elite side of the population
	size_t champ_pos = pop->champion - pop->pool;
	Specimen t;
	if (champ_pos >= pop->elite_size)
	{
		swap_ext(pop->pool[0], *pop->champion, t);
		pop->champion = pop->pool;
	}
	// final full sort
	quicksort_elites(pop, 0, pop->pool_size-1);
	pop->champion = pop->pool;

	free(remaining_nodes);
}

/* **************************************************************************************************
*				NEW GENERATION: mating->mutation->2-OPT->change of generation
************************************************************************************************** */
void new_generation(Population* pop, void (*crossover)(int*, int*, int*, int),
	Graph* g, double mutation_p, double timelim)
{
	double start_time = second();

	// generate a new population of offsprings
	Population* pop_offsprings = POP_new_like(pop);

	// produce offsprings
	char time_expired = 0;
	for (int i = 0; i < pop->pool_size && !time_expired; i++)
	{
		// choose parents
		Specimen* parent1 = POP_choose(pop);
		Specimen* parent2 = POP_choose(pop);
		Specimen* child = pop_offsprings->pool + i;
		arr_malloc_s(child->chromo, pop->nnodes, int);

		// crossover
		crossover(parent1->chromo, parent2->chromo, child->chromo, pop->nnodes);

		// mutate
		if (random() < mutation_p)
			mutation_rms(child->chromo, pop->nnodes);

		if (random() < 0.001)
		{
			double delta = remove_crossings_c(child->chromo, g, 10);
		}

		// check correctness
		log_line_ext(VERBOSITY, LOGLVL_PEDANTIC, "CORRECTNESS AFTER MUTATION: %d", check_chromo(child->chromo, pop->nnodes));

		// compute cost of child
		child->cost = compute_cost(child->chromo, g);

		// nominate champion if needed
		if (child->cost < POP_fitness_best(pop_offsprings))
			pop_offsprings->champion = child;

		time_expired = second() - start_time >= timelim;

		log_line_ext(VERBOSITY, LOGLVL_EXTINFO, "[INFO+] %d/%d done", i + 1, pop->pool_size);
	}

	if (!time_expired)
	{
		// merge the two populations
		merge_pops(pop, pop_offsprings);
	}
	// kill remaining elements
	POP_free(pop_offsprings);
}