#include "../include/opt_utility.h"

/* **************************************************************************************************
*						2-OPT MOVE: removes a crossing if it exists
************************************************************************************************** */
double move_2opt(int* succ, Graph* g, char allow_unimproving)
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
	for (int i = 0; i < g->nnodes; i++)
	{
		b = succ[a];
		ab_dist = dist(g, a, b);

		// c starts after b
		c = succ[b];
		// for each remaining edge in the tour
		// apart from (?, a), (a,b), (b, ?)
		for (int j = 0; j < g->nnodes - 3; j++)
		{
			d = succ[c];
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
			c = d;
			d = succ[d];

		}
		// get to the next edge
		a = b;
		b = succ[b];
	}

	// if this is an improving move (delta<0)
	// or non improving moves are allowed
	if (delta_min < 0 || allow_unimproving)
	{
		// remember next of b_min
		int new_next = b_min;
		int new_curr = succ[b_min];
		int new_prev;

		succ[a_min] = c_min;
		succ[b_min] = d_min;

		// reverse path from b_min to c_min
		do
		{
			// swap direction
			new_prev = succ[new_curr];
			succ[new_curr] = new_next;

			// shift
			new_next = new_curr;
			new_curr = new_prev;

		} while (new_next != c_min);


		return delta_min;
	}

	return 0;
}

double remove_crossings(int* succ, Graph* g)
{
	double delta_sum = 0;
	double delta;
	while ((delta = move_2opt(succ, g, 0)) < 0) delta_sum += delta;
	return delta_sum;
}