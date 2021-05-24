#include "../include/graphs.h"

/* **************************************************************************************************
*						CONNECTED COMPONENTS FINDER FUNCTION (using DFS)
************************************************************************************************** */
int find_conncomps_dfs(graph* g, const double* xstar, int* succ, int* comp, int* ncomp)
{
	*ncomp = 0;
	// setup succ and comp arrays
	for (int i = 0; i < g->nnodes; i++)
	{
		succ[i] = -1;
		comp[i] = -1;
	}

	// start DFS search
	for (int start = 0; start < g->nnodes; start++)
	{
		if (comp[start] >= 0) continue;  // node "start" was already visited, just skip it

		// a new component is found
		int i = start;
		int done = 0;
		while (!done)  // go and visit the current component
		{
			comp[i] = *ncomp;
			done = 1;
			for (int j = 0; j < g->nnodes; j++)
			{
				if (i != j && xstar[xpos(i, j, g->nnodes)] > 0.5 && comp[j] == -1) // the edge [i,j] is selected in xstar and j was not visited before 
				{
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}
		succ[i] = start;  // last arc to close the cycle

		// go to the next component...
		(*ncomp)++;
	}
	
	return *ncomp;
}

char counterclockwise(graph* g, int p, int q, int r)
{
	double
		a = ggetx(g, p),
		b = ggety(g, p),
		c = ggetx(g, q),
		d = ggety(g, q),
		e = ggetx(g, r),
		f = ggety(g, r);
	return (f - b) * (c - a) > (d - b) * (e - a);
}

void build_convex_hull(graph* g, Solution* sol)
{
	// find the leftmost node
	int l = 0;
	for (int i = 1; i < g->nnodes; i++)
	{
		if (g->xcoord[i] < g->xcoord[l])
			l = i;
	}

	sol->handle_node = l;
	// iterate through each node and wrap
	int p = l, q;
	do
	{
		// select a different node from p
		q = p + 1;
		// iterate through every node
		for (int i = 0; i < g->nnodes; i++)
		{
			// if node i is more to the right wrt q looking from p: select q=i
			if (counterclockwise(g, p, i, q))
				q = i;
		}
		add_edge_solution(sol, p, q, dist(p, q, g));
		p = q;
	} while (p != l);

}

void furthest_nodes_sol(graph* g, Solution* sol)
{
	int a_idx = 0, b_idx = 0;
	double max_dist = 0.0, curr_dist;
	for (int i = 0; i < g->nnodes; i++)
	{
		for (int j = i + 1; j < g->nnodes; j++)
		{
			curr_dist = dist(i, j, g);
			if (curr_dist > max_dist)
			{
				max_dist = curr_dist;
				a_idx = i;
				b_idx = j;
			}
		}
	}
	// make cycle
	add_edge_solution(sol, a_idx, b_idx, max_dist);
	add_edge_solution(sol, b_idx, a_idx, max_dist);

	sol->handle_node = a_idx;
}