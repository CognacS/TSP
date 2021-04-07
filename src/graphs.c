#include "../graphs.h"

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
