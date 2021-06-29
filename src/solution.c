#include "../include/solution.h"

/* **************************************************************************************************
*							SOLUTION HANDLING
************************************************************************************************** */
void print_solution(Solution* sol)
{
	int nnodes = sol->nnodes;
	SolFormat format = sol->format;

	// momentarily collapse format
	sol->format = solformat_collapse(format);

	// initialize iterator from node 0
	SolutionIterator iter;
	initialize_sol_iterator(&iter, sol, 0);

	// iterate the solution and print each node
	int node;
	do
	{
		node = curr_node_sol_iterator(&iter);
		printf("%d\n", node);

	} while (next_node_in_solution(&iter));

}

void create_solution(Solution* sol, SolFormat format, int nnodes)
{
	sol->nnodes = nnodes;
	sol->format = format;

	if (format & SOLFORMAT_SUCC)
		arr_malloc_s(sol->succ, sol->nnodes, int);
	else sol->succ = NULL;

	if (format & SOLFORMAT_XSTAR)
		arr_malloc_s(sol->xstar, complete_graph_edges(nnodes), double);
	else sol->xstar = NULL;

	if (format & SOLFORMAT_CHROMO)
		arr_malloc_s(sol->chromo, sol->nnodes, int);
	else sol->chromo = NULL;

	// use a big number for cost
	sol->cost = INFINITY;
	// set optimal
	sol->optimal = 0;
}

void erase_solution(Solution* sol)
{
	int nnodes = sol->nnodes;
	SolFormat format = sol->format;

	if (format & SOLFORMAT_SUCC)
		for (int i = 0; i < nnodes; i++) sol->succ[i] = -1;

	if (format & SOLFORMAT_XSTAR)
		for (int i = 0; i < complete_graph_edges(nnodes); i++) sol->xstar[i] = 0.0;

	if (format & SOLFORMAT_CHROMO)
		for (int i = 0; i < nnodes; i++) sol->chromo[i] = -1;

	// reset cost
	sol->cost = 0;
	// set optimal
	sol->optimal = 0;
	// set handle node
	sol->handle_node = -1;
}

void shallow_copy_solution(Solution* dst, Solution* src)
{
	// deep copy params
	dst->cost = src->cost;
	dst->format = src->format;
	dst->nnodes = src->nnodes;
	dst->optimal = src->optimal;
	dst->handle_node = src->handle_node;

	// free all that was allocated
	if (dst->succ != NULL) free(dst->succ);
	if (dst->xstar != NULL) free(dst->xstar);
	if (dst->chromo != NULL) free(dst->chromo);
	// shallow copy arrays
	dst->succ = src->succ;
	dst->xstar = src->xstar;
	dst->chromo = src->chromo;

}

void deep_copy_solution(Solution* dst, Solution* src)
{
	// deep copy params
	dst->cost = src->cost;
	dst->format = src->format;
	dst->nnodes = src->nnodes;
	dst->optimal = src->optimal;
	dst->handle_node = src->handle_node;

	int xstar_size = complete_graph_edges(src->nnodes);
	// free all that was allocated
	if (dst->succ != NULL) DEEP_COPY_ARR(dst->succ, src->succ, src->nnodes);
	if (dst->xstar != NULL) DEEP_COPY_ARR(dst->xstar, src->xstar, xstar_size);
	if (dst->chromo != NULL) DEEP_COPY_ARR(dst->chromo, src->chromo, src->nnodes);

}

void free_solution(Solution* sol)
{
	free(sol->xstar);
	free(sol->succ);
	free(sol->chromo);
}

void print_directed_sol(Graph* g, double* xstar)
{
	int curr_node = 0;
	for (int i = 0; i < g->nnodes; i++)
	{
		for (int j = 0; j < g->nnodes; j++)
		{
			if (curr_node == j) continue;
			// test if value is 0 or 1
			if (is_one(xstar[xxpos(curr_node, j, g->nnodes)]))
			{
				log_line_ext(VERBOSITY, LOGLVL_MSG, "%d -> %d with dist %f", curr_node + 1, j + 1, dist(g, curr_node, j));
				curr_node = j;
				break;
			}
		}
	}
}

void print_undirected_sol(Graph* g, double* xstar)
{
	int curr_node = 0;
	int prev_node = 0;
	for (int i = 0; i < g->nnodes; i++)
	{
		for (int j = 0; j < g->nnodes; j++)
		{
			if (curr_node == j || prev_node == j) continue;
			// test if value is 0 or 1
			if (is_one(xstar[xpos(curr_node, j, g->nnodes)]))
			{
				log_line_ext(VERBOSITY, LOGLVL_MSG, "%d <-> %d with dist %f", curr_node + 1, j + 1, dist(g, curr_node, j));
				prev_node = curr_node;
				curr_node = j;
				break;
			}
		}
	}
}

void print_xstar(int size, double* xstar)
{
	for (int i = 0; i < size; i++)
	{
		log_line_ext(VERBOSITY, LOGLVL_PEDANTIC, "Edge %d: %f", i, xstar[i]);
	}
}

/* **************************************************************************************************
*							SOLUTION FORMAT FUNCTIONS
************************************************************************************************** */

int solformat2num(SolFormat format)
{
	int shift = 0;
	for (shift = 0; shift < NUM_SOLFORMATS && !(format & 1); shift++, format = format >> 1);
	return shift;
}

SolFormat solformat_collapse(SolFormat format)
{
	return  (1 << solformat2num(format)) & SOLFORMAT_ALLFORMAT;
}

/* **************************************************************************************************
*							COMPUTE COSTS
************************************************************************************************** */
double compute_cost(int* chromo, Graph* g)
{
	double cost = 0;
	for (int i = 0; i < g->nnodes - 1; i++)
	{
		cost += dist(g, chromo[i], chromo[i + 1]);
	}
	cost += dist(g, chromo[g->nnodes-1], chromo[0]);
	return cost;
}

/* **************************************************************************************************
*							CONVERSION
************************************************************************************************** */
int xstar2succ(double* xstar, int* succ, int nnodes)
{
	// array for visited nodes
	char* visited = NULL;	calloc_s(visited, nnodes, char);
	// identify tour of FEASABLE solution
	int curr = 0; char feasable = 1;
	int active_edges = 0;
	for (int i = 0; i < nnodes - 1 && feasable; i++)
	{
		visited[curr] = 1;
		int h;
		feasable = 0;
		for (h = 0; h < nnodes; h++)
		{
			// if node h form an edge with curr, and it is not the previous one
			if (curr != h && !visited[h] && is_one(xstar[xpos(curr, h, nnodes)]))
			{
				succ[curr] = h;
				curr = h;
				// feasable if a node h successor of curr was found
				feasable = 1;
				active_edges++;
				break;
			}
		}
	}
	// set last one's successor as the start 0 if feasable!
	feasable = feasable && is_one(xstar[xpos(curr, 0, nnodes)]);
	if (feasable) succ[curr] = 0;

	free(visited);
	return feasable && active_edges == nnodes;
}

int succ2xstar(int* succ, double* xstar, int nnodes)
{
	// clear xstar
	for (int i = 0; i < complete_graph_edges(nnodes); i++)	xstar[i] = 0.0;

	// array for visited nodes
	char* visited = NULL;	calloc_s(visited, nnodes, char);
	// activate edges of FEASABLE solution
	int curr = 0; char feasable = 1;
	int active_edges = 0;
	for (int i = 0; i < nnodes && feasable; i++)
	{
		if (succ[curr] >= 0)
		{
			xstar[xpos(curr, succ[curr], nnodes)] = 1.0;
			feasable = !visited[curr];
			visited[curr] = 1;
			curr = succ[curr];
			active_edges++;
		}
	}

	free(visited);
	return feasable && active_edges == nnodes;
}

int succ2chromo(int* succ, int* chromo, int nnodes)
{
	int curr = 0;
	int active_edges = 0;
	for (int i = 0; i < nnodes; i++)
	{
		if (succ[curr] >= 0)
		{
			chromo[i] = curr;
			curr = succ[curr];
			active_edges++;
		}
	}

	return active_edges == nnodes;
}

int chromo2succ(int* chromo, int* succ, int nnodes)
{
	int active_edges = 0;
	int node = 0;
	for (int i = 0; i < nnodes - 1; i++)
	{
		succ[chromo[i]] = chromo[i + 1];
	}
	succ[chromo[nnodes - 1]] = chromo[0];

	return active_edges == nnodes;
}

int chromo2index(int* chromo, int* index, int nnodes)
{
	for (int i = 0; i < nnodes; i++)
	{
		index[chromo[i]] = i;
	}
}


int convert_solution(Solution* sol, SolFormat format)
{
	// if the format is already correct, return correct
	if (sol->format == format) return 1;

	int result = 1;

	// let f = format
	// let s = sol->format

	// add = in f and not in s = f & ~s
	SolFormat add_formats = format & ~sol->format;
	// rem = not in f and in s = ~f & s
	SolFormat rem_formats = ~format & sol->format;
	// use a base format to convert
	SolFormat base_format = solformat_collapse(sol->format);
	if (base_format == SOLFORMAT_NOFORMAT) print_error(ERR_OPT_SOLFORMAT_NOTSET, NULL);

	// add formats
	SolFormat needed_format = 1;
	int feasable = 1;
	for (int i = 0; i < NUM_SOLFORMATS && feasable; i++)
	{
		if (add_formats & needed_format)
		{
			// here, surely needed_format != base_format
			// then apply allocation and transformation
			FORMAT2FORMAT(sol, base_format, needed_format, feasable);
		}
		needed_format = needed_format << 1;
	}

	SolFormat notneeded_format = 1;
	for (int i = 0; i < NUM_SOLFORMATS && feasable; i++)
	{
		if (rem_formats & notneeded_format)
		{
			// here, surely notneeded_format is not in "format"
			// then apply removal
			REMOVE_FORMAT(sol, notneeded_format);
		}
		notneeded_format = notneeded_format << 1;
	}

	sol->format = format;

	return result;
}

/* **************************************************************************************************
*							SOLUTION MANIPULATION
************************************************************************************************** */
void add_edge_solution(Solution* sol, int i, int j, double cost)
{
	SolFormat format = sol->format;

	// add edge
	if (format & SOLFORMAT_CHROMO)
	{
		// convert solution to succ and disable chromo
		convert_solution(sol, (format | SOLFORMAT_SUCC) & ~SOLFORMAT_CHROMO);
		// add the edge
		add_edge_solution(sol, i, j, cost);
		// convert back to the original format
		convert_solution(sol, format);
	}
	else
	{
		if (format & SOLFORMAT_SUCC) sol->succ[i] = j;
		if (format & SOLFORMAT_XSTAR) sol->xstar[xpos(i, j, sol->nnodes)] += 1.0;
	}

	// add cost
	sol->cost += cost;
}

char rem_edge_solution(Solution* sol, int i, int j, double cost)
{
	SolFormat format = sol->format;
	char exists = 1;

	if (format & SOLFORMAT_SUCC)
	{
		if (sol->succ[i] == j) sol->succ[i] = -1;
		else if (sol->succ[j] == i) sol->succ[j] = -1;
		else exists = 0;
	}
	if (format & SOLFORMAT_XSTAR)
	{
		int var_idx = xpos(i, j, sol->nnodes);
		exists = is_one(sol->xstar[var_idx]);
		if (exists) sol->xstar[var_idx] -= 1.0;
	}
	if (format & SOLFORMAT_CHROMO)
		print_error(ERR_GENERIC_NOT_IMPL, "rem_edge_solution for SOLFORMAT_CHROMO");

	if (exists)
	{
		sol->cost -= cost;
		sol->optimal = 0;
	}

	return exists;
}

/* **************************************************************************************************
*							FEASABILITY CHECK
************************************************************************************************** */
char check_succ(int* succ, int nnodes)
{
	char correct = 1;
	int curr = 0;
	for (int i = 0; (i < nnodes) && correct; i++)
	{
		correct = (succ[i] >= 0) && (succ[i] < nnodes);
		if (correct) curr = succ[curr];
	}
	return correct && curr == 0;
}

char check_chromo(int* chromo, int nnodes)
{
	char correct = 1;
	char* visited = NULL; calloc_s(visited, nnodes, char);
	char err1, err2;
	// check repetitions
	for (int i = 0; (i < nnodes) && correct; i++)
	{
		err1 = visited[chromo[i]];
		err2 = chromo[i] < 0 || chromo[i] >= nnodes;
		if (err1) log_line_ext(VERBOSITY, LOGLVL_WARN, "[WARN] node %d at pos %d already visited!", chromo[i], i);
		if (err2) log_line_ext(VERBOSITY, LOGLVL_WARN, "[WARN] node %d at pos %dout of bound!", chromo[i], i);

		correct = !(err1 || err2);
		visited[chromo[i]] = 1;
	}
	// check tour passing through all nodes
	for (int i = 0; (i < nnodes) && correct; i++)
	{
		correct = visited[chromo[i]];
		if (!correct) log_line_ext(VERBOSITY, LOGLVL_WARN, "[WARN] node %d at pos %d was not visited!", chromo[i], i);
	}

	free(visited);
	return correct;
}

/* **************************************************************************************************
*							GENERAL PURPOSE
************************************************************************************************** */
int count_active_edges(int size, double* xstar)
{
	int active_edges = 0;
	for (int i = 0; i < size; i++)
	{
		active_edges += !is_zero_strict(xstar[i]);
	}

	return active_edges;
}

SetOfNodes* compute_set_unvisited_nodes(Solution* sol)
{
	// start with a full set of nodes
	SetOfNodes* set = SETN_new_allnodes(sol->nnodes);

	SolutionIterator iter;
	initialize_sol_iterator(&iter, sol, sol->handle_node);
	// iterate through the tour and remove all nodes found
	do
	{
		SETN_remove(set, iter.curr);
	} while (next_node_in_solution(&iter));

	// return pos as list size
	return set;

}

/* **************************************************************************************************
*							SOLUTION ALGORITHMS
************************************************************************************************** */
void build_convex_hull(Graph* g, Solution* sol)
{
	// find the leftmost node
	int l = 0;
	for (int i = 1; i < g->nnodes; i++)
	{
		if (g->xcoord[i] < g->xcoord[l])
			l = i;
	}

	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Leftmost node: %d", l);
	sol->handle_node = l;
	// iterate through each node and wrap
	int p = l, q;
	do
	{
		// select a different node from p
		q = (p + 1) % g->nnodes;
		// iterate through every node
		for (int i = 0; i < g->nnodes; i++)
		{
			// if node i is more to the right wrt q looking from p: select q=i
			if (counterclockwise(g, p, i, q))
				q = i;
		}
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Adding edge %d->%d", p, q);
		add_edge_solution(sol, p, q, dist(g, p, q));
		p = q;
	} while (p != l);

}

void furthest_nodes_sol(Graph* g, Solution* sol)
{
	int a_idx = 0, b_idx = 0;
	double max_dist = 0.0, curr_dist;
	for (int i = 0; i < g->nnodes; i++)
	{
		for (int j = i + 1; j < g->nnodes; j++)
		{
			curr_dist = dist(g, i, j);
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

/* **************************************************************************************************
*							SOLUTION ITERATOR FUNCTIONS
************************************************************************************************** */
void initialize_sol_iterator(SolutionIterator* iter, Solution* sol, int start_node)
{
	iter->sol = sol;
	iter->curr = start_node;
	if (sol->format & SOLFORMAT_CHROMO)	// search for starting node
		for (iter->info = 0; iter->info < sol->nnodes && sol->chromo[iter->info] != start_node; iter->info++);
	// else just start from 0
	else iter->info = 0;

	iter->start = start_node;
}

char next_node_in_solution(SolutionIterator* iter)
{
	Solution* sol = iter->sol;

	// if there is next, produce next node
	if (sol->format & SOLFORMAT_SUCC)
	{
		// get succ to curr
		iter->curr = sol->succ[iter->curr];
	}
	else if (sol->format & SOLFORMAT_CHROMO)
	{
		// rotate position
		iter->info = (iter->info + 1) % sol->nnodes;
		// get next in position
		iter->curr = sol->chromo[iter->info];
	}
	else if (sol->format & SOLFORMAT_XSTAR)
	{
		// find active edge in xstar
		int i;
		for (i = 0; i < sol->nnodes; i++)
		{
			// if it is not the previous and the edge is active
			if (i != iter->info && is_one(sol->xstar[xpos(iter->curr, i, sol->nnodes)]))
			{
				iter->info = iter->curr;
				iter->curr = i;
				break;
			}
		}
		if (i == sol->nnodes) iter->curr = -1;
	}

	// return a next if not back at start and if not a dead end
	return iter->curr != iter->start && iter->curr != -1;
}