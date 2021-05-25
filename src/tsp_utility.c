#include "../include/tsp_utility.h"

void coord_transform(graph* g)
{
	switch (g->distance_type)
	{
	case EUC_2D:
	case ATT:
		break;
	case GEO:
		geo_transform(g);
		break;
	default:
		print_error_ext(ERR_GENERIC_NOT_IMPL, "distance %d", g->distance_type);
	}
}

double ggetx(graph* g, int idx)
{
	switch (g->distance_type)
	{
	case EUC_2D:
	case ATT:
		return g->xcoord[idx];
	case GEO:
		return g->tr_xcoord[idx];
	default:
		print_error_ext(ERR_GENERIC_NOT_IMPL, "distance %d", g->distance_type);
	}

}
double ggety(graph* g, int idx)
{
	switch (g->distance_type)
	{
	case EUC_2D:
	case ATT:
		return g->ycoord[idx];
	case GEO:
		return g->tr_ycoord[idx];
	default:
		print_error_ext(ERR_GENERIC_NOT_IMPL, "distance %d", g->distance_type);
	}
}

void geo_transform(graph* g)
{
	// allocate arrays
	calloc_s(g->tr_xcoord, g->nnodes, double);
	calloc_s(g->tr_ycoord, g->nnodes, double);

	// transform each (x,y) coordinates from degrees to radians
	for (int i = 0; i < g->nnodes; i++)
	{
		g->tr_xcoord[i] = degrees_to_radians(g->xcoord[i]);
		g->tr_ycoord[i] = degrees_to_radians(g->ycoord[i]);
	}
}

double degrees_to_radians(double coord)
{
	int deg = (int)coord;
	double min = coord - deg;
	double rad = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
	return rad;
}

double dist(int i, int j, graph* g)
{
	double distance = 0;
	switch (g->distance_type)
	{
	case EUC_2D:
		distance = euc_dist(i, j, g);
		break;
	case ATT:
		distance = att_dist(i, j, g);
		break;
	case GEO:
		distance = geo_dist(i, j, g);
		break;
	}

	// if using integer costs, round the distance to
	// closer integer
	if (g->integer_costs) distance = dist_to_int(distance);

	return distance;
}

double euc_dist(int i, int j, graph* g)
{
	// get difference of coordinates
	double dx = g->xcoord[i] - g->xcoord[j];
	double dy = g->ycoord[i] - g->ycoord[j];
	// compute euclidean distance
	double dist = sqrt(dx * dx + dy * dy);
	// return distance
	return dist;
}

double att_dist(int i, int j, graph* g)
{
	// get difference of coordinates
	double dx = g->xcoord[i] - g->xcoord[j];
	double dy = g->ycoord[i] - g->ycoord[j];
	// compute pseudo-euclidean distance
	double dist = sqrt((dx * dx + dy * dy) / 10.0);
	return dist;
}

double geo_dist(int i, int j, graph* g)
{
	double q1 = cos(g->tr_ycoord[i] - g->tr_ycoord[j]);
	double q2 = cos(g->tr_xcoord[i] - g->tr_xcoord[j]);
	double q3 = cos(g->tr_xcoord[i] + g->tr_xcoord[j]);
	int dij = (int)(RRR * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
	return dij;
}

int dist_to_int(double distance)
{
	// return the integer cost as the nearest integer
	int int_dist = (int)(distance + 0.5 - EPSILON);
	return int_dist;
}

void compute_dists(double* dists, graph* g)
{
	int pos = 0;
	for (int i = 0; i < g->nnodes; i++)
	{
		for (int j = i + 1; j < g->nnodes; j++)
		{
			dists[pos++] = dist(i, j, g);
		}
	}
}

double delta_cost(int a, int b, int c, graph* g)
{
	return dist(a, c, g) + dist(c, b, g) - dist(a, b, g);
}

int xpos(int i, int j, int nnodes)
{
	// if indices are equal throw error
	if (i == j) print_error(ERR_INVALID_FUNC_ARGS, "got equal indices in xpos");
	// if indices are inverted return the inverted xpos
	if (i > j) return xpos(j, i, nnodes);
	// compute xpos
	return (int)(i * (nnodes - 3.0 / 2.0) - i * i / 2.0 + j - 1);
}

int xxpos(int i, int j, int nnodes)
{
	// if indices are equal throw error
	if (i == j) print_error(ERR_INVALID_FUNC_ARGS, "got equal indices in xxpos");
	// else return xxpos
	if (i < j) return 2 * xpos(i, j, nnodes);
	else return 2 * xpos(j, i, nnodes) + 1;
}

int upos(int i, int nnodes)
{
	// if index is of first node throw error
	if (i == 0) print_error(ERR_INVALID_FUNC_ARGS, "got node 0 in upos");
	// else return upos
	return i - 1 + (nnodes-1) * nnodes;
}

int ypos(int i, int j, int nnodes)
{
	// if indices are equal throw error
	if (i == j) print_error(ERR_INVALID_FUNC_ARGS, "got equal indices in ypos");
	// else return ypos
	return xxpos(i, j, nnodes) + (nnodes - 1) * nnodes;
}

int is_one(double num)
{
	return num > 0.5;
}

int is_zero_strict(double num)
{
	return num < XSMALL;
}
int is_one_strict(double num)
{
	return num > 1.0 - XSMALL;
}

size_t count_active_edges(size_t size, double* xstar)
{
	size_t active_edges = 0;
	for (size_t i = 0; i < size; i++)
	{
		active_edges += !is_zero_strict(xstar[i]);
	}

	return active_edges;
}


void extract_active_edges(int nnodes, double* xstar, int* ecount, int** elist)
{
	int original_edge_num = (int)(nnodes * (nnodes - 1) / 2);
	// compute edge count
	*ecount = 0;
	for (int i = 0; i < original_edge_num; i++)
	{
		if (xstar[i] > XSMALL) (*ecount)++;
	}
	// allocate elist
	arr_malloc_s(*elist, *ecount, int);
	// fill elist with active edges
	int arr_idx = 0;
	int x_idx;
	for (int i = 0; i < nnodes; i++)
	{
		for (int j = i + 1; j < nnodes; j++)
		{
			x_idx = xpos(i, j, nnodes);
			if (xstar[x_idx] > XSMALL)
			{
				(*elist)[2*arr_idx] = i;
				(*elist)[2*arr_idx+1] = j;
				arr_idx++;
			}
		}
	}
}

void compute_idx(int nnodes, int ecount, int* elist, int** idxlist)
{
	// allocate arrays
	arr_malloc_s(*idxlist, ecount, int);

	// iterate through elist
	for (int i = 0; i < ecount; i++)
	{
		(*idxlist)[i] = xpos(elist[2 * i], elist[2 * i + 1], nnodes);
	}
}

void compute_idx_dst(graph* g, int ecount, int* elist, int** idxlist, unsigned int** dstlist)
{
	// allocate arrays
	arr_malloc_s(*idxlist, ecount, int);
	arr_malloc_s(*dstlist, ecount, unsigned int);

	// iterate through elist
	for (int i = 0; i < ecount; i++)
	{
		(*idxlist)[i] = xpos(elist[2 * i], elist[2 * i + 1], g->nnodes);
		(*dstlist)[i] = (unsigned int)dist(elist[2 * i], elist[2 * i + 1], g);
	}
}

void radix_sort(int ecount, int** idxlist, unsigned int** dstlist)
{

	// define main arrays
	int*			idx_main_arr = *idxlist;
	unsigned int*	dst_main_arr = *dstlist;
	// define aux arrays
	int*			idx_aux_arr = NULL;	arr_malloc_s(idx_aux_arr, ecount, int);
	unsigned int*	dst_aux_arr = NULL;	arr_malloc_s(dst_aux_arr, ecount, unsigned int);
	// define swap arrays
	int*			idx_swap_arr = NULL;
	unsigned int*	dst_swap_arr = NULL;

	// find the max number
	unsigned int maxnum = 0;
	for (int i = 0; i < ecount; i++)
		if (dst_main_arr[i] > maxnum) maxnum = dst_main_arr[i];

	// use radix of 256
	int count[256];
	int idx;

	// iterate through digits of 8 bits until the max number has a digit
	for (int s = 0; maxnum >> s; s += 8)
	{
		// reset counts
		for (int i = 0; i < 256; i++) count[i] = 0;
		// count occurrences of each digit using shift and a mask of 8 bits
		for (int i = 0; i < ecount; i++)
			count[(dst_main_arr[i] >> s) & 0xff]++;

		// compute prefix sums
		for (int i = 1; i < 256; i++)
			count[i] += count[i - 1];

		// build output array iterating in inverse order
		for (int i = ecount - 1; i >= 0; i--)
		{
			// calculate index of interest
			idx = (dst_main_arr[i] >> s) & 0xff;
			count[idx]--;
			// copy value to auxiliary array
			idx_aux_arr[count[idx]] = idx_main_arr[i];
			dst_aux_arr[count[idx]] = dst_main_arr[i];
		}

		// swap pointers of arrays to make aux->main
		idx_swap_arr = idx_main_arr;	dst_swap_arr = dst_main_arr;
		idx_main_arr = idx_aux_arr;		dst_main_arr = dst_aux_arr;
		idx_aux_arr = idx_swap_arr;		dst_aux_arr = dst_swap_arr;
	}

	// free auxiliary arrays as the output arrays are contained in the main arrays
	free(idx_aux_arr);
	free(dst_aux_arr);
	// assign main arrays
	*idxlist = idx_main_arr;
	*dstlist = dst_main_arr;
}


int xstar2succ(double* xstar, int* succ, int nnodes)
{
	// array for visited nodes
	char* visited = NULL;	calloc_s(visited, nnodes, char);
	// identify tour of FEASABLE solution
	int curr = 0; char feasable = 1;
	int active_edges = 0;
	for (int i = 0; i < nnodes-1 && feasable; i++)
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
	memset(xstar, 0, (int)(nnodes * (nnodes - 1) / 2.0) * sizeof(double));

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

int chromo2succ(int* succ, int* chromo, int nnodes)
{
	int active_edges = 0;
	for (int i = 0; i < nnodes-1; i++)
	{
		succ[chromo[i]] = chromo[i+1];
	}
	succ[chromo[nnodes-1]] = chromo[0];

	return active_edges == nnodes;
}


int convert_solution(Solution* sol, solformat format)
{
	// if the format is already correct, return correct
	if (sol->format == format) return 1;

	int result = 1;

	// let f = format
	// let s = sol->format

	// add = in f and not in s = f & ~s
	solformat add_formats = format & ~sol->format;
	// rem = not in f and in s = ~f & s
	solformat rem_formats = ~format & sol->format;
	// use a base format to convert
	solformat base_format = solformat_collapse(sol->format);
	if (base_format == SOLFORMAT_NOFORMAT) print_error(ERR_OPT_SOLFORMAT_NOTSET, NULL);

	// add formats
	solformat needed_format = 1;
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

	solformat notneeded_format = 1;
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


void log_datastruct(void* object, int type, int runlvl, int loglvl)
{
	if (runlvl >= loglvl)
	{
		switch (type)
		{
		case TYPE_GRAPH:
			print_graph((graph*) object);
			break;
		case TYPE_PARAM:
			print_params((params*)object);
			break;
		case TYPE_GLOB:
			print_global_data((global_data*)object);
			break;
		case TYPE_MODEL:
			print_model((model*)object);
			break;
		case TYPE_INST:
			print_instance((instance*)object);
			break;
		default:
			print_warn(WARN_WRONG_DATASTRUCT, NULL);
		}
	}
}

void empty_solution(Solution* sol, solformat format, int nnodes)
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

void clear_solution(Solution* sol)
{
	int nnodes = sol->nnodes;
	solformat format = sol->format;

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
	sol->handle_node = -1;
}

void add_edge_solution(Solution* sol, int i, int j, double cost)
{
	solformat format = sol->format;

	// add edge
	if (format & SOLFORMAT_CHROMO)
	{
		// convert solution to succ and disable chromo
		convert_solution(sol, (format | SOLFORMAT_SUCC ) & ~SOLFORMAT_CHROMO);
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
	solformat format = sol->format;
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

int compute_list_unvisited_nodes(Solution* sol, int* nodelist)
{
	char* visited = NULL;	calloc_s(visited, sol->nnodes, char);

	SolutionIterator iter;
	initialize_sol_iterator(&iter, sol, sol->handle_node);
	// iterate through the tour
	do
	{
		visited[iter.curr] = 1;
	} while (next_node_in_solution(&iter));

	// fill the list with unvisited nodes
	int pos = 0;
	for (int i = 0; i < sol->nnodes; i++)
	{
		if (!visited[i]) nodelist[pos++] = i;
	}
	// return 0 if all nodes are visited
	if (pos == 0) return 0;

	// fill remaining positions with -1
	for (int i = pos; i < sol->nnodes; i++)
	{
		nodelist[i] = -1;
	}
	// return pos as list size
	return pos;

}

void print_directed_sol(graph* g, double* xstar)
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
				log_line_ext(VERBOSITY, LOGLVL_MSG, "%d -> %d with dist %f", curr_node + 1, j + 1, dist(curr_node, j, g));
				curr_node = j;
				break;
			}
		}
	}
}

void print_undirected_sol(graph* g, double* xstar)
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
				log_line_ext(VERBOSITY, LOGLVL_MSG, "%d <-> %d with dist %f", curr_node + 1, j + 1, dist(curr_node, j, g));
				prev_node = curr_node;
				curr_node = j;
				break;
			}
		}
	}
}

void print_xstar(size_t size, double* xstar)
{
	for (size_t i = 0; i < size; i++)
	{
		log_line_ext(VERBOSITY, LOGLVL_PEDANTIC, "Edge %llu: %f", i, xstar[i]);
	}
}