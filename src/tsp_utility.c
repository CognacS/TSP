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
				break;
			}
		}
	}
	// set last one's successor as the start 0 if feasable!
	feasable = feasable && is_one(xstar[xpos(curr, 0, nnodes)]);
	if (feasable) succ[curr] = 0;

	free(visited);
	return feasable;
}

int succ2xstar(int* succ, double* xstar, int nnodes)
{
	// array for visited nodes
	char* visited = NULL;	calloc_s(visited, nnodes, char);
	// activate edges of FEASABLE solution
	int curr = 0; char feasable = 1;
	for (int i = 0; i < nnodes && feasable; i++)
	{
		xstar[xpos(curr, succ[curr], nnodes)] = 1.0;
		feasable = !visited[curr];
		visited[curr] = 1;
		curr = succ[curr];
	}

	free(visited);
	return feasable;
}


int convert_solution(Solution* sol, solformat format)
{
	// if the format is already correct, return correct
	if (sol->format == format) return 1;

	int result = 1;
	// if the needed format is not available, produce it
	if (sol->format != SOLFORMAT_BOTH)
	{
		char need_succ = format == SOLFORMAT_SUCC || (format == SOLFORMAT_BOTH && sol->format == SOLFORMAT_XSTAR);
		char need_xstar = format == SOLFORMAT_XSTAR || (format == SOLFORMAT_BOTH && sol->format == SOLFORMAT_SUCC);

		if (need_succ) {
			arr_malloc_s(sol->succ, sol->nnodes, int);
			result = xstar2succ(sol->xstar, sol->succ, sol->nnodes);
		}
		if (need_xstar) {
			arr_malloc_s(sol->xstar, (int)(sol->nnodes*(sol->nnodes-1)/2), double);
			result = succ2xstar(sol->succ, sol->xstar, sol->nnodes);
		}

	}
	// free the format not needed
	if		(format == SOLFORMAT_XSTAR) free_s(sol->succ);
	else if (format == SOLFORMAT_SUCC)	free_s(sol->xstar);
	sol->format = format;

	return result;
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

	if (format == SOLFORMAT_SUCC || format == SOLFORMAT_BOTH)
		arr_malloc_s(sol->succ, sol->nnodes, int);

	if (format == SOLFORMAT_XSTAR || format == SOLFORMAT_BOTH)
		arr_malloc_s(sol->xstar, (int)(sol->nnodes * (sol->nnodes - 1) / 2.0), double);

	// use a big number for cost
	sol->cost = 1e100;
	

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