#include "../include/graph.h"

/* **************************************************************************************************
*							GRAPH HANDLING
************************************************************************************************** */
void print_graph(Graph* g)
{
	printf("Graph of %d nodes, distances are int=%d, coords:\n", g->nnodes, g->integer_costs);
}

void default_graph(Graph* g)
{
	g->nnodes = -1;
	g->distance_type = -1;
	g->xcoord = NULL;
	g->ycoord = NULL;
	g->tr_xcoord = NULL;
	g->tr_ycoord = NULL;
	g->distance_matrix = NULL;
	g->integer_costs = DEF_INTEGER_COSTS;
}

void free_graph(Graph* g)
{
	free_s(g->xcoord);
	free_s(g->ycoord);
	free_s(g->tr_xcoord);
	free_s(g->tr_ycoord);

	if (g->distance_matrix)
	{
		// free each row of the matrix
		for (int i = 0; i < g->nnodes - 1; i++) free(g->distance_matrix[i]);
	}
	free_s(g->distance_matrix);
}

/* **************************************************************************************************
*							COORDINATES
************************************************************************************************** */
void coord_transform(Graph* g)
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

double ggetx(Graph* g, int node)
{
	switch (g->distance_type)
	{
	case EUC_2D:
	case ATT:
		return g->xcoord[node];
	case GEO:
		return g->tr_xcoord[node];
	default:
		print_error_ext(ERR_GENERIC_NOT_IMPL, "distance %d", g->distance_type);
	}

}
double ggety(Graph* g, int node)
{
	switch (g->distance_type)
	{
	case EUC_2D:
	case ATT:
		return g->ycoord[node];
	case GEO:
		return g->tr_ycoord[node];
	default:
		print_error_ext(ERR_GENERIC_NOT_IMPL, "distance %d", g->distance_type);
	}
}

void geo_transform(Graph* g)
{
	// allocate arrays
	if (g->tr_xcoord == NULL) calloc_s(g->tr_xcoord, g->nnodes, double);
	if (g->tr_ycoord == NULL) calloc_s(g->tr_ycoord, g->nnodes, double);

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

/* **************************************************************************************************
*							DISTANCE
************************************************************************************************** */
double dist(Graph* g, int i, int j)
{
	if (i == j) return 0.0;

	// if the distance matrix is defined
	if (g->distance_matrix != NULL)
	{
		// if i is bigger, swap indices
		if (i > j)
		{
			int c = i; i = j; j = c;
		}
		// return element from the distance matrix
		return g->distance_matrix[i][j-i-1];
	}

	// if not using the distance matrix, compute the distance
	return calc_dist(g, i, j);
}

double calc_dist(Graph* g, int i, int j)
{
	// compute distance based on distance type used
	double distance = 0;
	switch (g->distance_type)
	{
	case EUC_2D:
		distance = euc_dist(g, i, j);
		break;
	case ATT:
		distance = att_dist(g, i, j);
		break;
	case GEO:
		distance = geo_dist(g, i, j);
		break;
	}

	// if using integer costs, round the distance to
	// closer integer
	if (g->integer_costs) distance = dist_to_int(distance);

	return distance;
}

double euc_dist(Graph* g, int i, int j)
{
	// get difference of coordinates
	double dx = g->xcoord[i] - g->xcoord[j];
	double dy = g->ycoord[i] - g->ycoord[j];
	// compute euclidean distance
	double dist = sqrt(dx * dx + dy * dy);
	// return distance
	return dist;
}

double att_dist(Graph* g, int i, int j)
{
	// get difference of coordinates
	double dx = g->xcoord[i] - g->xcoord[j];
	double dy = g->ycoord[i] - g->ycoord[j];
	// compute pseudo-euclidean distance
	double dist = sqrt((dx * dx + dy * dy) / 10.0);
	return dist;
}

double geo_dist(Graph* g, int i, int j)
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
	return nearest_integer(distance);
}

void compute_distance_matrix(Graph* g)
{
	// allocate rows of distance matrix
	arr_malloc_s(g->distance_matrix, g->nnodes - 1, double*);

	for (int i = 0; i < g->nnodes-1; i++)
	{
		// allocate columns of distance matrix
		arr_malloc_s(g->distance_matrix[i], g->nnodes-i-1, double);

		// fill cells of distance matrix
		for (int j = i + 1; j < g->nnodes; j++)
		{
			g->distance_matrix[i][j-i-1] = calc_dist(g, i, j);
		}
	}
}

double delta_cost(Graph* g, int a, int b, int c)
{
	return dist(g, a, c) + dist(g, c, b) - dist(g, a, b);
}

/* **************************************************************************************************
*								GENERAL PURPOSE
************************************************************************************************** */
char counterclockwise(Graph* g, int p, int q, int r)
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

/* **************************************************************************************************
*							INDEXING
************************************************************************************************** */
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

/* **************************************************************************************************
*								ALGORITHMS
************************************************************************************************** */
int find_conncomps_dfs(Graph* g, const double* xstar, int* succ, int* comp, int* ncomp)
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