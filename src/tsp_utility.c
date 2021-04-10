#include "../tsp_utility.h"

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