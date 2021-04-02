#include "../graph_plot.h"

void plot_tsp_solution_undirected(graph* g, double* xstar)
{

	h_GPC_Plot* plot_handle = gpc_init_xy(
		"TSP solution plot",
		"X coords",
		"Y coords",
		GPC_AUTO_SCALE,
		GPC_KEY_DISABLE
	);

	// *************** set gnuplot parameters ***************
	// set window size
	fprintf(plot_handle->pipe, "set term qt size 1000, 800 position 0, 0\n");
	// set relative circle radius
	fprintf(plot_handle->pipe, "set style circle radius character 1.1\n");
	// set offset
	fprintf(plot_handle->pipe, "set offsets graph 0.01, graph 0.01, graph 0.01, graph 0.01\n");
	// set x label
	fprintf(plot_handle->pipe, "set xlabel offset screen 0, 0.07\n");
	// set y label
	fprintf(plot_handle->pipe, "set ylabel offset screen 0.07, 0\n");
	

	// *************** start data block ***************
	fprintf(plot_handle->pipe, "$data << EOD\n");

	for (int i = 0; i < g->nnodes; i++)
	{
		for (int j = i + 1; j < g->nnodes; j++)
		{
			if (is_one(xstar[xpos(i, j, g->nnodes)]))
			{
				// print edge nodes
				fprintf(plot_handle->pipe, "%f %f %d\n", g->xcoord[i], g->ycoord[i], i+1);
				fprintf(plot_handle->pipe, "%f %f %d\n", g->xcoord[j], g->ycoord[j], j+1);
				// print edge separator
				fprintf(plot_handle->pipe, "\n");
			}
		}
	}

	// *************** end data block ***************
	fprintf(plot_handle->pipe, "EOD\n");

	// *************** plot data block ***************
	fprintf(plot_handle->pipe,
		"plot "
		"$data using 1:2		with lines lc rgb \"red\" lw 2 notitle, "
		"$data using 1:2		with circles fill solid lc rgb \"black\" notitle, "
		"$data using 1:2:3		with labels tc rgb \"white\" font \"arial, 8\" notitle\n"
	);

	fprintf(plot_handle->pipe, "pause mouse close\n");

	gpc_close(plot_handle);
}

void plot_tsp_solution_directed(graph* g, double* xstar)
{

	h_GPC_Plot* plot_handle = gpc_init_xy(
		"TSP solution plot",
		"X coords",
		"Y coords",
		GPC_AUTO_SCALE,
		GPC_KEY_DISABLE
	);

	// *************** set gnuplot parameters ***************
	// set window size
	fprintf(plot_handle->pipe, "set term qt size 1000, 800 position 0, 0\n");
	// set relative circle radius
	fprintf(plot_handle->pipe, "set style circle radius character 1.1\n");
	// set offset
	fprintf(plot_handle->pipe, "set offsets graph 0.01, graph 0.01, graph 0.01, graph 0.01\n");
	// set x label
	fprintf(plot_handle->pipe, "set xlabel offset screen 0, 0.07\n");
	// set y label
	fprintf(plot_handle->pipe, "set ylabel offset screen 0.07, 0\n");

	// *************** start nodes block ***************
	fprintf(plot_handle->pipe, "$data_nodes << EOD\n");

	for (int i = 0; i < g->nnodes; i++)
	{
		// print edge nodes
		fprintf(plot_handle->pipe, "%f %f %d\n", g->xcoord[i], g->ycoord[i], i + 1);
	}
	// *************** end data block ***************
	fprintf(plot_handle->pipe, "EOD\n");

	// *************** start vectors data block ***************
	fprintf(plot_handle->pipe, "$data_vectors << EOD\n");

	for (int i = 0; i < g->nnodes; i++)
	{
		for (int j = 0; j < g->nnodes; j++)
		{
			if (i == j) continue;
			if (is_one(xstar[xxpos(i, j, g->nnodes)]))
			{
				// print edge nodes
				fprintf(plot_handle->pipe, "%f %f %f %f\n", g->xcoord[i], g->ycoord[i], g->xcoord[j], g->ycoord[j]);
			}
		}
	}

	// *************** end data block ***************
	fprintf(plot_handle->pipe, "EOD\n");

	// *************** plot data block ***************
	fprintf(plot_handle->pipe,
		"plot "
		"$data_vectors	using 1:2:($3-$1):($4-$2)	with vectors lc rgb \"red\" lw 2 notitle, "
		"$data_nodes	using 1:2					with circles fill solid lc rgb \"black\" notitle, "
		"$data_nodes	using 1:2:3					with labels tc rgb \"white\" font \"arial, 8\" notitle\n"
	);

	fprintf(plot_handle->pipe, "pause mouse close\n");

	gpc_close(plot_handle);
}
