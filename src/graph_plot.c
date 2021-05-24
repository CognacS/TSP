#include "../graph_plot.h"

h_GPC_Plot* setup_tsp_gnuplot()
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

	return plot_handle;
}

/* **************************************************************************************************
*						XSTAR PLOTS
************************************************************************************************** */
void plot_tsp_xstar_undirected(graph* g, double* xstar)
{
	// *************** setup gnuplot ***************
	h_GPC_Plot* plot_handle = setup_tsp_gnuplot();

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

void plot_tsp_xstar_directed(graph* g, double* xstar)
{
	// *************** setup gnuplot ***************
	h_GPC_Plot* plot_handle = setup_tsp_gnuplot();

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

/* **************************************************************************************************
*						SUCC PLOTS
************************************************************************************************** */
void plot_tsp_succ_undirected(graph* g, int* succ)
{
	// *************** setup gnuplot ***************
	h_GPC_Plot* plot_handle = setup_tsp_gnuplot();

	// *************** start data block ***************
	fprintf(plot_handle->pipe, "$data << EOD\n");

	for (int i = 0; i < g->nnodes; i++)
	{
		if (succ[i] >= 0)
		{
			// print edge nodes
			fprintf(plot_handle->pipe, "%f %f %d\n", g->xcoord[i], g->ycoord[i], i + 1);
			fprintf(plot_handle->pipe, "%f %f %d\n", g->xcoord[succ[i]], g->ycoord[succ[i]], succ[i] + 1);
			// print edge separator
			fprintf(plot_handle->pipe, "\n");
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


void plot_tsp_succ_directed(graph* g, int* succ)
{
	// *************** setup gnuplot ***************
	h_GPC_Plot* plot_handle = setup_tsp_gnuplot();

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
		// print edge nodes
		fprintf(plot_handle->pipe, "%f %f %f %f\n", g->xcoord[i], g->ycoord[i], g->xcoord[succ[i]], g->ycoord[succ[i]]);
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

/* **************************************************************************************************
*						SOLUTION PLOTS
************************************************************************************************** */

void plot_tsp_solution_undirected(graph* g, Solution* sol)
{
	if (sol->format & SOLFORMAT_SUCC) plot_tsp_succ_undirected(g, sol->succ);
	else if (sol->format & SOLFORMAT_XSTAR) plot_tsp_xstar_undirected(g, sol->xstar);
}

void plot_tsp_solution_directed(graph* g, Solution* sol)
{
	if (sol->format & SOLFORMAT_SUCC) plot_tsp_succ_directed(g, sol->succ);
	else if (sol->format & SOLFORMAT_XSTAR) plot_tsp_xstar_directed(g, sol->xstar);
}


void plot_tsp_hardfixing_undirected(graph* g, int* succ, char* fixed, int* bridges)
{
	// *************** setup gnuplot ***************
	h_GPC_Plot* plot_handle = setup_tsp_gnuplot();

	// *************** start data block ***************
	fprintf(plot_handle->pipe, "$data_edges << EOD\n");

	for (int i = 0; i < g->nnodes; i++)
	{
		if (fixed[i])
		{
			// print edge nodes
			fprintf(plot_handle->pipe, "%f %f %d\n", g->xcoord[i], g->ycoord[i], i + 1);
			fprintf(plot_handle->pipe, "%f %f %d\n", g->xcoord[succ[i]], g->ycoord[succ[i]], succ[i] + 1);
			// print edge separator
			fprintf(plot_handle->pipe, "\n");
		}
	}

	// *************** end data block ***************
	fprintf(plot_handle->pipe, "EOD\n");

	// *************** start data block ***************
	fprintf(plot_handle->pipe, "$data_bridges << EOD\n");

	// ***************
	for (int i = 0; i < g->nnodes; i++)
	{
		if (bridges[i] >= 0)
		{
			// print bridge nodes
			fprintf(plot_handle->pipe, "%f %f %d\n", g->xcoord[i], g->ycoord[i], i + 1);
			fprintf(plot_handle->pipe, "%f %f %d\n", g->xcoord[bridges[i]], g->ycoord[bridges[i]], bridges[i] + 1);
			// print bridge separator
			fprintf(plot_handle->pipe, "\n");
		}
	}
	// *************** end data block ***************
	fprintf(plot_handle->pipe, "EOD\n");

	// *************** plot data block ***************
	fprintf(plot_handle->pipe,
		"plot "
		"$data_bridges using 1:2	with lines lc rgb \"red\" lw 2 dt 2 notitle, "
		"$data_edges using 1:2		with lines lc rgb \"red\" lw 2 notitle, "
		"$data_edges using 1:2		with circles fill solid lc rgb \"black\" notitle, "
		"$data_edges using 1:2:3	with labels tc rgb \"white\" font \"arial, 8\" notitle\n"
	);

	fprintf(plot_handle->pipe, "pause mouse close\n");

	gpc_close(plot_handle);
}

void plot_heuristic_perflog(LinkedList* ll)
{
	h_GPC_Plot* plot_handle = setup_tsp_gnuplot();
	// set axis labels
	fprintf(plot_handle->pipe, "set xlabel 'Time (s)'\n");
	fprintf(plot_handle->pipe, "set ylabel 'Objective value'\n");

	// *************** start data block ***************
	fprintf(plot_handle->pipe, "$data << EOD\n");

	for (Cell* cell = ll->start; cell != NULL; cell = cell->next)
	{
		// print point
		fprintf(plot_handle->pipe, "%f %f\n", cell->x, cell->y);
	}

	// *************** end data block ***************
	fprintf(plot_handle->pipe, "EOD\n");

	// *************** plot data block ***************
	fprintf(plot_handle->pipe,
		"plot "
		"$data using 1:2 with lp pt 7 ps 1.5 title 'Heuristic objective', "
		"$data using 1:2 with i lt 0 lw 0.5 lc rgb \"#ff0000\"\n"
	);

	fprintf(plot_handle->pipe, "pause mouse close\n");

	gpc_close(plot_handle);

}