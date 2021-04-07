#include "../tsp.h"

/* ***************************************************************************************************
*						TSP OPTIMIZATION PROCEDURE
*************************************************************************************************** */
// optimization functions
int TSPopt(instance* inst, double** xstar_ptr)
{

	// *************************** SETUP ***************************
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	// configure environment parameters
	params* p = &(inst->inst_params);
	graph* g = &(inst->inst_graph);
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, p->randomseed);
	CPXsetintparam(env, CPX_PARAM_TILIM, p->timelimit);
	CPXsetintparam(env, CPX_PARAM_NODELIM, p->max_nodes);
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 0.00001);
	CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);

	// apply transformation to coordinates if required
	coord_transform(g);

	// ************************ OPTIMIZATION ************************
	switch (model_tsptype(p->model_type))
	{
	case TSP_ASYMM:
		solve_asymmetric_tsp(g, p->model_type, env, lp);
		break;
	case TSP_SYMM:
		solve_symmetric_tsp(g, p->model_type, env, lp);
		break;
	default:
		print_error("TSP variant", ERR_MODEL_NOT_IMPL);
	}

	if (VERBOSITY >= LOGLVL_INFO) CPXwriteprob(env, lp, "model.lp", NULL);

	// ************************ USE SOLUTION ************************
	int nnodes = inst->inst_graph.nnodes;
	int ncols = CPXgetnumcols(env, lp);

	// get the optimal solution
	double* xstar = 0;
	if (!(xstar = (double*)calloc(ncols, sizeof(double)))) print_error("xstar", ERR_NO_MEM_FOR_ALLOC);
	*xstar_ptr = xstar;
	if (error = CPXgetx(env, lp, xstar, 0, ncols - 1))
	{
		printf("CPX error %d\n", error);
		print_error("CPXgetx()", ERR_CPLEX);
	}

	double objvalue = 0;
	CPXgetobjval(env, lp, &objvalue);
	if (error = CPXgetx(env, lp, xstar, 0, ncols - 1))
	{
		printf("CPX error %d\n", error);
		print_error("CPXgetobjval()", ERR_CPLEX);
	}
	printf("Final objective value = %f\n", objvalue);

	// iterate through all variables
	

	// ************************ CLEAN-UP ************************

	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return 0; // or an appropriate nonzero error code

}

void fill_inst_default(instance* inst)
{
	// define the parameters for the instance
	graph* g = &inst->inst_graph;
	params* p = &inst->inst_params;

	// ************ DEFAULT PARAMETERS DEFINITION ************
	g->integer_costs = DEF_INTEGER_COSTS;
	g->xcoord = NULL;
	g->ycoord = NULL;
	g->tr_xcoord = NULL;
	g->tr_ycoord = NULL;
	p->model_type = DEF_MODEL_TYPE;
	strcpy(p->input_file, DEF_INPUT_FILE);
	strcpy(p->batch_file, DEF_BATCH_FILE);
	p->timelimit = DEF_TIMELIMIT;
	p->randomseed = DEF_RANDOMSEED;
	p->max_nodes = DEF_MAX_NODES;
	// *******************************************************
	
}


/* ***************************************************************************************************
*						DATASTRUCTURES TOSTRING FUNCTIONS
*************************************************************************************************** */


void tostring_graph(char* buffer, graph* inst_graph)
{
	sprintf(buffer, "Graph of %d nodes, distances are int=%d, coords:\n", inst_graph->nnodes, inst_graph->integer_costs);
	char line[100];
	for (int i = 0; i < inst_graph->nnodes; i++)
	{
		sprintf(line, "(%lf, %lf), ", inst_graph->xcoord[i], inst_graph->ycoord[i]);
		strcat(buffer, line);
	}
	strcat(buffer, "\n");
}

void tostring_params(char* buffer, params* inst_p)
{
	sprintf(buffer, "Parameters:\n"
		"\tmodel_type = %d\n"
		"\ttimelimit  = %f\n"
		"\tinput_file = %s\n"
		"\tbatch_file = %s\n"
		"\trandomseed = %d\n"
		"\tcutoff = %f\n"
		"\tmax_nodes = %d\n",
		inst_p->model_type,
		inst_p->timelimit,
		inst_p->input_file,
		inst_p->batch_file,
		inst_p->randomseed,
		inst_p->cutoff,
		inst_p->max_nodes
	);
}


void tostring_global_data(char* buffer, global_data* inst_global)
{
	sprintf(buffer, "Global Data:\n"
		"\tz_best = %f\n",
		inst_global->z_best
	);
}
void tostring_model(char* buffer, model* inst_model)
{
	sprintf(buffer, "Model:\n"
		"\tplc_holder = %d\n",
		inst_model->plc_holder
	);
}

void tostring_instance(char* buffer, instance* inst)
{
	char mini_buffers[4][2000];
	tostring_graph(mini_buffers[0], &inst->inst_graph);
	tostring_params(mini_buffers[1], &inst->inst_params);
	tostring_global_data(mini_buffers[2], &inst->inst_global_data);
	tostring_model(mini_buffers[3], &inst->inst_model);

	sprintf(buffer, "********** INSTANCE REPORT **********\n%s\n%s\n%s\n%s\n*************************************\n",
		mini_buffers[0], mini_buffers[1], mini_buffers[2], mini_buffers[3]);

}

/* ***************************************************************************************************
*						DATASTRUCTURES DESTROYERS
*************************************************************************************************** */

void free_graph(graph* inst_graph)
{
	free(inst_graph->xcoord);
	free(inst_graph->ycoord);
	free(inst_graph->tr_xcoord);
	free(inst_graph->tr_ycoord);
}
void free_instance(instance* inst)
{
	free_graph(&inst->inst_graph);
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