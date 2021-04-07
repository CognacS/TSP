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
	CPXsetintparam(env, CPX_PARAM_TILIM, (int)p->timelimit);
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