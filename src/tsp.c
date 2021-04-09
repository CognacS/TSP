#include "../tsp.h"

/* ***************************************************************************************************
*						TSP OPTIMIZATION PROCEDURE
*************************************************************************************************** */
// optimization functions
int TSPopt(instance* inst)
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
	CPXsetdblparam(env, CPX_PARAM_TILIM, p->timelimit);
	CPXsetintparam(env, CPX_PARAM_NODELIM, p->max_nodes);
	CPXsetdblparam(env, CPX_PARAM_EPGAP, p->cutoff);
	CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);

	// apply transformation to coordinates if required
	coord_transform(g);

	// ************************ OPTIMIZATION ************************
	switch (model_tsptype(p->model_type))
	{
	case MODEL_TSP_ASYMM:
		solve_asymmetric_tsp(inst, env, lp);
		break;
	case MODEL_TSP_SYMM:
		solve_symmetric_tsp(inst, env, lp);
		break;
	default:
		print_error("TSP variant", ERR_MODEL_NOT_IMPL);
	}

	if (VERBOSITY >= LOGLVL_INFO) CPXwriteprob(env, lp, "model.lp", NULL);

	// ************************ USE SOLUTION ************************
	int nnodes = inst->inst_graph.nnodes;
	int ncols = CPXgetnumcols(env, lp);
	global_data* gd = &inst->inst_global_data;

	// get the optimal solution
	if (gd->xstar != NULL) free(gd->xstar);
	calloc_s(gd->xstar, ncols, double);
	if (error = CPXgetx(env, lp, gd->xstar, 0, ncols - 1))
	{
		printf("CPX error %d\n", error);
		print_error("CPXgetx()", ERR_CPLEX);
	}

	// get the obj value
	if (error = CPXgetobjval(env, lp, &gd->zbest))
	{
		printf("CPX error %d\n", error);
		print_error("CPXgetobjval()", ERR_CPLEX);
	}
	printf("Final objective value = %f\n", gd->zbest);

	// get the lower bound
	if (error = CPXgetbestobjval(env, lp, &gd->lbbest))
	{
		printf("CPX error %d\n", error);
		print_error("CPXgetobjval()", ERR_CPLEX);
	}


	// ************************ CLEAN-UP ************************

	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return 0; // or an appropriate nonzero error code

}


