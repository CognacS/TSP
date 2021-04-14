#include "../tsp.h"

/* ***************************************************************************************************
*						TSP OPTIMIZATION PROCEDURE
*************************************************************************************************** */
// optimization functions
int TSPopt(instance* inst, char getsol)
{
	double starting_time = second();

	// *************************** SETUP ***************************
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	// configure environment parameters
	params* p = &(inst->inst_params);
	graph* g = &(inst->inst_graph);
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, p->randomseed);
	CPXsetintparam(env, CPX_PARAM_NODELIM, p->max_nodes);
	CPXsetdblparam(env, CPX_PARAM_EPGAP, p->cutoff);
	CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0);

	// configure starting parameters
	global_data* gd = &inst->inst_global_data;
	gd->tstart = starting_time;
	gd->texec = 0;

	// setup time limit
	mip_timelimit(env, p->timelimit, inst);

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
		print_error(ERR_MODEL_NOT_IMPL, "TSP variant");
	}

	if (VERBOSITY >= LOGLVL_INFO) CPXwriteprob(env, lp, "model.lp", NULL);

	// get run execution time
	gd->texec = second() - gd->tstart;

	// ************************* CHECK RUN *************************
	// setup code to return after optimization procedure
	opt_result output_code = OPT_OK;
	if (time_limit_expired(inst))
	{
		output_code = OPT_TL_EXPIRED;
	}
	// ************************ USE SOLUTION ************************
	else
	{
		if (getsol)
		{
			int nnodes = inst->inst_graph.nnodes;
			int ncols = CPXgetnumcols(env, lp);

			// get the optimal solution
			if (gd->xstar != NULL) free(gd->xstar);
			calloc_s(gd->xstar, ncols, double);
			if (error = CPXgetx(env, lp, gd->xstar, 0, ncols - 1))
				print_error_ext(ERR_CPLEX, "CPXgetx() Finalized, CPX error: %d", error);

			// get the obj value
			if (error = CPXgetobjval(env, lp, &gd->zbest))
				print_error_ext(ERR_CPLEX, "CPXgetobjval(), CPX error: %d", error);

			log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE]: Final objective value = %f\n", gd->zbest);

			// get the lower bound
			if (error = CPXgetbestobjval(env, lp, &gd->lbbest))
				print_error_ext(ERR_CPLEX, "CPXgetbestobjval(), CPX error: %d", error);
		}
	}
	// ************************ CLEAN-UP ************************

	// free and close cplex model 
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return output_code; // or an appropriate nonzero error code

}


