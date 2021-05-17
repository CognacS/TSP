#include "../include/tsp.h"


/* ***************************************************************************************************
*						TSP OPTIMIZATION PROCEDURE
*************************************************************************************************** */
// optimization functions
int TSPopt(instance* inst)
{
	double starting_time = second();

	// *************************** SETUP ***************************
	// initialize optimization procedure wrapping of data
	OptData optdata;
	optdata.inst = inst;
	optdata.cpx = NULL;

	// configure environment parameters
	params* p = &(inst->inst_params);
	graph* g = &(inst->inst_graph);
	model* m = &(inst->inst_model);
	global_data* gd = &inst->inst_global_data;

	// configure starting parameters
	gd->tstart = starting_time;
	gd->texec = 0;
	gd->perf_measure = 0;

	// apply transformation to coordinates if required
	coord_transform(g);

	// setup random seed
	srand(DEFAULT_SEED);

	// setup cplex support
	if (model_tsptype(p->model_type) == MODEL_TSP_ASYMM ||
		model_tsptype(p->model_type) == MODEL_TSP_SYMM)
	{
		mip_setup_cplex(&optdata);
	}

	// ************************ OPTIMIZATION ************************
	switch (model_tsptype(p->model_type))
	{
	case MODEL_TSP_ASYMM:
		solve_asymmetric_tsp(&optdata);
		break;
	case MODEL_TSP_SYMM:
		solve_symmetric_tsp(&optdata);
		break;
	case MODEL_HEURISTICS:
		solve_heuristically(&optdata);
		break;
	default:
		print_error(ERR_MODEL_NOT_IMPL, "TSP variant");
	}

	// print LP model if required
	if (using_cplex(&optdata) && VERBOSITY >= LOGLVL_DEBUG)
		CPXwriteprob(optdata.cpx->env, optdata.cpx->lp, "model.lp", NULL);

	// get run execution time
	gd->texec = second() - gd->tstart;

	// ************************* CHECK RUN AND FINALIZE *************************
	// setup code to return after optimization procedure
	opt_result output_code = OPT_OK;
	switch (model_tsptype(p->model_type))
	{
	case MODEL_TSP_ASYMM:
	case MODEL_TSP_SYMM:
		output_code = OPT_OK;
		// IF TIMELIMIT EXPIRED
		if (time_limit_expired(inst))
		{
			output_code = OPT_TL_EXPIRED;
			gd->perf_measure *= OPT_TL_PENALIZATION_MULT;
			break;
		}
		// IF ALL WENT WELL
		// set performance measure as execution time
		gd->perf_measure = gd->texec;
		// extract solution
		if (gd->xstar != NULL) free_s(gd->xstar);
		arr_malloc_s(gd->xstar, m->ncols, double);
		mip_extract_sol_obj_lb(&optdata, "Finalized");
		break;
	case MODEL_HEURISTICS:
		output_code = OPT_HEUR_OK;
		// set performance measure as 
		gd->perf_measure = gd->texec;
		// no need to extract solution -> it is already provided by the heuristic
		break;
	}
	// ************************ CLEAN-UP ************************
	// free and close cplex model 
	if (optdata.cpx != NULL)
	{
		mip_close_cplex(&optdata);
	}
	
	return output_code; // or an appropriate nonzero error code

}


