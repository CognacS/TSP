#include "../include/heuristics.h"

/* **************************************************************************************************
*						SOLUTION USING HEURISTICS (MASTER FUNCTION)
************************************************************************************************** */
void solve_heuristically(OptData* optdata)
{
	// ***************************** SETUP HEURISTIC *****************************
	// build an heuristic datastructure and fill it with parameters
	Heuristic heur;
	decode_heuristic(optdata, &heur);

	// if cplex is required during the heuristic procedure, build the best cplex
	// model available
	callback_instance cb_inst;
	if (heur.requires_cplex)
	{
		// ********************** SETUP BASE MODEL **********************
		// setup cplex for use
		mip_setup_cplex(optdata);

		// build naive model
		build_model_base_undirected(optdata);

		// ********************** BUILD DATASTRUCTURES **********************
		// unpack
		instance* inst = optdata->inst;
		CPXENVptr env = optdata->cpx->env;
		CPXLPptr lp = optdata->cpx->lp;

		// build callback
		int nnodes = inst->inst_graph.nnodes;
		int nedges = nnodes * (nnodes - 1) / 2;
		int* elist;	arr_malloc_s(elist, 2 * nedges, int);
		int arr_idx = 0;
		for (int i = 0; i < nnodes; i++)
		{
			for (int j = i + 1; j < nnodes; j++)
			{
				elist[arr_idx++] = i;
				elist[arr_idx++] = j;
			}
		}
		// make callback instance
		cb_inst.inst = inst;
		cb_inst.args = elist;
		cb_inst.sep_procedure = CC_add_sec_on_subtours;
		cb_inst.rej_procedure = add_sec_on_subtours;
		cb_inst.prob_decay = CUT_HI_COEFF_DECAY;
		cb_inst.prob_function = exp_decay_prob;

		// ********************** SETUP CALLBACK FUNCTION **********************
		CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION;
		int error;
		// set the callback function
		if (error = CPXcallbacksetfunc(env, lp, contextid, sec_callback, &cb_inst))
			print_error_ext(ERR_CPLEX, "CPXcallbacksetfunc() TSP solver, CPX error: %d", error);
	}

	// ********************** SETUP SOLUTION DATASTRUCTURES **********************
	// setup starting solution
	Solution final_sol;
	// if need to log progress during the heuristic, build a performance log
	if (VERBOSITY >= LOGLVL_PLOTSOL)
	{
		optdata->perflog = newLinkedList();
	}

	// ********************** HEURISTIC OPTIMIZATION PROCEDURE ***********************
	// invoke backbone method and optimize
	heur.backbone_heur(optdata, &final_sol, &heur, heur.backbone_data);
	// *******************************************************************************


	// ********************** WORK WITH FOUND SOLUTION **********************
	// pass solution to instance
	convert_solution(&final_sol, SOLFORMAT_XSTAR);
	global_data* gd = &optdata->inst->inst_global_data;
	gd->xstar = final_sol.xstar;
	gd->zbest = final_sol.cost;

	// if need to plot the progress during heuristic
	if (VERBOSITY >= LOGLVL_PLOTSOL)
	{
		if (!LL_is_empty(optdata->perflog)) plot_heuristic_perflog(optdata->perflog);
		LL_free(optdata->perflog);
		optdata->perflog = NULL;
	}

	// **************************** CLEANUP ****************************
	// deallocate structures
	free_s(heur.backbone_data);
	free_s(heur.construct_data);
	free_s(heur.refine_data);

	// cplex will be deallocated outside together with all other methods
}

/* **************************************************************************************************
*					DECODE HEURISTIC STRING GIVEN IN INPUT INSIDE INSTANCE
************************************************************************************************** */
void decode_heuristic(OptData* optdata, Heuristic* heur)
{
	// unpack
	instance* inst = optdata->inst;
	CPXENVptr env = optdata->cpx->env;
	CPXLPptr lp = optdata->cpx->lp;
	char* code = inst->inst_params.heuristic_code;

	// local parameters
	char method = 0;
	char* params;
	char* token;
	void* data;

	heur->requires_cplex = 0;

	// ******************************* DECODE BACKBONE *******************************
	method = *code;		// method is first character of code
	code += 2;
	params = strtok(code, ")");	// params are comprised in '('...')' and separated by ','
	code += strlen(params);		// bring code to end of section

	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "Backbone method: %c", method);
	switch (method)
	{
	case 'i':
		heur->backbone_heur = backbone_iter_local_search;
		malloc_s(data, IterativeData);
		// construct time
		token = strtok(params, ",");
		((IterativeData*)data)->construct_timelimit = atof(token);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "construct_time: \"%s\"", token);
		// refine time
		token = strtok(NULL, ",");
		((IterativeData*)data)->refine_timelimit = atof(token);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "refine_time: \"%s\"", token);
		// assign data
		heur->backbone_data = data;
		heur->backbone_sol_format = SOLFORMAT_XSTAR;
		break;
	default:
		print_error_ext(ERR_HEUR_DECODE, "unknown backbone method: %c", method);
	}

	code++;
	// ******************************* DECODE CONSTRUCT ******************************
	method = *code;		// method is first character of code
	code += 2;
	params = strtok(code, ")");	// params are comprised in '('...')' and separated by ','
	code += strlen(params);		// bring code to end of section

	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "Construct method: %c", method);
	switch (method)
	{
	case 's':	// construct through solver
		heur->construct_heur = construct_tspsolver;
		heur->requires_cplex = 1;
		heur->construct_data = NULL;
		heur->construct_sol_format = SOLFORMAT_XSTAR;
		break;
	case 'g':	// construct through greedy algorithm
		heur->construct_heur = construct_greedy;
		malloc_s(data, GreedyAlgoData);
		// deviation probability
		token = strtok(params, ",");
		((GreedyAlgoData*)data)->grasp.p_dev = atof(token);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "p_dev: \"%s\"", token);
		// candidate pool size
		token = strtok(NULL, ",");
		((GreedyAlgoData*)data)->grasp.candpool_size = atoi(token);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "candpool_size: \"%s\"", token);
		// assign data
		heur->construct_data = data;
		heur->construct_sol_format = SOLFORMAT_BOTH;
		break;
	case 'e':	// construct through extra mileage algorithm
		heur->construct_heur = construct_extramileage;
		malloc_s(data, ExtraMileageData);
		// deviation probability
		token = strtok(params, ",");
		((ExtraMileageData*)data)->grasp.p_dev = atof(token);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "p_dev: \"%s\"", token);
		// candidate pool size
		token = strtok(NULL, ",");
		((ExtraMileageData*)data)->grasp.candpool_size = atoi(token);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "candpool_size: \"%s\"", token);
		// variant for starting edges choice
		token = strtok(NULL, ",");
		((ExtraMileageData*)data)->variant_startingchoice = *token;
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "variant_startingchoice: \"%s\"", token);
		// assign data
		heur->construct_data = data;
		heur->construct_sol_format = SOLFORMAT_BOTH;
		break;
	default:
		print_error_ext(ERR_HEUR_DECODE, "unknown constructive heuristic: %c", *code);
	}
	code++;
	// ******************************** DECODE REFINE ********************************
	method = *code;		// method is first character of code
	code += 2;
	params = strtok(code, ")");	// params are comprised in '(...')' and separated by ','
	code += strlen(params);		// bring code to end of section

	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "Refine method: %c", method);
	switch (method)
	{
	case 'h':	// refine through hardfixing scheme
		heur->refine_heur = refine_hardfixing;
		heur->requires_cplex = 1;
		malloc_s(data, HardfixingData);
		// variant
		token = strtok(params, ",");
		((HardfixingData*)data)->variant = *token;
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "variant: \"%s\"", token);
		// values
		for (int i = 0; i < HEUR_HARDFIX_VALUES; i++)
		{
			token = strtok(NULL, ",");
			((HardfixingData*)data)->values[i] = atof(token);
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "values[%d]: \"%s\"", i, token);
		}
		// assign data
		((HardfixingData*)data)->progress = 0;
		heur->refine_data = data;
		heur->refine_sol_format = SOLFORMAT_XSTAR;
		break;
	case 'l':	// refine through local branching scheme
		heur->refine_heur = construct_greedy;
		heur->requires_cplex = 1;
		malloc_s(data, LocalbranchingData);
		// variant
		token = strtok(params, ",");
		((LocalbranchingData*)data)->variant = *token;
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "variant: \"%s\"", token);
		// values
		for (int i = 0; i < HEUR_HARDFIX_VALUES; i++)
		{
			token = strtok(NULL, ",");
			((LocalbranchingData*)data)->values[i] = atoi(token);
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "values[%d]: \"%s\"", i, token);
		}
		// assign data
		((LocalbranchingData*)data)->progress = 0;
		heur->refine_data = data;
		heur->refine_sol_format = SOLFORMAT_XSTAR;
		break;
	default:
		print_error_ext(ERR_HEUR_DECODE, "unknown refining heuristic: %c", *code);
	}

}

/* **************************************************************************************************
*					RESOLVE CONVERSIONS INSIDE THE BACKBONE HEURISTIC METHOD
************************************************************************************************** */
void resolve_solformats(Heuristic* heur)
{
	if (heur->construct_sol_format == SOLFORMAT_BOTH)
	{
		if (heur->refine_sol_format == SOLFORMAT_BOTH)
		{
			// if both are free, choose the preferred
			heur->construct_sol_format = heur->refine_sol_format = HEUR_PREFERRED_FORMAT;
		}
		else
		{
			// if constr is free but refine isn't, use the refine format
			heur->construct_sol_format = heur->refine_sol_format;
		}
	}

	if (heur->refine_sol_format == SOLFORMAT_BOTH &&
		heur->construct_sol_format != SOLFORMAT_BOTH)
	{
		// if refine is free but constr isn't, use the constr format
		heur->refine_sol_format = heur->construct_sol_format;
	}
}

/* **************************************************************************************************
*					BACKBONE METHOD: ITERATIVE LOCAL SEARCH
*			- construct a starting solution
*			- refine iteratively the solution until time expires
************************************************************************************************** */
void backbone_iter_local_search(OptData* optdata, Solution* sol, Heuristic* heur, void* data)
{
	// ******************************** SETUP ********************************
	// unpack
	instance* inst = optdata->inst;
	IterativeData* bb_data = (IterativeData*)data;

	// resolve solution format problems (xstar, succ)
	resolve_solformats(heur);
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Solution formats: bb=%d ct=%d rf=%d",
		heur->backbone_sol_format, heur->construct_sol_format, heur->refine_sol_format);

	// startup solution and format it for construction
	empty_solution(sol, heur->construct_sol_format, inst->inst_graph.nnodes);

	// ************************* CONSTRUCT STARTING SOLUTION *************************
	double timelimit, time_passed = 0, time_phase, t_start;
	// compute timelimit for constructive heuristic
	timelimit = min(residual_time(inst), bb_data->construct_timelimit);
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Construction should run for %f sec.s", timelimit);

	// construct starting heuristic solution
	t_start = second();
	// *********** CONSTRUCT ***********
	heur->construct_heur(optdata, sol, heur->construct_data, timelimit);
	// *********************************
	time_phase = second() - t_start;
	time_passed += time_phase;

	// log
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Construction done in %f sec.s", time_phase);
	log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE]: Obj value after construction = %f", sol->cost);
	if (VERBOSITY >= LOGLVL_PLOTSOL) LL_add_value(optdata->perflog, time_passed, sol->cost);
	if (VERBOSITY >= LOGLVL_DEBUGPLOT) plot_tsp_solution_undirected(&inst->inst_graph, sol);


	// ****************************** START ITERATIONS *******************************
	// convert solution to the refinement format
	convert_solution(sol, heur->refine_sol_format);

	while (!time_limit_expired(inst))
	{
		// ***************************** REFINE SOLUTION *****************************
		// compute timelimit for refinement heuristic
		timelimit = min(residual_time(inst), bb_data->refine_timelimit);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] 1 iteration of refinement should run for %f sec.s", timelimit);

		// refine heuristic solution
		t_start = second();
		// *********** REFINE ***********
		heur->refine_heur(optdata, sol, heur->refine_data, timelimit);
		// ******************************
		time_phase = second() - t_start;
		time_passed += time_phase;

		// log
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] 1 iteration of refinement done in %f sec.s", time_phase);
		log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE]: Obj value after refinement = %f", sol->cost);
		if (VERBOSITY >= LOGLVL_PLOTSOL) LL_add_value(optdata->perflog, time_passed, sol->cost);
		if (VERBOSITY >= LOGLVL_DEBUGPLOT) plot_tsp_solution_undirected(&inst->inst_graph, sol);
	}
}

/* **************************************************************************************************
*					CONSTRUCTIVE HEURISTIC: CPLEX TSP SOLVER
*			- run the best cplex tsp solver for a fixed amount of time
*			- make cplex emphasize finding a good heuristic solution
************************************************************************************************** */
void construct_tspsolver(OptData* optdata, Solution* sol, void* data, double timelim)
{
	// unpack
	instance* inst = optdata->inst;
	CPXENVptr env = optdata->cpx->env;
	CPXLPptr lp = optdata->cpx->lp;

	// ********************** OPTIMIZE **********************
	// solve the problem with the callback
	mip_timelimit(optdata, timelim);
	CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_HEURISTIC);
	int error;
	if (error = CPXmipopt(env, lp))
		print_error_ext(ERR_CPLEX, "CPXmipopt() TSP solver, CPX error: %d", error);

	// ********************** EXTRACT **********************
	// get the optimal solution
	mip_extract_sol_obj(optdata, sol, "TSP solver");

}

void construct_greedy(OptData* optdata, Solution* sol, void* data, double timelim)
{

}
void construct_extramileage(OptData* optdata, Solution* sol, void* data, double timelim)
{

}

/* **************************************************************************************************
*					REFINING HEURISTIC: HARD FIXING
*			- fix a certain % of variables (edges) to 1 (open paths are formed)
*			- fix the bridges linking the extremes of paths to 0
*			- re-optimize the simplified problem using cplex
*			- unfix all variables
*			VARIANTS: scheduled %, fixed %, random %
************************************************************************************************** */
void refine_hardfixing(OptData* optdata, Solution* sol, void* data, double timelim)
{
	// ************************************ SETUP ************************************
	// unpack optdata
	instance* inst = optdata->inst;
	CPXENVptr env = optdata->cpx->env;
	CPXLPptr lp = optdata->cpx->lp;

	// unpack rfndata
	HardfixingData* rfn_data = (HardfixingData*)data;
	char variant = rfn_data->variant;
	double* values = rfn_data->values;
	int progress = rfn_data->progress;

	// extract data
	int nnodes = sol->nnodes;

	// extract current solution
	double* xstar = sol->xstar;
	int ncols = inst->inst_model.ncols;

	// ******************************** % COMPUTATION ********************************
	// compute % of variables to fix
	double fix_ratio;
	switch (variant)
	{
	case 's':
		fix_ratio = values[progress];
		break;
	case 'f':
		fix_ratio = values[0];
		break;
	case 'u':
		fix_ratio = ((double)rand() / RAND_MAX * (values[1] - values[0])) + values[0];
		break;
	default:
		print_error_ext(ERR_HEUR_CONSTR_PARAM, "variant, got  %c, should be s/f/u", variant);
	}
	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Will fix %f%% of variables", fix_ratio);
	// compute number of variables to fix
	int n_fix = (int)(fix_ratio * nnodes);

	// ************************ SELECTION OF VARIABLES TO FIX ************************
	// allocate arrays
	int* succ;		arr_malloc_s(succ, nnodes, int);	// succ node in tour
	int* nodes;		arr_malloc_s(nodes, nnodes, int);	// nodes in tour
	for (int i = 0; i < nnodes; i++) nodes[i] = i;
	int* bridges = nodes;								// bridges between extremes
	char* fixed;	calloc_s(fixed, nnodes, char);		// fixed or not nodes
	
	// avoid computations if no variables are to be fixed
	if (n_fix > 0)
	{

		// *************** identify tour of FEASABLE solution **************
		if (!xstar2succ(xstar, succ, nnodes)) print_error(ERR_HEUR_INFEASABLE_SOL, NULL);


		// *************** shuffle n_fix variables and fix them to 1 **************
		log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Fixing %d variables", n_fix);
		int idx, var_idx, start_node;
		char lu = 'L';
		double bd = 1.0;
		for (int i = 0; i < n_fix; i++)
		{
			// select a random edge/variable from the remaining ones
			idx = (int)((double)rand() / RAND_MAX * (nnodes - i)) + i;
			start_node = nodes[idx];
			// fix variable to 1
			fixed[start_node] = 1;
			var_idx = xpos(start_node, succ[start_node], nnodes);
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Fixing edge %d->%d, idx: %d", start_node, succ[start_node], var_idx);
			CPXchgbds(env, lp, 1, &var_idx, &lu, &bd);
			// stack the variable on top
			nodes[idx] = nodes[i];
			nodes[i] = start_node;
		}

		// *************** find bridges and fix them to 0 **************
		int path_len = 0;
		int curr = succ[nodes[n_fix]]; // start tour from edge next to a non-fixed edge
		start_node = -1;
		lu = 'U';
		bd = 0.0;
		for (int i = 0; i < nnodes; i++) bridges[i] = -1; // reset bridge from node i

		for (int i = 0; i < nnodes; i++)
		{
			// if the edge is fixed -> increase length of fixed path
			if (fixed[curr]) path_len++;

			// if not started and the next edge is fixed -> start of path
			if (start_node == -1 && fixed[curr])
			{
				start_node = curr; // start path from the current node
			}
			// else if started and the next edge is not fixed -> end of path
			else if (start_node >= 0 && !fixed[curr])
			{
				// fix to 0 the bridge only if its path is bigger than just an edge
				if (path_len > 1)
				{
					bridges[start_node] = curr; // link bridge from start to curr
					// fix variable to 0
					var_idx = xpos(start_node, curr, nnodes);
					CPXchgbds(env, lp, 1, &var_idx, &lu, &bd);
					log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Fixing bridge %d->%d, idx: %d", start_node, curr, var_idx);
				}
				else log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Skipped bridge %d->%d because it is on a path of length 1", start_node, curr);

				path_len = 0;		// reset path
				start_node = -1;	// reset start_node
			}
			curr = succ[curr];
		}

		// plot the fixed edges and bridges
		if (VERBOSITY >= LOGLVL_DEBUGPLOT_P) plot_tsp_hardfixing_undirected(&inst->inst_graph, succ, fixed, bridges);
	}

	// ******************************** OPTIMIZATION *********************************
	// optimize with cplex
	int error;
	mip_timelimit(optdata, timelim);
	mip_warmstart(optdata, xstar);
	CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_BALANCED);
	if (error = CPXmipopt(env, lp))
		print_error_ext(ERR_CPLEX, "CPXmipopt() Hardfixing, CPX error: %d", error);

	// ********************** EXTRACT SOLUTION AND UPDATE METHOD *********************
	// extract best solution (expected to be already allocated in the construction phase!)
	double obj_before = sol->cost;
	mip_extract_sol_obj(optdata, sol, "Hardfixing");
	double obj_after = sol->cost;
	// if using schedule and obj value is not good enough
	if (variant == 's' && obj_after >= obj_before * IMPROVEMENT_RATIO)
	{
		// increase progress in schedule (capped)
		rfn_data->progress = min(progress + 1, HEUR_HARDFIX_VALUES-1);
		log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Hardfix schedule progress to %d", rfn_data->progress);
	}
	else
	{
		log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Solution improved by %f%%", (obj_before- obj_after)/ obj_before);

	}
	
	// ******************************* UNFIX VARIABLES *******************************
	if (n_fix > 0)
	{
		int var_idx;
		char lu;
		double bd;
		// unfix variables
		for (int i = 0; i < nnodes; i++)
		{
			// unfix tour edges
			if (fixed[i])
			{
				var_idx = xpos(i, succ[i], nnodes);
				log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Unfixing edge %d", var_idx);
				lu = 'L'; bd = 0.0;
				CPXchgbds(env, lp, 1, &var_idx, &lu, &bd);
			}
			// unfix bridges
			if (bridges[i] >= 0)
			{
				var_idx = xpos(i, bridges[i], nnodes);
				log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Unfixing bridge %d", var_idx);
				lu = 'U'; bd = 1.0;
				CPXchgbds(env, lp, 1, &var_idx, &lu, &bd);
			}
		}
	}

	// *********************************** CLEANUP ***********************************
	// cleanup
	free(succ);
	free(nodes);
	free(fixed);
	
}

void refine_localbranching(OptData* optdata, Solution* sol, void* data, double timelim)
{

}