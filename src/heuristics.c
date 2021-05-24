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
		log_line(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Cplex needed: setting up"),

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

	// resolve solution format problems (xstar, succ)
	resolve_solformats(&heur);
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Solution formats: bb=%x ct=%x rf=%x",
		heur.backbone_sol_format, heur.construct_sol_format, heur.refine_sol_format);
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Possible formats are: \n"
	"- %x : succ\n"
	"- %x : xstar\n"
	"- %x : chromo\n"
	"and combinations",
		SOLFORMAT_SUCC, SOLFORMAT_XSTAR, SOLFORMAT_CHROMO);

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
	heur->to_timelimit = 1;

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
		heur->to_timelimit = 0;
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
		heur->construct_sol_format = SOLFORMAT_XSTAR | SOLFORMAT_SUCC;
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
		heur->construct_sol_format = SOLFORMAT_SUCC;
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
		// max_rounds
		token = strtok(NULL, ",");
		((HardfixingData*)data)->max_rounds = atoi(token);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "max_rounds: \"%s\"", token);
		// values
		for (int i = 0; i < HEUR_HARDFIX_VALUES; i++)
		{
			token = strtok(NULL, ",");
			((HardfixingData*)data)->values[i] = atof(token);
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "values[%d]: \"%s\"", i, token);
		}
		// assign data
		((HardfixingData*)data)->rounds_nonimprove = 0;
		((HardfixingData*)data)->progress = 0;
		heur->refine_data = data;
		heur->refine_sol_format = SOLFORMAT_XSTAR;
		break;
	case 'l':	// refine through local branching scheme
		heur->refine_heur = refine_localbranching;
		heur->requires_cplex = 1;
		malloc_s(data, LocalbranchingData);
		// max_rounds
		token = strtok(params, ",");
		((LocalbranchingData*)data)->max_rounds = atoi(token);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "max_rounds: \"%s\"", token);
		// values
		for (int i = 0; i < HEUR_LOCALBRANCH_VALUES; i++)
		{
			token = strtok(NULL, ",");
			((LocalbranchingData*)data)->values[i] = atoi(token);
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "values[%d]: \"%s\"", i, token);
		}
		// assign data
		((LocalbranchingData*)data)->rounds_nonimprove = 0;
		((LocalbranchingData*)data)->progress = 0;
		heur->refine_data = data;
		heur->refine_sol_format = SOLFORMAT_XSTAR;
		break;
	case '2':	// refine through 2-OPT moves
		heur->refine_heur = refine_2opt;
		heur->refine_data = NULL;
		heur->refine_sol_format = SOLFORMAT_SUCC;
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
	// emphasis on common format between backbone and refinement
	solformat bb_rf_common_format =
		heur->backbone_sol_format &
		heur->refine_sol_format;
	
	// collapse common format in backbone and refinement
	bb_rf_common_format = solformat_collapse(bb_rf_common_format);

	if (bb_rf_common_format)
	{
		heur->backbone_sol_format = bb_rf_common_format;
		heur->refine_sol_format = bb_rf_common_format;
	}

	// finally find common format between backbone and construction
	solformat cn_bb_common_format =
		heur->construct_sol_format &
		heur->backbone_sol_format;

	// collapse common format in backbone and refinement
	cn_bb_common_format = solformat_collapse(cn_bb_common_format);

	if (cn_bb_common_format)
	{
		heur->backbone_sol_format = cn_bb_common_format;
		heur->construct_sol_format = cn_bb_common_format;
	}

	// make a final collapse in order to avoid multiple formats
	heur->backbone_sol_format = solformat_collapse(heur->backbone_sol_format);
	heur->construct_sol_format = solformat_collapse(heur->construct_sol_format);
	heur->refine_sol_format = solformat_collapse(heur->refine_sol_format);
	
}

/* **************************************************************************************************
*				ADD A NODE TO THE SOLUTION'S TOUR THROUGH THE extramileage PROCEDURE
************************************************************************************************** */
char extramileage_move(Solution* sol, graph* g, GraspData* grasp)
{
	// get unvisited nodes
	int* unvisited_nodes = NULL;	arr_malloc_s(unvisited_nodes, sol->nnodes, int);
	int unvisited_nodes_num = compute_list_unvisited_nodes(sol, unvisited_nodes);

	// if an extramileage move can be made
	if (unvisited_nodes_num > 0)
	{
		// setup iterator
		SolutionIterator iter;
		initialize_sol_iterator(&iter, sol, sol->handle_node);

		// ******************************* ITERATE THROUGH *******************************
		// allocate cand pool if using grasp
		IndexedValue* candpool = NULL;
		int candpool_size = grasp != NULL ? grasp->candpool_size : 1;
		arr_malloc_s(candpool, candpool_size, IndexedValue);

		// initialize pool of candidates
		OIA_clear(candpool, candpool_size);

		int prev = iter.curr;
		int next;
		int ext;
		char next_ready;
		double delta;
		// iterate each edge of the tour
		do
		{
			next_ready = next_node_in_solution(&iter);
			next = iter.curr;
			// iterate unvisited nodes
			for (int i = 0; i < unvisited_nodes_num; i++)
			{
				ext = unvisited_nodes[i];
				// if the triplet (prev, next, ext) has a low enough delta cost
				if (OIA_eligible(candpool, (delta = delta_cost(prev, next, ext, g))))
				{
					// insert it
					OIA_insert(candpool, candpool_size, OIA_pack3(prev, next, ext, delta));
				}
			}
			prev = next;

		} while (next_ready);

		IndexedValue cand_chosen;
		if (grasp->p_dev)
		{
			cand_chosen = random() < grasp->p_dev ?
				OIA_choose(candpool, candpool_size, 0) :
				OIA_best(candpool, candpool_size);
		}
		else cand_chosen = *candpool;

		int
			a = cand_chosen.index.arr[0],
			b = cand_chosen.index.arr[1],
			c = cand_chosen.index.arr[2];

		log_line_ext(VERBOSITY, LOGLVL_PEDANTIC, "Found triangle %d, %d, %d with delta cost = %f",
			a, b, c, cand_chosen.value);

		rem_edge_solution(sol, a, b, dist(a, b, g));
		add_edge_solution(sol, a, c, dist(a, c, g));
		add_edge_solution(sol, c, b, dist(c, b, g));

		free(candpool);
	}

	free(unvisited_nodes);

	return unvisited_nodes_num > 0;
}

double move_2opt(int* succ, graph* g, char allow_unimproving)
{
	// define edges ab and cd as those to replace
	// define edges ac and bd as those to insert
	double ab_dist, cd_dist, ac_dist, bd_dist;
	int a, b, c, d;
	int a_min = 0, b_min = 0, c_min = 0, d_min = 0;
	double delta_min = INFINITY;
	double delta_curr;

	a = 0;
	// for each node of the tour
	for (int i = 0; i < g->nnodes; i++)
	{
		b = succ[a];
		ab_dist = dist(a, b, g);
		
		// c starts after b
		c = succ[b];
		// for each remaining edge in the tour
		// apart from (?, a), (a,b), (b, ?)
		for (int j = 0; j < g->nnodes - 3; j++)
		{
			d = succ[c];
			cd_dist = dist(c, d, g);
			ac_dist = dist(a, c, g);
			bd_dist = dist(b, d, g);

			delta_curr = (ac_dist + bd_dist) - (ab_dist + cd_dist);
			if (delta_curr < delta_min)
			{
				delta_min = delta_curr;
				a_min = a;
				b_min = b;
				c_min = c;
				d_min = d;
			}
			c = d;
			d = succ[d];

		}
		// get to the next edge
		a = b;
		b = succ[b];
	}

	// if this is an improving move (delta<0)
	// or non improving moves are allowed
	if (delta_min < 0 || allow_unimproving)
	{
		// remember next of b_min
		int new_next = b_min;
		int new_curr = succ[b_min];
		int new_prev;

		succ[a_min] = c_min;
		succ[b_min] = d_min;

		// reverse path from b_min to c_min
		do
		{
			// swap direction
			new_prev = succ[new_curr];
			succ[new_curr] = new_next;

			// shift
			new_next = new_curr;
			new_curr = new_prev;

		} while (new_next != c_min);


		return delta_min;
	}

	return 0;
}

/* **************************************************************************************************
*					BACKBONE METHOD: ITERATIVE LOCAL SEARCH
*			- construct a starting solution
*			- refine iteratively the solution until time expires (or if needed until
*				a local optimum is found)
************************************************************************************************** */
void backbone_iter_local_search(OptData* optdata, Solution* sol, Heuristic* heur, void* data)
{
	// ******************************** SETUP ********************************
	// unpack
	instance* inst = optdata->inst;
	IterativeData* bb_data = (IterativeData*)data;

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
	if (!convert_solution(sol, heur->refine_sol_format)) print_error(ERR_HEUR_INFEASABLE_SOL, NULL);

	// iterate until time expires or the refining heuristic cannot find a better solution (local opt)
	while (!time_limit_expired(inst) && (heur->to_timelimit || !sol->optimal))
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

	if (sol->optimal) log_line(VERBOSITY, LOGLVL_INFO, "[INFO]: Local optimum was found!");

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
	// extract stat and look for optimality
	int stat = CPXgetstat(env, lp);
	sol->optimal = stat == CPXMIP_OPTIMAL;

	// ********************** EXTRACT **********************
	// get the optimal solution
	mip_extract_sol_obj(optdata, sol, "TSP solver");

}

/* **************************************************************************************************
*					CONSTRUCTIVE HEURISTIC: GREEDY METHOD
*			- starting from each node, find the tour following the greedy choice:
*				"choose the closest non-visited node"
*			- return the best solution
************************************************************************************************** */
void construct_greedy(OptData* optdata, Solution* sol, void* data, double timelim)
{
	// ************************************ SETUP ************************************
	// start time
	double start_time = second();
	
	// unpack optdata
	instance* inst = optdata->inst;
	graph* g = &inst->inst_graph;

	// unpack cnstdata
	GreedyAlgoData* cnstdata = (GreedyAlgoData*)data;
	double p_dev = cnstdata->grasp.p_dev;
	int candpool_size = cnstdata->grasp.candpool_size;
	if (p_dev) log_line(VERBOSITY, LOGLVL_INFO, "[INFO] Running greedy method with GRASP");
	else log_line(VERBOSITY, LOGLVL_INFO, "[INFO] Running standard greedy (no GRASP)");

	// setup best and current solutions
	Solution* best_sol = sol;
	sol->cost = INFINITY;
	Solution st_curr_sol;
	Solution* curr_sol = &st_curr_sol;
	empty_solution(curr_sol, sol->format, sol->nnodes);
	Solution* swap_sol;

	// ******************************* ITERATE THROUGH *******************************
	char* visited = NULL;			arr_malloc_s(visited, g->nnodes, char);
	IndexedValue* candpool = NULL;	arr_malloc_s(candpool, candpool_size, IndexedValue);
	IndexedValue cand_chosen;
	double h_dist;
	int curr;
	// for every start point (or until time expires if using GRASP)
	for (int start = 0; start < g->nnodes || p_dev > 0; start++)
	{
		// restart start if needed
		start %= g->nnodes;
		// select curr
		curr = start;
		// reset visited nodes
		for (int i = 0; i < g->nnodes; i++) visited[i] = 0;
		// clear solution
		clear_solution(curr_sol);

		// for every nodes in the path
		for (int n = 0; n < g->nnodes-1; n++)
		{
			// set current as visited
			visited[curr] = 1;

			// initialize pool of candidates
			OIA_clear(candpool, candpool_size);

			// search for candidate nodes for curr
			for (int h = 0; h < g->nnodes; h++)
			{
				// if not visited and is eligible to be a candidate
				if (!visited[h] && OIA_eligible(candpool, (h_dist = dist(curr, h, g))))
				{
					// insert h in the pool of candidates
					OIA_insert(candpool, candpool_size, OIA_pack1(h, h_dist));
				}
			}

			// choose whether to make the greedy choice or allow some deviation
			cand_chosen = random() < p_dev ?
				OIA_choose(candpool, candpool_size, 0) :
				OIA_best(candpool, candpool_size);

			// set closest node as successor of curr
			add_edge_solution(curr_sol, curr, cand_chosen.index.arr[0], cand_chosen.value);
			// go on with current node
			curr = cand_chosen.index.arr[0];
		}
		// link last node with start
		add_edge_solution(curr_sol, curr, start, dist(curr, start, g));

		// check if current solution is the best until now, if so swap them
		if (curr_sol->cost < best_sol->cost)
		{
			swap_sol = best_sol;
			best_sol = curr_sol;
			curr_sol = swap_sol;
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Found a new best solution with cost %f", best_sol->cost);
		}
		
		// if timelimit is reached, break
		if (second() - start_time >= timelim) break;
	}

	// copy the best known solution into the final solution (freeing what is not needed)
	if (sol != best_sol) shallow_copy_solution(sol, best_sol);
	else // free the other solution
	{
		free(curr_sol->succ); free(curr_sol->xstar);
	}

	// *********************************** CLEANUP ***********************************
	// cleanup
	free(visited);
	free(candpool);

}

/* **************************************************************************************************
*					CONSTRUCTIVE HEURISTIC: EXTRA MILEAGE METHOD
*			- build a starting incomplete solution with a tour
*			- add each external node c by replacing a tour edge (a,b) with (a,c) and (c,b)
*			- select node c which minimizes:
*					delta_cost(a,b,c) = dist(a,c) + dist(c,b) - dist(a,b)
*			- if randomizing, compute more solutions
*			- return the best solution
************************************************************************************************** */
void construct_extramileage(OptData* optdata, Solution* sol, void* data, double timelim)
{
	// ************************************ SETUP ************************************
	// start time
	double start_time = second();

	// unpack optdata
	instance* inst = optdata->inst;
	graph* g = &inst->inst_graph;

	// unpack cnstdata
	ExtraMileageData* cnstdata = (ExtraMileageData*)data;
	double p_dev = cnstdata->grasp.p_dev;
	int candpool_size = cnstdata->grasp.candpool_size;
	if (p_dev) log_line(VERBOSITY, LOGLVL_INFO, "[INFO] Running extramileage method with GRASP");
	else log_line(VERBOSITY, LOGLVL_INFO, "[INFO] Running standard extramileage (no GRASP)");
	char variant_startingchoice = cnstdata->variant_startingchoice;

	// setup best and current solutions
	Solution* best_sol = sol;
	sol->cost = INFINITY;
	Solution st_curr_sol;
	Solution* curr_sol = &st_curr_sol;
	// setup solution as an empty solution (no edges)
	empty_solution(curr_sol, sol->format, sol->nnodes);
	Solution* swap_sol;

	do
	{
		// ************************** BUILD INITIAL CONFIGURATION ************************
		// clear solution
		clear_solution(curr_sol);
		int a, b;
		double ab_dist;
		switch (variant_startingchoice)
		{
		case 'h':	// hull
			build_convex_hull(g, curr_sol);
			break;
		case 'r':	// random
			a = (int)(random() * g->nnodes);
			b = (int)(random() * g->nnodes);
			if (a == b) b = (b + 1) % g->nnodes;
			ab_dist = dist(a, b, g);
			add_edge_solution(curr_sol, a, b, ab_dist);
			curr_sol->handle_node = a;
			break;
		case 'f':	// furthest
			furthest_nodes_sol(g, curr_sol);
			break;
		default:
			print_error_ext(ERR_HEUR_CONSTR_PARAM, "variant, got  %c, should be h/r/f", variant_startingchoice);
			break;

		}
		// plot solution if needed
		if (VERBOSITY >= LOGLVL_DEBUGPLOT_P) plot_tsp_solution_undirected(g, curr_sol);

		// ******************************** FILL SOLUTION ********************************

		// for every remaining node, do an extra milage move
		int i = 0;
		int log_period = 100;
		do
		{
			i++;
			if (VERBOSITY >= LOGLVL_DEBUGPLOT_P && i % log_period == 0) plot_tsp_solution_undirected(g, curr_sol);

		} while (extramileage_move(curr_sol, g, &cnstdata->grasp));

		// check if current solution is the best until now, if so swap them
		if (curr_sol->cost < best_sol->cost)
		{
			swap_sol = best_sol;
			best_sol = curr_sol;
			curr_sol = swap_sol;
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Found a new best solution with cost %f", best_sol->cost);
		}

		// if timelimit is reached, break
		if (second() - start_time >= timelim) break;

	} while (p_dev > 0.0);

	// copy the best known solution into the final solution (freeing what is not needed)
	if (sol != best_sol) shallow_copy_solution(sol, best_sol);
	else // free the other solution
	{
		free(curr_sol->succ); free(curr_sol->xstar);
	}
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
		fix_ratio = (random() * (values[1] - values[0])) + values[0];
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
			idx = (int)(random() * (nnodes - i)) + i;
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

	// extract stat and look for optimality
	int stat = CPXgetstat(env, lp);
	// solution is optimal if cplex returned an optimal solution for the whole problem
	sol->optimal = stat == CPXMIP_OPTIMAL && n_fix == 0;

	// ********************** EXTRACT SOLUTION AND UPDATE METHOD *********************
	// extract best solution (expected to be already allocated in the construction phase!)
	double obj_before = sol->cost;
	mip_extract_sol_obj(optdata, sol, "Hardfixing");
	double obj_after = sol->cost;
	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Solution improved by %f%%", (obj_before - obj_after) / obj_before);

	// if using schedule and obj value is not good enough
	if (obj_after >= obj_before * IMPROVEMENT_RATIO)
	{
		rfn_data->rounds_nonimprove++;
		log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Rounds of non-improving solution = %d", rfn_data->rounds_nonimprove);
		// if exceeded rounds of non improving solutions
		if (rfn_data->rounds_nonimprove >= rfn_data->max_rounds)
		{
			rfn_data->rounds_nonimprove = 0;
			if (variant == 's')
			{
				// if patience ran out and progress is ended, make it optimal
				if (rfn_data->progress == HEUR_HARDFIX_VALUES - 1) sol->optimal = 1;
				// increase progress in schedule (capped)
				rfn_data->progress = min(progress + 1, HEUR_HARDFIX_VALUES - 1);
				log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Hardfix schedule progress to %d", rfn_data->progress);

			}
			else
			{
				// if patience ran out and progress is not provided, make it optimal
				sol->optimal = 1;
			}
		}
	}
	else // reset counter
	{
		log_line(VERBOSITY, LOGLVL_INFO, "[INFO] Non-improving solution rounds reset to 0");
		rfn_data->rounds_nonimprove = 0;
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


/* **************************************************************************************************
*					REFINING HEURISTIC: LOCAL BRANCHING
*			- add the k-opt moves neighborhood constraint
*			- re-optimize the simplified problem using cplex
*			- remove the constraint
************************************************************************************************** */
void refine_localbranching(OptData* optdata, Solution* sol, void* data, double timelim)
{
	// ************************************ SETUP ************************************
	// unpack optdata
	instance* inst = optdata->inst;
	CPXENVptr env = optdata->cpx->env;
	CPXLPptr lp = optdata->cpx->lp;

	// unpack rfndata
	LocalbranchingData* rfn_data = (LocalbranchingData*)data;
	int* values = rfn_data->values;
	int progress = rfn_data->progress;

	// extract data
	int nnodes = sol->nnodes;

	// extract current solution
	double* xstar = sol->xstar;
	int ncols = inst->inst_model.ncols;

	// extract number k of edges to replace
	int k_replace = values[progress];
	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] %d free variables", k_replace);

	// ******************** BUILD CONSTRAINT FOR K-OPT NEIGHBORHOOD ******************
	// add constraint: sum of active edges >= n-k
	int nnz = inst->inst_graph.nnodes;
	double rhs = inst->inst_graph.nnodes - k_replace;
	char sense = 'G';
	int pos = 0;
	int* index;		arr_malloc_s(index, nnz, int);
	double* value;	arr_malloc_s(value, nnz, double);

	// iterate through each edge and find those which are active
	for (int i = 0; i < ncols; i++)
	{
		if (is_one(xstar[i]))
		{
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Active edge of idx: %d", i);
			if (pos >= nnz) print_error(ERR_HEUR_INFEASABLE_SOL, NULL);
			index[pos] = i;
			value[pos] = 1.0;
			pos++;
		}
	}

	// add the constraint
	int constr_idx = CPXgetnumrows(env, lp);
	mip_add_cut(env, lp, nnz, rhs, sense, index, value, CUT_STATIC, "k-OPT neighborhood constraint", -1);

	// ******************************** OPTIMIZATION *********************************
	// optimize with cplex
	int error;
	mip_timelimit(optdata, timelim);
	mip_warmstart(optdata, xstar);
	CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_BALANCED);
	if (error = CPXmipopt(env, lp))
		print_error_ext(ERR_CPLEX, "CPXmipopt() Localbranching, CPX error: %d", error);

	// extract stat and look for optimality
	int stat = CPXgetstat(env, lp);
	// solution is local optimal if cplex returned an optimal solution for the whole problem
	sol->optimal = stat == CPXMIP_OPTIMAL && k_replace == inst->inst_graph.nnodes;

	// ********************** EXTRACT SOLUTION AND UPDATE METHOD *********************
	// extract best solution (expected to be already allocated in the construction phase!)
	double obj_before = sol->cost;
	mip_extract_sol_obj(optdata, sol, "Localbranching");
	double obj_after = sol->cost;
	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Solution improved by %f%%", (obj_before - obj_after) / obj_before);

	// if using schedule and obj value is not good enough
	if (obj_after >= obj_before * IMPROVEMENT_RATIO)
	{
		rfn_data->rounds_nonimprove++;
		log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Rounds of non-improving solution = %d", rfn_data->rounds_nonimprove);
		// if exceeded rounds of non improving solutions
		if (rfn_data->rounds_nonimprove >= rfn_data->max_rounds)
		{
			rfn_data->rounds_nonimprove = 0;

			// if patience ran out and progress is ended, make it optimal
			if (rfn_data->progress == HEUR_LOCALBRANCH_VALUES - 1) sol->optimal = 1;
			// increase progress in schedule (capped)
			rfn_data->progress = min(progress + 1, HEUR_LOCALBRANCH_VALUES - 1);
			log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Local branching schedule progress to %d", rfn_data->progress);

		}
	}
	else // reset counter
	{
		log_line(VERBOSITY, LOGLVL_INFO, "[INFO] Non-improving solution rounds reset to 0");
		rfn_data->rounds_nonimprove = 0;
	}

	// ******************************* REMOVE CONSTRAINT *****************************
	if (error = CPXdelrows(env, lp, constr_idx, constr_idx))
		print_error_ext(ERR_CPLEX, "CPXdelrows() Localbranching, CPX error: %d", error);
	log_line(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Removed k-opt constraint");

	// *********************************** CLEANUP ***********************************
	// cleanup
	free(index);
	free(value);

}


void refine_2opt(OptData* optdata, Solution* sol, void* data, double timelim)
{
	// compute move
	double delta = move_2opt(sol->succ, &optdata->inst->inst_graph, 0);
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Solution cost changed by: %f through a 2-OPT move", delta);
	// change cost
	sol->cost += delta;

	// if not improving
	if (delta >= 0)
	{
		// mark solution as local optimum
		sol->optimal = 1;
		log_line(VERBOSITY, LOGLVL_INFO, "[INFO] Local minimum in 2-OPT neighborhood reached");	
	}

}