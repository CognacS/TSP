#include "../include/heuristics.h"

/* **************************************************************************************************
*						HEURISTIC CODE SECTION
************************************************************************************************** */
// structure of whole code:
//	s1s2s3
// structure of a single section:
//	m(p1,p2,p3,...)
void decompose_heuristic_code(HeuristicCode* decomp, char* heuristic_code)
{
	char* section_token = NULL;		char* section_context = NULL;
	char* param_token = NULL;		char* param_context = NULL;
	
	char section_sep[] = ")";
	char param_sep[] = ",";

	// start sections tokenization
	section_token = strtok_u(heuristic_code, section_sep, &section_context);
	// for each section
	int s = 0;
	for (s = 0; (s < HEUR_SECTIONS) && (section_token != NULL); s++)
	{
		// extract method
		decomp->methods[s] = section_token[0];

		// start parameters tokenization
		param_token = strtok_u(section_token+2, param_sep, &param_context);
		// for each parameter
		for (int p = 0; (p < HEUR_NUMPARAM_MAX) && (param_token != NULL); p++)
		{
			// extract parameter
			decomp->params[s][p] = param_token;
			// next param token
			param_token = strtok_u(NULL, param_sep, &param_context);
		}

		// next section token
		section_token = strtok_u(NULL, section_sep, &section_context);
	}
	// iff missing some sections, set the method as 0
	if (s < HEUR_SECTIONS) decomp->methods[s] = 0;
}

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
	CallbackInstance cb_inst;
	if (heur.requires_cplex)
	{
		log_line(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Cplex needed: setting up"),

		// ********************** SETUP BASE MODEL **********************
		// setup cplex for use
		cpx_setup_cplex(optdata);

		// build naive model
		build_model_base_undirected(optdata);

		// ********************** BUILD DATASTRUCTURES **********************
		// unpack
		Instance* inst = optdata->inst;
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
		optdata->perflog[0] = LL_new();
		for (int i = 1; i < MAX_PERFLOGS; i++)
			optdata->perflog[i] = NULL;
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
	GlobalData* gd = &optdata->inst->inst_global_data;
	gd->xstar = final_sol.xstar;
	gd->zbest = final_sol.cost;

	// if need to plot the progress during heuristic
	if (VERBOSITY >= LOGLVL_PLOTSOL)
	{
		for (int i = 0; i < MAX_PERFLOGS; i++)
		{
			if (optdata->perflog[i] != NULL)
			{
				if (!LL_is_empty(optdata->perflog[i])) plot_heuristic_perflog(optdata->perflog[i]);
				LL_free(optdata->perflog[i]);
				optdata->perflog[i] = NULL;
			}
		}
		
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
	Instance* inst = optdata->inst;
	CPXENVptr env = optdata->cpx->env;
	CPXLPptr lp = optdata->cpx->lp;
	char* code = inst->inst_params.heuristic_code;

	// decompose heuristic code
	HeuristicCode decomp_code;
	decompose_heuristic_code(&decomp_code, code);
	HeuristicSection s;
	char m;

	// setup heuristic
	heur->requires_cplex = 0;
	heur->to_timelimit = 1;

	// ******************************* DECODE BACKBONE *******************************
	s = HEUR_BACKBONE;
	m = decomp_code.methods[s];
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "Backbone method: %c", m);

	// construct time
	heur->construct_timelimit = atof(decomp_code.params[s][0]);
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "construct_time: %f", heur->construct_timelimit);

	// refine time
	heur->refine_timelimit = atof(decomp_code.params[s][1]);
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "refine_time: %f", heur->refine_timelimit);

	if (m == 'i')			// iterative local search
	{
		heur->backbone_heur = backbone_iter_local_search;

		// assign data
		heur->backbone_data = NULL;
		heur->backbone_sol_format = SOLFORMAT_XSTAR | SOLFORMAT_SUCC;
		heur->to_timelimit = 1;
	}
	else if (m == 'v')		// variable neighborhood search
	{
		heur->backbone_heur = backbone_var_neighborhood_search;
		VNSData* data = NULL;	malloc_s(data, VNSData);

		// max_k
		data->max_k = atoi(decomp_code.params[s][2]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "max_k: %d", data->max_k);

		// assign data
		heur->backbone_data = data;
		heur->backbone_sol_format = SOLFORMAT_SUCC;
		heur->to_timelimit = 1;
	}
	else if (m == 't')		// tabu search
	{
		heur->backbone_heur = backbone_tabu_search;
		TabuSearchData* data = NULL;	malloc_s(data, TabuSearchData);

		// min_tenure
		data->min_tenure = atoi(decomp_code.params[s][2]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "min_tenure: %d", data->min_tenure);

		// max_tenure
		data->max_tenure = atoi(decomp_code.params[s][3]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "max_tenure: %d", data->max_tenure);

		// period
		data->period = atoi(decomp_code.params[s][4]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "period: %d", data->period);

		// interm_levels
		data->interm_levels = atoi(decomp_code.params[s][5]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "interm_levels: %d", data->interm_levels);

		// assign data
		heur->backbone_data = data;
		heur->backbone_sol_format = SOLFORMAT_SUCC;
		heur->to_timelimit = 1;
	}
	else if (m == 'g')		// genetic algorithm
	{
		heur->backbone_heur = backbone_genetic_algorithm;
		GeneticAlgData* data = NULL;	malloc_s(data, GeneticAlgData);

		// crossover_variant
		data->crossover_variant = *decomp_code.params[s][2];
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "crossover_variant: %c", data->crossover_variant);

		// pool_size
		data->pool_size = atoi(decomp_code.params[s][3]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "pool_size: %d", data->pool_size);

		// elite_ratio
		data->elite_ratio = atof(decomp_code.params[s][4]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "elite_ratio: %f", data->elite_ratio);

		// mutation_p
		data->mutation_p = atof(decomp_code.params[s][5]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "mutation_p: %f", data->mutation_p);

		// assign data
		heur->backbone_data = data;
		heur->construct_data = NULL;
		heur->refine_data = NULL;
		heur->backbone_sol_format = SOLFORMAT_CHROMO;
		heur->construct_sol_format = SOLFORMAT_CHROMO;
		heur->refine_sol_format = SOLFORMAT_CHROMO;
		heur->to_timelimit = 1;
	}
	else
	{
		print_error_ext(ERR_HEUR_DECODE, "unknown backbone method: %c", m);
	}

	// ******************************* DECODE CONSTRUCT ******************************
	s = HEUR_CONSTRUCT;
	m = decomp_code.methods[s];
	if (m == 0) return;
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "Construct method: %c", m);

	if (m == 's')			// construct through cplex solver
	{
		heur->construct_heur = construct_tspsolver;

		// assign data
		heur->construct_data = NULL;
		heur->requires_cplex = 1;
		heur->construct_sol_format = SOLFORMAT_XSTAR;
	}
	else if (m == 'g')		// construct through greedy algorithm
	{
		heur->construct_heur = construct_greedy;
		GreedyAlgoData* data = NULL;  malloc_s(data, GreedyAlgoData);

		// deviation probability
		data->grasp.p_dev = atof(decomp_code.params[s][0]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "p_dev: %f", data->grasp.p_dev);

		// candidate pool size
		data->grasp.candpool_size = atoi(decomp_code.params[s][1]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "candpool_size: %d", data->grasp.candpool_size);

		// assign data
		heur->construct_data = data;
		heur->construct_sol_format = SOLFORMAT_XSTAR | SOLFORMAT_SUCC;
	}

	else if (m == 'e')		// construct through extra mileage algorithm
	{
		heur->construct_heur = construct_extramileage;
		ExtraMileageData* data = NULL;	malloc_s(data, ExtraMileageData);

		// deviation probability
		data->grasp.p_dev = atof(decomp_code.params[s][0]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "p_dev: %f", data->grasp.p_dev);

		// candidate pool size
		data->grasp.candpool_size = atoi(decomp_code.params[s][1]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "candpool_size: %d", data->grasp.candpool_size);

		// variant for starting edges choice
		data->variant_startingchoice = *decomp_code.params[s][2];
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "variant_startingchoice: %c", data->variant_startingchoice);

		// assign data
		heur->construct_data = data;
		heur->construct_sol_format = SOLFORMAT_SUCC;
	}
	else
	{
		print_error_ext(ERR_HEUR_DECODE, "unknown constructive heuristic: %c", m);
	}

	// ******************************** DECODE REFINE ********************************
	s = HEUR_REFINE;
	m = decomp_code.methods[s];
	if (m == 0) return;
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "Refine method: %c", m);

	if (m == 'h')			// refine through hardfixing scheme
	{
		heur->refine_heur = refine_hardfixing;
		HardfixingData* data = NULL;	malloc_s(data, HardfixingData);

		// variant
		data->variant = *decomp_code.params[s][0];
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "variant: %c", data->variant);

		// max_rounds
		data->max_rounds = atoi(decomp_code.params[s][1]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "max_rounds: %d", data->max_rounds);

		// values
		for (int i = 0; i < HEUR_HARDFIX_VALUES; i++)
		{
			data->values[i] = atof(decomp_code.params[s][2+i]);
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "values[%d]: %f", i, data->values[i]);
		}

		// assign data
		data->rounds_nonimprove = 0;
		data->progress = 0;
		data->num_values = HEUR_HARDFIX_VALUES;
		heur->refine_data = data;
		heur->requires_cplex = 1;
		heur->refine_sol_format = SOLFORMAT_XSTAR;
	}
	else if (m == 'l')		// refine through local branching scheme
	{
		heur->refine_heur = refine_localbranching;
		LocalbranchingData* data = NULL;	malloc_s(data, LocalbranchingData);

		// max_rounds
		data->max_rounds = atoi(decomp_code.params[s][0]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "max_rounds: %d", data->max_rounds);

		// values
		for (int i = 0; i < HEUR_LOCALBRANCH_VALUES; i++)
		{
			data->values[i] = atoi(decomp_code.params[s][1 + i]);
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "values[%d]: %f", i, data->values[i]);
		}

		// assign data
		data->rounds_nonimprove = 0;
		data->progress = 0;
		data->num_values = HEUR_LOCALBRANCH_VALUES;
		heur->refine_data = data;
		heur->requires_cplex = 1;
		heur->refine_sol_format = SOLFORMAT_XSTAR;
	}
	else if (m == 'm')		// refine through a mix of hardfix and local branching
	{
		heur->refine_heur = refine_mixhardsoft;
		MixHardSoftData* data = NULL;	malloc_s(data, MixHardSoftData);

		// max_rounds
		data->max_rounds = atoi(decomp_code.params[s][0]);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "max_rounds: %d", data->max_rounds);

		// values
		for (int i = 0; i < HEUR_MIXHARDSOFT_VALUES; i++)
		{
			data->values[i] = atof(decomp_code.params[s][1 + i]);
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "values[%d]: %f", i, data->values[i]);
		}

		// assign data
		data->rounds_nonimprove = 0;
		data->progress = 0;
		data->num_values = HEUR_MIXHARDSOFT_VALUES;
		heur->refine_data = data;
		heur->requires_cplex = 1;
		heur->refine_sol_format = SOLFORMAT_XSTAR;
	}
	else if (m == '2')		// refine through 2-OPT moves
	{
		heur->refine_heur = refine_2opt;
		Opt2MoveData* data = NULL;	malloc_s(data, Opt2MoveData);
		data->tabu = NULL;
		data->allow_worsening = 0;

		// assign data
		heur->refine_data = data;
		heur->refine_sol_format = SOLFORMAT_SUCC;
	}
	else
	{
		print_error_ext(ERR_HEUR_DECODE, "unknown refining heuristic: %c", m);
	}
}

/* **************************************************************************************************
*					RESOLVE CONVERSIONS INSIDE THE BACKBONE HEURISTIC METHOD
************************************************************************************************** */
void resolve_solformats(Heuristic* heur)
{
	// emphasis on common format between backbone and refinement
	SolFormat bb_rf_common_format =
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
	SolFormat cn_bb_common_format =
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
char extramileage_move(Solution* sol, Graph* g, SetOfNodes* ext_nodes, GraspData* grasp)
{

	// if an extramileage move can be made (nodes can be added to the tour)
	char move_doable = !SETN_isempty(ext_nodes);
	if (move_doable)
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
			for (int i = 0; i < ext_nodes->curr_size; i++)
			{
				// extract a node from the set of unvisited nodes
				ext = SETN_get(ext_nodes, i);
				// if the triplet (prev, next, ext) has a low enough delta cost
				if (OIA_eligible(candpool, (delta = delta_cost(g, prev, next, ext))))
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

		rem_edge_solution(sol, a, b, dist(g, a, b));
		add_edge_solution(sol, a, c, dist(g, a, c));
		add_edge_solution(sol, c, b, dist(g, c, b));

		free(candpool);

		// remove chosen external node c from the set
		SETN_remove(ext_nodes, c);

	}

	return move_doable;
}

/* **************************************************************************************************
*				KICK SOLUTION IN A RANDOM POSITION IN ITS K-OPT NEIGHBORHOOD:
*			- graph with N nodes
*			- select randomly k<=N nodes: n1, n2,... nk
*			- pick their respective succ: s1, s2,... sk
*			- shift all succ such that succ[ni] = s_{(i+1) mod N}
************************************************************************************************** */
void kopt_kick(Solution* sol, Graph* g, int k)
{
	int nnodes = g->nnodes;

	// choose k edges to be removed
	SetOfNodes* chosen_nodes = SETN_new(k, nnodes);
	int* succ_nodes = NULL;		arr_malloc_s(succ_nodes, k, int);
	int node;

	for (int i = 0; i < k; i++)
	{
		// continue generating random nodes until it can be added
		do
		{
			node = (int)(random() * nnodes);
		} while (!SETN_add(chosen_nodes, node));

		log_line_ext(VERBOSITY, LOGLVL_PEDANTIC, "[PEDANTIC] Chosen edge %d->%d to change", node, sol->succ[node]);
		
	}

	// reorder chosen nodes in tour order of appearance
	int curr = 0;
	int nodes_found = 0;
	for (int i = 0; i < nnodes; i++)
	{
		// if curr node was chosen, then swap it with the first in line
		if (SETN_exists(chosen_nodes, curr))
		{
			// reposition node on top
			SETN_reposition(chosen_nodes, curr, nodes_found);
			// also set succ
			succ_nodes[nodes_found] = sol->succ[curr];
			nodes_found++;

			log_line_ext(VERBOSITY, LOGLVL_PEDANTIC, "[PEDANTIC] Found edge %d->%d on the way", curr, sol->succ[curr]);
		}

		curr = sol->succ[curr];
	}

	// accumulate cost and shift successors in reverse
	double cost_add = 0;
	for (int i = 0; i < k; i++)
	{
		node = SETN_get(chosen_nodes, i);

		log_line_ext(VERBOSITY, LOGLVL_PEDANTIC, "[PEDANTIC] changing %d->%d with %d->%d",
			node, sol->succ[node], node, succ_nodes[(k + i - 2) % k]);

		cost_add -= dist(g, node, sol->succ[node]);
		sol->succ[node] = succ_nodes[(k + i - 2) % k];
		cost_add += dist(g, node, sol->succ[node]);
	}

	// change overall cost of solution
	sol->cost += cost_add;

	free(succ_nodes);
	SETN_free(chosen_nodes);
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
	Instance* inst = optdata->inst;

	// startup solution and format it for construction
	create_solution(sol, heur->construct_sol_format, inst->inst_graph.nnodes);

	// ************************* CONSTRUCT STARTING SOLUTION *************************
	double timelimit, time_passed = 0, time_phase, t_start;
	// compute timelimit for constructive heuristic
	timelimit = min(residual_time(inst), heur->construct_timelimit);
	log_line(VERBOSITY, LOGLVL_INFO, "************ STARTING CONSTRUCTION PHASE ************");
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Construction should run for %f sec.s", timelimit);

	// construct starting heuristic solution
	t_start = second();
	// *********** CONSTRUCT ***********
	heur->construct_heur(optdata, sol, heur->construct_data, timelimit);
	if (time_limit_expired(inst))
	{
		log_line(VERBOSITY, LOGLVL_WARN, "[WARN] Could not construct a starting solution, returning without a result");
		sol->cost = INFINITY;
		return;
	}
	// *********************************
	time_phase = second() - t_start;
	time_passed += time_phase;

	// log
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Construction done in %f sec.s", time_phase);
	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO]: Obj value after construction = %f", sol->cost);
	if (VERBOSITY >= LOGLVL_PLOTSOL) LL_add_value(optdata->perflog[0], time_passed, sol->cost);
	if (VERBOSITY >= LOGLVL_DEBUGPLOT) plot_tsp_solution_undirected(&inst->inst_graph, sol);

	// ****************************** START ITERATIONS *******************************
	// convert solution to the refinement format
	if (!convert_solution(sol, heur->refine_sol_format)) print_error(ERR_HEUR_INFEASABLE_SOL, NULL);

	int epoch = 0;

	// iterate until time expires or the refining heuristic cannot find a better solution (local opt)
	while (!time_limit_expired(inst) && (heur->to_timelimit || !sol->optimal))
	{
		// ***************************** REFINE SOLUTION *****************************
		// compute timelimit for refinement heuristic
		timelimit = min(residual_time(inst), heur->refine_timelimit);
		log_line_ext(VERBOSITY, LOGLVL_INFO, "************ STARTING EPOCH %d ************", epoch++);
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
		log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO]: Obj value after refinement = %f", sol->cost);
		if (VERBOSITY >= LOGLVL_PLOTSOL) LL_add_value(optdata->perflog[0], time_passed, sol->cost);
		if (VERBOSITY >= LOGLVL_DEBUGPLOT) plot_tsp_solution_undirected(&inst->inst_graph, sol);
	}

	if (sol->optimal) log_line(VERBOSITY, LOGLVL_INFO, "[INFO]: Local optimum was found!");
	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Total number of epochs: %d", epoch);

}

/* **************************************************************************************************
*					BACKBONE METHOD: VARIABLE NEIGHBORHOOD SEARCH
*			- construct a starting solution
*			- refine iteratively the solution
*			- if the refined solution gets stuck in a local minima, shake it in a bigger
*				neighborhood than what is reachable by the refinement procedure
************************************************************************************************** */
void backbone_var_neighborhood_search(OptData* optdata, Solution* sol, Heuristic* heur, void* data)
{
	// ******************************** SETUP ********************************
	// unpack
	Instance* inst = optdata->inst;
	VNSData* bbdata = (VNSData*)data;

	// startup solution and format it for construction
	create_solution(sol, heur->construct_sol_format, inst->inst_graph.nnodes);

	// ************************* CONSTRUCT STARTING SOLUTION *************************
	double timelimit, time_passed = 0, time_phase, t_start;
	// compute timelimit for constructive heuristic
	timelimit = min(residual_time(inst), heur->construct_timelimit);
	log_line(VERBOSITY, LOGLVL_INFO, "************ STARTING CONSTRUCTION PHASE ************");
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
	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO]: Obj value after construction = %f", sol->cost);
	if (VERBOSITY >= LOGLVL_PLOTSOL) LL_add_value(optdata->perflog[0], time_passed, sol->cost);
	if (VERBOSITY >= LOGLVL_DEBUGPLOT) plot_tsp_solution_undirected(&inst->inst_graph, sol);

	// ****************************** START ITERATIONS *******************************
	// convert solution to the refinement format
	if (!convert_solution(sol, heur->backbone_sol_format)) print_error(ERR_HEUR_INFEASABLE_SOL, NULL);

	int epoch = 0;
	int kick_k = 3;

	// setup best and current solutions
	Solution* best_sol = sol;
	Solution st_curr_sol;
	Solution* curr_sol = &st_curr_sol;
	// initialize and copy best solution to the current solution
	create_solution(curr_sol, best_sol->format, best_sol->nnodes);
	deep_copy_solution(curr_sol, best_sol);

	// iterate until time expires or the refining heuristic cannot find a better solution (local opt)
	while (!time_limit_expired(inst) && (heur->to_timelimit || !sol->optimal))
	{
		// ***************************** REFINE SOLUTION *****************************
		log_line_ext(VERBOSITY, LOGLVL_INFO, "************ STARTING EPOCH %d ************", epoch++);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "CORRECTNESS BEFORE REFINEMENT: %d", check_succ(curr_sol->succ, curr_sol->nnodes));
		// convert solution to the refinement format
		if (!convert_solution(curr_sol, heur->refine_sol_format)) print_error(ERR_HEUR_INFEASABLE_SOL, NULL);
		// compute timelimit for refinement heuristic
		timelimit = min(residual_time(inst), heur->refine_timelimit);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] 1 iteration of refinement should run for %f sec.s", timelimit);

		// refine heuristic solution
		t_start = second();
		// *********** REFINE ***********
		heur->refine_heur(optdata, curr_sol, heur->refine_data, timelimit);
		// ******************************
		time_phase = second() - t_start;
		time_passed += time_phase;

		// log
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "CORRECTNESS AFTER REFINEMENT: %d", check_succ(curr_sol->succ, curr_sol->nnodes));
		if (VERBOSITY >= LOGLVL_DEBUGPLOT_P) plot_tsp_solution_undirected(&inst->inst_graph, curr_sol);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] 1 iteration of refinement done in %f sec.s", time_phase);

		// ***************************** CHECK SOLUTION ******************************
		// convert solution to the refinement format
		if (!convert_solution(curr_sol, heur->backbone_sol_format)) print_error(ERR_HEUR_INFEASABLE_SOL, NULL);
		// check if current solution is the best until now, if so swap them
		if (curr_sol->cost < best_sol->cost)
		{
			// save curr solution as the new best solution
			deep_copy_solution(best_sol, curr_sol);
			// log new best solution
			log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO]: New best solution found, cost = %f", best_sol->cost);
			if (VERBOSITY >= LOGLVL_PLOTSOL) LL_add_value(optdata->perflog[0], time_passed, best_sol->cost);
			if (VERBOSITY >= LOGLVL_DEBUGPLOT) plot_tsp_solution_undirected(&inst->inst_graph, best_sol);
			// reset kick size
			kick_k = 3;
		}
		else
		{
			// copy back the best solution into the current solution
			deep_copy_solution(curr_sol, best_sol);
		}

		// ****************************** KICK SOLUTION ******************************
		// if the current solution is in a local optimum
		if (curr_sol->optimal)
		{
			// kick the solution in its k-opt neighborhood
			kopt_kick(curr_sol, &inst->inst_graph, kick_k);

			log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO]: Kicked solution in its %d-OPT neighborhood", kick_k);
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "CORRECTNESS AFTER K-OPT KICK: %d", check_succ(curr_sol->succ, curr_sol->nnodes));
			if (VERBOSITY >= LOGLVL_DEBUGPLOT_P) plot_tsp_solution_undirected(&inst->inst_graph, curr_sol);
			// increase k
			kick_k = min(kick_k+1, bbdata->max_k);
		}
	}

	// free the other solution
	free(curr_sol->succ); free(curr_sol->xstar); free(curr_sol->chromo);

	if (sol->optimal) log_line(VERBOSITY, LOGLVL_INFO, "[INFO]: Local optimum was found!");
	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Total number of epochs: %d", epoch);

}

/* **************************************************************************************************
*					BACKBONE METHOD: TABU SEARCH
*			- construct a starting solution
*			- refine iteratively the solution, allowing worsening moves
*			- impede going back to local minima by making "tabu" some moves
************************************************************************************************** */
void backbone_tabu_search(OptData* optdata, Solution* sol, Heuristic* heur, void* data)
{
	// ******************************** SETUP ********************************
	// unpack
	Instance* inst = optdata->inst;
	TabuSearchData* bbdata = (TabuSearchData*)data;

	// startup solution and format it for construction
	create_solution(sol, heur->construct_sol_format, inst->inst_graph.nnodes);

	// ************************* CONSTRUCT STARTING SOLUTION *************************
	double timelimit, time_passed = 0, time_phase, t_start;
	// compute timelimit for constructive heuristic
	timelimit = min(residual_time(inst), heur->construct_timelimit);
	log_line(VERBOSITY, LOGLVL_INFO, "************ STARTING CONSTRUCTION PHASE ************");
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
	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO]: Obj value after construction = %f", sol->cost);
	if (VERBOSITY >= LOGLVL_PLOTSOL) LL_add_value(optdata->perflog[0], time_passed, sol->cost);
	if (VERBOSITY >= LOGLVL_DEBUGPLOT) plot_tsp_solution_undirected(&inst->inst_graph, sol);

	// ****************************** START ITERATIONS *******************************
	// convert solution to the refinement format
	if (!convert_solution(sol, heur->refine_sol_format)) print_error(ERR_HEUR_INFEASABLE_SOL, NULL);

	int epoch = 0;

	// create tabu list and assign it to the refinement data
	TabuList* tabu = TABU_new(sol->nnodes, bbdata->min_tenure);
	((Opt2MoveData*)heur->refine_data)->tabu = tabu;
	// compute parameters
	int tenure_halfperiod = (bbdata->period / 2);
	int tenure_update_time = tenure_halfperiod / (bbdata->interm_levels + 2);
	int tenure_change = (bbdata->max_tenure - bbdata->min_tenure) / (bbdata->interm_levels + 1);

	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG]: halfperiod = %d; up_time = %d; change = %d",
		tenure_halfperiod, tenure_update_time, tenure_change);

	// setup best and current solutions
	Solution* best_sol = sol;
	Solution st_curr_sol;
	Solution* curr_sol = &st_curr_sol;
	// initialize and copy best solution to the current solution
	create_solution(curr_sol, best_sol->format, best_sol->nnodes);
	deep_copy_solution(curr_sol, best_sol);

	// iterate until time expires or the refining heuristic cannot find a better solution (local opt)
	while (!time_limit_expired(inst) && (heur->to_timelimit || !sol->optimal))
	{
		// ***************************** REFINE SOLUTION *****************************
		log_line_ext(VERBOSITY, LOGLVL_INFO, "************ STARTING EPOCH %d ************", epoch++);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "CORRECTNESS BEFORE REFINEMENT: %d", check_succ(curr_sol->succ, curr_sol->nnodes));
		// // convert solution to the refinement format
		if (!convert_solution(curr_sol, heur->refine_sol_format)) print_error(ERR_HEUR_INFEASABLE_SOL, NULL);
		// compute timelimit for refinement heuristic
		timelimit = min(residual_time(inst), heur->refine_timelimit);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] 1 iteration of refinement should run for %f sec.s", timelimit);

		// refine heuristic solution
		t_start = second();
		// *********** REFINE ***********
		heur->refine_heur(optdata, curr_sol, heur->refine_data, timelimit);
		// ******************************
		time_phase = second() - t_start;
		time_passed += time_phase;

		// log
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "CORRECTNESS AFTER REFINEMENT: %d", check_succ(curr_sol->succ, curr_sol->nnodes));
		if (VERBOSITY >= LOGLVL_DEBUGPLOT_P) plot_tsp_solution_undirected(&inst->inst_graph, curr_sol);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] 1 iteration of refinement done in %f sec.s", time_phase);

		// ***************************** CHECK SOLUTION ******************************
		// convert solution to the refinement format
		if (!convert_solution(curr_sol, heur->backbone_sol_format)) print_error(ERR_HEUR_INFEASABLE_SOL, NULL);
		// check if current solution is the best until now, if so swap them
		if (curr_sol->cost < best_sol->cost)
		{
			// save curr solution as the new best solution
			deep_copy_solution(best_sol, curr_sol);
			// log new best solution
			log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO]: New best solution found, cost = %f", best_sol->cost);
			if (VERBOSITY >= LOGLVL_PLOTSOL) LL_add_value(optdata->perflog[0], time_passed, best_sol->cost);
			if (VERBOSITY >= LOGLVL_DEBUGPLOT) plot_tsp_solution_undirected(&inst->inst_graph, best_sol);
		}

		// ****************************** UPDATE TABU VALUES ******************************
		if (curr_sol->optimal)
		{
			((Opt2MoveData*)heur->refine_data)->allow_worsening = 1;
			log_line(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Found an optimal solution => allow a worsening move");
		}

		if (VERBOSITY >= LOGLVL_DEBUGPLOT_P) TABU_print(tabu);
		if ((epoch+1) % tenure_halfperiod == 0)
		{ 
			// if hitting the half period, invert change
			tenure_change = -tenure_change;
			log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Tenure peak reached: %d", tabu->tenure);
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Inverted tenure change to %d", tenure_change);
		}
		else if ((epoch + 1) % tenure_update_time == 0)
		{
			if (VERBOSITY >= LOGLVL_PLOTSOL) LL_add_value(optdata->perflog[0], time_passed, (double)tabu->tenure * 10);
			int newtenure = min(bbdata->max_tenure, max(bbdata->min_tenure, tabu->tenure + tenure_change));
			TABU_set_tenure(tabu, newtenure);
			log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Tenure set to: %d", tabu->tenure);
			if (VERBOSITY >= LOGLVL_PLOTSOL) LL_add_value(optdata->perflog[0], time_passed, (double)tabu->tenure * 10);
		}
		
	}

	// free tabu list
	TABU_free(tabu);

	// free the other solution
	free(curr_sol->succ); free(curr_sol->xstar); free(curr_sol->chromo);

	if (sol->optimal) log_line(VERBOSITY, LOGLVL_INFO, "[INFO]: Local optimum was found!");
	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Total number of epochs: %d", epoch);

}

/* **************************************************************************************************
*					BACKBONE METHOD: GENETIC ALGORITHM
*			- produce a population of random feasable solutions
*			- generate new solutions through a mating and mutation process
************************************************************************************************** */
void backbone_genetic_algorithm(OptData* optdata, Solution* sol, Heuristic* heur, void* data)
{
	// ************************************ SETUP ************************************
	// add a log to account for averages
	optdata->perflog[1] = LL_new();
	
	// unpack
	Instance* inst = optdata->inst;
	Graph* g = &inst->inst_graph;
	GeneticAlgData* bbdata = (GeneticAlgData*)data;

	// startup solution and format it for chromosomes
	create_solution(sol, SOLFORMAT_CHROMO, g->nnodes);

	// compute parameters
	int pool_size = bbdata->pool_size;
	int elites_size = (int)(pool_size * bbdata->elite_ratio);
	double mutation_p = bbdata->mutation_p;

	void (*crossover)(int*, int*, int*, int) = NULL;
	switch (bbdata->crossover_variant)
	{
	case 'a':
		crossover = crossover_aex;
		break;
	case 'p':
		crossover = crossover_pmx;
		break;
	case 'm':
		crossover = crossover_mix;
		break;
	default:
		print_error_ext(ERR_HEUR_CONSTR_PARAM, "variant, got  %c, should be a/p/m", bbdata->crossover_variant);
	}

	// ***************************** INITIAL POPULATION ******************************
	// allocate population
	Population* pop = POP_new(bbdata->pool_size, elites_size, g->nnodes);
	double timelimit, time_passed = 0, time_phase, t_start;
	// compute timelimit for constructive heuristic
	timelimit = min(residual_time(inst), heur->construct_timelimit);
	log_line(VERBOSITY, LOGLVL_INFO, "************ STARTING POPULATION PHASE ************");

	// construct starting heuristic solution
	t_start = second();
	// *********** POPULATE ************
	populate(pop, g);
	// *********************************
	time_phase = second() - t_start;
	time_passed += time_phase;
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] population done in %f sec.s", time_phase);

	// compute statistics
	double best_fitness, avg_fitness;
	best_fitness = POP_fitness_best(pop); avg_fitness = POP_fitness_avg(pop);

	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Best fitness = %f", best_fitness);
	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Average fitness = %f", avg_fitness);
	if (VERBOSITY >= LOGLVL_PLOTSOL) LL_add_value(optdata->perflog[0], time_passed, best_fitness);
	if (VERBOSITY >= LOGLVL_PLOTSOL) LL_add_value(optdata->perflog[1], time_passed, avg_fitness);

	// ****************************** START ITERATIONS *******************************
	int epoch = 0;
	double long_range_prev_cost = INFINITY;

	while (!time_limit_expired(inst))
	{
		log_line_ext(VERBOSITY, LOGLVL_INFO, "************ GENERATION %d ************", epoch++);

		t_start = second();
		// *********** NEW GENERATION ***********
		new_generation(pop, crossover, g, mutation_p, residual_time(inst));
		// **************************************
		time_phase = second() - t_start;
		time_passed += time_phase;

		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] new generation in %f sec.s", time_phase);

		// compute statistics
		best_fitness = POP_fitness_best(pop); avg_fitness = POP_fitness_avg(pop);

		log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Best fitness = %f", best_fitness);
		log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Average fitness = %f", avg_fitness);
		//log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Entropy = %f", POP_entropy(pop));
		if (VERBOSITY >= LOGLVL_PLOTSOL && (epoch % 100 == 0)) LL_add_value(optdata->perflog[0], time_passed, best_fitness);
		if (VERBOSITY >= LOGLVL_PLOTSOL && (epoch % 100 == 0)) LL_add_value(optdata->perflog[1], time_passed, avg_fitness);
	}
	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Total number of generations: %d", epoch);
	log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Final fitness of champion: %f", POP_fitness_best(pop));

	// set champion as solution
	Specimen* champion = POP_getchampion(pop);
	sol->chromo = champion->chromo;
	champion->chromo = NULL;	// to avoid freeing the solution chromosome
	sol->cost = champion->cost;
	convert_solution(sol, SOLFORMAT_SUCC);

	// free datastructure
	POP_free(pop);
	
}

/* **************************************************************************************************
*					CONSTRUCTIVE HEURISTIC: CPLEX TSP SOLVER
*			- run the best cplex tsp solver for a fixed amount of time
*			- make cplex emphasize finding a good heuristic solution
************************************************************************************************** */
void construct_tspsolver(OptData* optdata, Solution* sol, void* data, double timelim)
{
	// unpack
	Instance* inst = optdata->inst;
	CPXENVptr env = optdata->cpx->env;
	CPXLPptr lp = optdata->cpx->lp;

	// ********************** OPTIMIZE **********************
	int solution_found = 0;
	while (!solution_found && !time_limit_expired(inst))
	{
		// solve the problem with the callback
		cpx_timelimit(optdata, timelim);
		CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_HEURISTIC);
		int error;
		if (error = CPXmipopt(env, lp))
			print_error_ext(ERR_CPLEX, "CPXmipopt() TSP solver, CPX error: %d", error);
		// extract stat and look for optimality
		int stat = CPXgetstat(env, lp);
		solution_found = stat != CPXMIP_TIME_LIM_INFEAS;
		sol->optimal = stat == CPXMIP_OPTIMAL;
		if (!solution_found) log_line(VERBOSITY, LOGLVL_INFO, "[INFO] No integer solution found, retrying");
		timelim *= 1.5;
	}

	// ********************** EXTRACT **********************
	// get the optimal solution
	if (solution_found) cpx_extract_sol_obj(optdata, sol, "TSP solver");

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
	Instance* inst = optdata->inst;
	Graph* g = &inst->inst_graph;

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
	create_solution(curr_sol, sol->format, sol->nnodes);
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
		erase_solution(curr_sol);

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
				if (!visited[h] && OIA_eligible(candpool, (h_dist = dist(g, curr, h))))
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
		add_edge_solution(curr_sol, curr, start, dist(g, curr, start));

		// check if current solution is the best until now, if so swap them
		if (curr_sol->cost < best_sol->cost)
		{
			swap_sol = best_sol;
			best_sol = curr_sol;
			curr_sol = swap_sol;
			log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Found a new best solution with cost %f", best_sol->cost);
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
	Instance* inst = optdata->inst;
	Graph* g = &inst->inst_graph;

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
	create_solution(curr_sol, sol->format, sol->nnodes);
	Solution* swap_sol;

	// ************************** BUILD INITIAL CONFIGURATION ************************
	Solution st_starting_sol;
	Solution* starting_sol = &st_starting_sol;
	create_solution(starting_sol, sol->format, sol->nnodes);
	// clear solution
	erase_solution(starting_sol);
	int a, b;
	double ab_dist;
	switch (variant_startingchoice)
	{
	case 'h':	// hull
		build_convex_hull(g, starting_sol);
		break;
	case 'r':	// random
		a = (int)(random() * g->nnodes);
		b = (int)(random() * g->nnodes);
		if (a == b) b = (b + 1) % g->nnodes;
		ab_dist = dist(g, a, b);
		add_edge_solution(starting_sol, a, b, ab_dist);
		starting_sol->handle_node = a;
		break;
	case 'f':	// furthest
		furthest_nodes_sol(g, starting_sol);
		break;
	default:
		print_error_ext(ERR_HEUR_CONSTR_PARAM, "variant, got  %c, should be h/r/f", variant_startingchoice);
		break;

	}
	// plot solution if needed
	if (VERBOSITY >= LOGLVL_DEBUGPLOT_P) plot_tsp_solution_undirected(g, starting_sol);

	// compute set of unvisited nodes
	SetOfNodes* set_unv_start = compute_set_unvisited_nodes(starting_sol);
	SetOfNodes* set_unv_curr = SETN_new(set_unv_start->max_size, set_unv_start->nnodes);

	do
	{
		// copy starting solution on the current solution
		deep_copy_solution(curr_sol, starting_sol);
		// copy unvisited nodes on the current set of unvisited nodes
		SETN_deepcopy(set_unv_curr, set_unv_start);

		// ******************************** FILL SOLUTION ********************************
		// for every remaining node, do an extra milage move
		int i = 0;
		int log_period = 100;
		do
		{
			i++;
			if ((VERBOSITY >= LOGLVL_DEBUGPLOT_PP) && (i % log_period == 0)) plot_tsp_solution_undirected(g, curr_sol);

		} while (extramileage_move(curr_sol, g, set_unv_curr , &cnstdata->grasp));

		// check if current solution is the best until now, if so swap them
		if (curr_sol->cost < best_sol->cost)
		{
			swap_sol = best_sol;
			best_sol = curr_sol;
			curr_sol = swap_sol;
			log_line_ext(VERBOSITY, LOGLVL_INFO, "[INFO] Found a new best solution with cost %f", best_sol->cost);
		}

		// if timelimit is reached, break
		if (second() - start_time >= timelim) break;

	} while (p_dev > 0.0);

	// free sets
	SETN_free(set_unv_start);
	SETN_free(set_unv_curr);

	// copy the best known solution into the final solution (freeing what is not needed)
	if (sol != best_sol) shallow_copy_solution(sol, best_sol);
	else // free the other solution
	{
		free(curr_sol->succ); free(curr_sol->xstar);
	}

	// free starting solution
	free(starting_sol->succ); free(starting_sol->xstar);

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
	Instance* inst = optdata->inst;
	CPXENVptr env = optdata->cpx->env;
	CPXLPptr lp = optdata->cpx->lp;

	// unpack rfndata
	HardfixingData* rfn_data = (HardfixingData*)data;
	char variant = rfn_data->variant;
	double* values = rfn_data->values;
	int progress = rfn_data->progress;
	int num_values = rfn_data->num_values;

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
	cpx_timelimit(optdata, timelim);
	cpx_warmstart(optdata, xstar);
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
	cpx_extract_sol_obj(optdata, sol, "Hardfixing");
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
				if (rfn_data->progress == num_values - 1) sol->optimal = 1;
				// increase progress in schedule (capped)
				rfn_data->progress = min(progress + 1, num_values - 1);
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
	Instance* inst = optdata->inst;
	CPXENVptr env = optdata->cpx->env;
	CPXLPptr lp = optdata->cpx->lp;

	// unpack rfndata
	LocalbranchingData* rfn_data = (LocalbranchingData*)data;
	double* values = rfn_data->values;
	int progress = rfn_data->progress;
	int num_values = rfn_data->num_values;

	// extract data
	int nnodes = sol->nnodes;

	// extract current solution
	double* xstar = sol->xstar;
	int ncols = inst->inst_model.ncols;

	// extract number k of edges to replace
	int k_replace = (int)values[progress];
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
	cpx_add_cut(env, lp, nnz, rhs, sense, index, value, CUT_STATIC, "k-OPT neighborhood constraint", -1);

	// ******************************** OPTIMIZATION *********************************
	// optimize with cplex
	int error;
	cpx_timelimit(optdata, timelim);
	cpx_warmstart(optdata, xstar);
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
	cpx_extract_sol_obj(optdata, sol, "Localbranching");
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
			if (rfn_data->progress == num_values - 1) sol->optimal = 1;
			// increase progress in schedule (capped)
			rfn_data->progress = min(progress + 1, num_values - 1);
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

/* **************************************************************************************************
*					REFINING HEURISTIC: MIX OF HARDFIXING AND LOCAL BRANCHING
*			- schedule whether to use hardfixing (first phases) or local branching (in the end)
************************************************************************************************** */
void refine_mixhardsoft(OptData* optdata, Solution* sol, void* data, double timelim)
{

	// unpack rfndata
	MixHardSoftData* rfn_data = (MixHardSoftData*)data;
	int progress = rfn_data->progress;
	//printf("%d,%d,%d,%d\n", rfn_data->num_values, rfn_data->rounds_nonimprove, rfn_data->max_rounds, rfn_data->progress);

	if (progress < HEUR_MIXHARDSOFT_VALUES / 2.0)
	{
		rfn_data->variant = 's';
		refine_hardfixing(optdata, sol, data, timelim);
	}
	else
	{
		refine_localbranching(optdata, sol, data, timelim);
	}
}

/* **************************************************************************************************
*					REFINING HEURISTIC: 2-OPT moves
*			- make 2-OPT moves until a local minima or timelimit is reached
************************************************************************************************** */
void refine_2opt(OptData* optdata, Solution* sol, void* data, double timelim)
{
	// start time
	double start_time = second();

	Opt2MoveData* rfn_data = (Opt2MoveData*)data;
	TabuList* tabu = rfn_data->tabu;
	double delta;

	do
	{
		// compute move
		delta = move_2opt(sol->succ, &optdata->inst->inst_graph, tabu, rfn_data->allow_worsening);
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Solution cost changed by: %f through a 2-OPT move", delta);
		// change cost
		sol->cost += delta;

		if (using_tabu(tabu) && !(delta >= 0 && !rfn_data->allow_worsening)) TABU_advance(tabu);

		// allow 1 worsening movement if it was already allowed
		if (delta > 0)
		{
			log_line(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Allowed a worsening move");
			delta = -1;
		}

		rfn_data->allow_worsening = 0;

	} while ((second() - start_time < timelim) && (delta < 0));

	// if plateau and worsening is now allowed, solution is local optimum
	if (delta == 0 && !rfn_data->allow_worsening)
	{
		// mark solution as local optimum
		sol->optimal = 1;
		log_line(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Local minimum in 2-OPT neighborhood reached");	
	}

}