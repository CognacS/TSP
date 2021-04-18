#include "../cpx_models.h"

/* **************************************************************************************************
*						BASE MODEL FOR UNDIRECTED GRAPHS
************************************************************************************************** */
void build_model_base_undirected(instance* inst, CPXENVptr env, CPXLPptr lp)
{

	// ********************************* SETUP *********************************
	// extract values
	graph* g = &inst->inst_graph;
	model* m = &inst->inst_model;
	// define constants
	double zero = 0.0;
	char binary = 'B';
	char cname[100];
	char* ptr_cname = cname;
	// define bounds of binary variables
	double lb = 0.0;
	double ub = 1.0;

	// ************************ ADD COLUMNS (VARIABLES) ************************
	// add binary var.s x(i,j) for i < j  
	for (int i = 0; i < g->nnodes; i++)
	{
		for (int j = i + 1; j < g->nnodes; j++)
		{
			// define name of variable
			sprintf(cname, "x(%d,%d)", i + 1, j + 1);
			// define coefficients
			double obj = dist(i, j, g); // cost == distance 
			// add variable in CPX
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, &ptr_cname))
				print_error(ERR_CPLEX, "wrong CPXnewcols on x var.s");
			// check correctness of xpos
			if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, g->nnodes))
				print_error(ERR_INCORRECT_FUNCTION_IMPL, "wrong position for x var.s using function \"xpos\"");
		}
	}
	// update number of columns
	m->ncols = CPXgetnumcols(env, lp);

	// ************************ ADD ROWS (CONSTRAINTS) ************************
	// set degree of contraint
	double rhs = 2.0;
	// 'E' for equality constraint 
	char sense = 'E';
	// prepare arrays for constraints
	int nnz = 0;
	int* index;		arr_malloc_s(index, g->nnodes - 1, int);
	double* value;	arr_malloc_s(value, g->nnodes - 1, double);

	// add the degree constraints
	for (int h = 0; h < g->nnodes; h++)  // degree constraints
	{
		// define name of constraint
		sprintf(cname, "degree(%d)", h + 1);

		nnz = 0;
		// fill the added constraint
		for (int i = 0; i < g->nnodes; i++)
		{
			// no constraint with itself
			if (i == h) continue;
			// constraint with node i
			index[nnz] = xpos(i, h, g->nnodes);
			value[nnz] = 1.0;
			nnz++;
		}

		mip_add_cut(env, lp, -1, nnz, rhs, sense, index, value, cname, CUT_STATIC);
	}
	// CLEANUP
	free(index);
	free(value);

}

/* **************************************************************************************************
*				ADD SUBTOUR ELIMINATION CONSTRAINTS ON AN INFEASABLE SOLUTION
************************************************************************************************** */
static int CPXPUBLIC sec_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle)
{
	// get instance handle
	callback_instance* cb_inst = (callback_instance*)userhandle;
	instance* inst = cb_inst->inst;
	model* m = &inst->inst_model;

	// get stats on context
	int cpx_node = -1;
	CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &cpx_node);
	int cpx_depth = -1;

	// get current solution and objective value
	double* xstar = NULL;
	double objval = CPX_INFBOUND;

	// define some values
	double rval;

	// apply separation or reject based on contextid
	switch (contextid)
	{
	case CPX_CALLBACKCONTEXT_CANDIDATE:
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "Entered callback candidate at node %d", cpx_node);

		// allocate xstar
		arr_malloc_s(xstar, m->ncols, double);
		// get candidate solution
		if (CPXcallbackgetcandidatepoint(context, xstar, 0, m->ncols - 1, &objval))
			print_error(ERR_CPLEX, "CPXcallbackgetcandidatepoint error");
		// execute rejection procedure
		if (cb_inst->rej_procedure)
			cb_inst->rej_procedure(context, NULL, inst, cb_inst->args, xstar, -1, CUT_CALLBACK_REJECT, 0);
		else
			print_error(ERR_CB_UNDEF_PROCEDURE, "candidate");
		// CLEANUP
		free(xstar);
		break;

	case CPX_CALLBACKCONTEXT_RELAXATION:
		
		// randomly choose whether to apply separation or not
		CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODEDEPTH, &cpx_depth);
		safe_rand(&rval);
		double decay_h = 0.3;
		double decay_e = 0.1;
		double decay_s = 0.1;
		double decay_r = 0.3;
		double hyperbolic_decay = 1.0 / (decay_h * cpx_depth + 1);
		double exponential_decay = exp(-decay_e * cpx_depth);
		double sigmoid_decay = 2.0 / (exp(decay_s * cpx_depth) + 1);
		double relu_decay = max(-0.1 * decay_r * cpx_depth + 1, 0);
		printf("%d -> %f, %f, %f, %f\n", cpx_depth, hyperbolic_decay, exponential_decay, sigmoid_decay, relu_decay);
		if (rval < (1 - relu_decay)) break;

		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "Entered callback relaxation at node %d", cpx_node);

		// allocate xstar
		arr_malloc_s(xstar, m->ncols, double);
		// get candidate solution
		if (CPXcallbackgetrelaxationpoint(context, xstar, 0, m->ncols - 1, &objval))
			print_error(ERR_CB_UNDEF_PROCEDURE, "CPXcallbackgetrelaxationpoint error");
		// execute separation procedure
		if (cb_inst->sep_procedure)
			cb_inst->sep_procedure(context, NULL, inst, cb_inst->args, xstar, -1, CPX_USECUT_FILTER, CUT_FLAGS_MINCUT_ENABLED);
		else
			print_error(ERR_CB_UNDEF_PROCEDURE, "candidate");
		// CLEANUP
		free(xstar);
		break;
		
	}
	
	return 0;
	
}

/* **************************************************************************************************
*				ADD SUBTOUR ELIMINATION CONSTRAINTS ON AN INFEASABLE SOLUTION
************************************************************************************************** */
int add_sec_on_subtours(void* env, void* cbdata, instance* inst, void* args, double* xstar, int wherefrom, int purgeable, int flags)
{
	// extract structures
	graph* g = &inst->inst_graph;
	model* m = &inst->inst_model;

	// ************************** FIND CONNECTED COMPONENTS **************************
	// allocate connected components arrays
	int ncomp = 0;
	int* succ = NULL, * comp = NULL;
	arr_malloc_s(succ, g->nnodes, int);
	arr_malloc_s(comp, g->nnodes, int);
	// find connected components
	find_conncomps_dfs(g, xstar, succ, comp, &ncomp);

	// if there is more than 1 connected component -> infeasable solution!
	// must add SEC's
	if (ncomp > 1)
	{
		// ********** ADD CUTS (1 SEC for each connected component) ***************
		int nnz = 0;
		double rhs = 0;
		char sense = 'L';
		double* value;		arr_malloc_s(value, m->ncols, double);
		int* index;			arr_malloc_s(index, m->ncols, int);
		// flag array for indicating a visited component
		char* visitedcomp;	calloc_s(visitedcomp, ncomp, char);
		// for each connected component add a constraint
		int h = 0;
		for (int c = 0; c < ncomp; c++)
		{
			// find a node from a non-visited component
			while (h < g->nnodes && visitedcomp[comp[h]]) h++;
			if (h >= g->nnodes)
				print_error(ERR_MODEL_INCONSISTENT, "needed more components but found no other node");

			// setup for current component
			rhs = 0;
			nnz = 0;
			// elaborate the component starting from node h
			// i is the current arrival point
			for (int i = succ[h]; i != h; i = succ[i])
			{
				rhs++;
				// j starts from the beginning of the subtour
				// and iterates all previous nodes with respect to i
				// and add the coefficients wrt x_ij
				for (int j = h; j != i; j = succ[j])
				{
					index[nnz] = xpos(i, j, g->nnodes);
					value[nnz] = 1.0;
					nnz++;
				}
			}
			// flag the component as visited and go on with nodes
			visitedcomp[comp[h++]] = 1;

			// add SEC
			mip_add_cut(env, cbdata, wherefrom, nnz, rhs, sense, index, value, "SEC", purgeable);
		}

		// CLEAN UP SEC ARRAYS
		free(index);
		free(value);
		free(visitedcomp);
	}

	// CLEAN UP CONN COMP ARRAYS
	free(comp);
	free(succ);

	// return how many SEC's where added
	return ncomp > 1? ncomp : 0;

}

/* **************************************************************************************************
*				ADD SUBTOUR ELIMINATION CONSTRAINTS USING CONCORDE
************************************************************************************************** */
int CC_add_sec_on_subtours(void* env, void* cbdata, instance* inst, void* args, double* xstar, int wherefrom, int purgeable, int flags)
{
	// extract structures
	graph* g = &inst->inst_graph;
	model* m = &inst->inst_model;

	// make input values
	int nnodes = g->nnodes;
	int nedges = nnodes * (nnodes - 1) / 2;
	int* elist = (int*)args;
	// make output values
	int ncomp;
	int* compscount;
	int* comps;

	// search for connected components using the Concorde method
	if (CCcut_connect_components(nnodes, nedges, elist, xstar, &ncomp, &compscount, &comps))
		print_error(ERR_CONCORDE, "CCcut_connect_components()");

	// if there are connected components, add the SEC
	if (ncomp > 1)
	{

		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "Adding SECs");
		// ********** ADD CUTS (1 SEC for each connected component) ***************
		int nnz = 0;
		double rhs = 0;
		char sense = 'L';
		double* value;		arr_malloc_s(value, m->ncols, double);
		int* index;			arr_malloc_s(index, m->ncols, int);
		// for each connected component add a constraint
		int comps_pos = 0;	// position in array comps to keep track of the conncomp
		for (int c = 0; c < ncomp; c++)
		{
			// setup for current component
			rhs = compscount[c] - 1;
			nnz = 0;
			// elaborate the component (compscount[c] * (compscount[c] - 1) / 2 edges!)
			for (int i = comps_pos; i < comps_pos + compscount[c]; i++)
			{
				for (int j = i + 1; j < comps_pos + compscount[c]; j++)
				{
					index[nnz] = xpos(comps[i], comps[j], g->nnodes);
					value[nnz] = 1.0;
					nnz++;
				}
			}
			// increase the position in comps
			comps_pos += compscount[c];

			// add SEC
			mip_add_cut(env, cbdata, wherefrom, nnz, rhs, sense, index, value, "SEC", purgeable);
		}

		free(value);
		free(index);

	}
	
	// if there is just one connected component, add mincut constraints
	if ((flags & CUT_FLAGS_MINCUT_ENABLED) && ncomp == 1)
	{
		// make concorde callback structure
		concorde_instance cc_inst;
		cc_inst.inst = inst;
		cc_inst.context = (CPXCALLBACKCONTEXTptr)env;
		arr_malloc_s(cc_inst.value, m->ncols, double);
		arr_malloc_s(cc_inst.index, m->ncols, int);
		calloc_s(cc_inst.labels, g->nnodes, int);

		// produce cuts
		log_line_ext(VERBOSITY, LOGLVL_DEBUG, "Adding mincut constraints");
		if (CCcut_violated_cuts(nnodes, nedges, elist, xstar, CC_CUTOFF, doit_fn_concorde2cplex, &cc_inst))
			print_error(ERR_CONCORDE, "CCcut_violated_cuts()");

		free(cc_inst.value);
		free(cc_inst.index);
		free(cc_inst.labels);
	}

	// CLEANUP
	CC_IFFREE(compscount, int);
	CC_IFFREE(comps, int);

	return ncomp > 1 ? ncomp : 0;
}


/* **************************************************************************************************
*						CONCORDE CALLBACK FUNCTION TO ADD A CUT TO CPLEX
************************************************************************************************** */
int doit_fn_concorde2cplex(double cutval, int cutcount, int* cut, void* args)
{
	concorde_instance* cc_inst = (concorde_instance*)args;
	instance* inst = cc_inst->inst;
	graph* g = &inst->inst_graph;
	model* m = &inst->inst_model;

	if (cutcount > 2 && cutcount < g->nnodes)
	{
		char sense = 'G';
		double rhs = 2.0;
		int nnz = 0;
		double* value = cc_inst->value;
		int* index = cc_inst->index;
		int* labels = cc_inst->labels;

		// label cut nodes
		for (int i = 0; i < cutcount; i++)
		{
			labels[cut[i]] = 1;
		}

		// make the indices for each edge in the cut
		for (int i = 0; i < cutcount; i++)
		{
			for (int j = 0; j < g->nnodes; j++)
			{
				// if the node j is not inside S
				if (labels[j] == 0)
				{
					index[nnz] = xpos(cut[i], j, g->nnodes);
					value[nnz] = 1.0;
					nnz++;
				}
			}
		}

		// add cut
		mip_add_cut(cc_inst->context, NULL, -1, nnz, rhs, sense, index, value, NULL, CPX_USECUT_FILTER);

		// reset labels for the next call
		for (int i = 0; i < cutcount; i++)
		{
			labels[cut[i]] = 0;
		}
	}

	return 0;
}

/* **************************************************************************************************
*						SOLUTION USING BENDERS' METHOD
************************************************************************************************** */
void solve_benders(instance* inst, CPXENVptr env, CPXLPptr lp)
{
	int error;
	// extract structures
	graph* g = &inst->inst_graph;

	// build naive model
	build_model_base_undirected(inst, env, lp);
	int ncols = CPXgetnumcols(env, lp);

	// build solution array
	double* xstar = NULL;
	calloc_s(xstar, ncols, double);

	// number of cuts added at each iteration
	int newcuts = 0;

	do
	{
		// solve model
		if (error = CPXmipopt(env, lp))
			print_error_ext(ERR_CPLEX, "CPXmipopt() Benders, CPX error: %d", error);

		// get solution xstar
		if (error = CPXgetx(env, lp, xstar, 0, ncols - 1))
			print_error_ext(ERR_CPLEX, "CPXgetx() Benders, CPX error: %d", error);

		// produce new cuts on violated constraints
		newcuts = add_sec_on_subtours(env, lp, inst, NULL, xstar, -1, CUT_STATIC, 0);
		log_line_ext(VERBOSITY, LOGLVL_INFO, "Added %d new SEC's", newcuts);

		// update time limit
		mip_timelimit(env, inst->inst_params.timelimit, inst);

	} while (newcuts && !time_limit_expired(inst));

	// CLEANUP
	free(xstar);

}

/* **************************************************************************************************
*						SOLUTION USING CALLBACK'S METHOD
************************************************************************************************** */
void solve_callback(instance* inst, CPXENVptr env, CPXLPptr lp)
{
	int error;
	modeltype mt = inst->inst_params.model_type;
	int nnodes = inst->inst_graph.nnodes;

	// ********************** SETUP STARTING MODEL **********************
	// build naive model
	build_model_base_undirected(inst, env, lp);

	// ********************** SETUP CALLBACK STRUCTURES **********************
	// make flags for model choice
	int use_cc_sep = use_cc_on_sep(mt);
	int use_cc_rej = use_cc_on_rej(mt);

	void* args = NULL;

	// if separation or rejection use CC method, then construct elist to be shared
	// as argument
	if (use_cc_sep || use_cc_rej)
	{
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
		args = (void*)elist;
	}

	// make callback instance
	callback_instance cb_inst;
	cb_inst.inst = inst;
	cb_inst.args = args;
	cb_inst.sep_procedure = use_cc_sep ? CC_add_sec_on_subtours : NULL;
	cb_inst.rej_procedure = use_cc_rej ? CC_add_sec_on_subtours : add_sec_on_subtours;

	// ********************** SETUP CALLBACK FUNCTION **********************
	// add lazy constraints when there is a new candidate
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
	// choose to use or not the additional separation
	if (use_cc_sep) contextid = contextid | CPX_CALLBACKCONTEXT_RELAXATION;

	// set the callback function
	if (error = CPXcallbacksetfunc(env, lp, contextid, sec_callback, &cb_inst))
		print_error_ext(ERR_CPLEX, "CPXcallbacksetfunc() Callback, CPX error: %d", error);

	// ********************** OPTIMIZE **********************
	// solve the problem with the callback
	if (error = CPXmipopt(env, lp))
		print_error_ext(ERR_CPLEX, "CPXmipopt() Callback, CPX error: %d", error);


	// CLEANUP
	free(args);

}

/* **************************************************************************************************
*						SOLUTION USING UNDIRECTED GRAPHS MODELS
************************************************************************************************** */
void solve_symmetric_tsp(instance* inst, CPXENVptr env, CPXLPptr lp)
{
	modeltype mt = inst->inst_params.model_type;
	// if the tsptype of the model is not asymmetric, throw error
	if (model_tsptype(mt) != MODEL_TSP_SYMM)
		print_error(ERR_WRONG_TSP_PROCEDURE, NULL);

	// build symmetric TSP model OR use symmetric TSP procedure
	switch (model_archetype(mt))
	{
	case MODEL_SY_BEND:
		solve_benders(inst, env, lp);
		break;
	case MODEL_SY_CLBCK:
		solve_callback(inst, env, lp);
		break;
	default:
		print_error(ERR_MODEL_NOT_IMPL, "symmetric variant");
	}
}


/* **************************************************************************************************
*						BASE MODEL FOR DIRECTED GRAPHS
************************************************************************************************** */
void build_model_base_directed(instance* inst, CPXENVptr env, CPXLPptr lp)
{

	// ********************************* SETUP *********************************
	// extract values
	graph* g = &inst->inst_graph;
	model* m = &inst->inst_model;
	// define constants
	double zero = 0.0;
	char binary = 'B';
	char cname[100];
	char* ptr_cname = cname;
	
	// define bounds of binary variables
	double lb = 0.0;
	double ub = 1.0;

	// ************************ ADD COLUMNS (VARIABLES) ************************
	// add binary var.s x(i,j) for all i=/=j
	for (int i = 0; i < g->nnodes; i++)
	{
		for (int j = i + 1; j < g->nnodes; j++)
		{
			// define cost
			double obj = dist(i, j, g); // cost == distance

			// 1 - define x(i,j)
			// define name of variable
			sprintf(cname, "x(%d,%d)", i + 1, j + 1);
			// add variable in CPX
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, &ptr_cname))
				print_error(ERR_CPLEX, "wrong CPXnewcols on x var.s");
			// check correctness of xxpos
			if (CPXgetnumcols(env, lp) - 1 != xxpos(i, j, g->nnodes))
				print_error(ERR_INCORRECT_FUNCTION_IMPL, "wrong position for x var.s using function \"xxpos\"");
			// 2 - define x(j,i)
			// define name of variable
			sprintf(cname, "x(%d,%d)", j + 1, i + 1);
			// add variable in CPX
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, &ptr_cname))
				print_error(ERR_CPLEX, "wrong CPXnewcols on x var.s");
			// check correctness of xxpos
			if (CPXgetnumcols(env, lp) - 1 != xxpos(j, i, g->nnodes))
				print_error(ERR_INCORRECT_FUNCTION_IMPL, "wrong position for x var.s using function \"xxpos\"");
		}
	}
	// update number of columns
	m->ncols = CPXgetnumcols(env, lp);

	// ************************ ADD ROWS (CONSTRAINTS) ************************
	// set degree of contraint
	double rhs = 1.0;
	// 'E' for equality constraint 
	char sense = 'E';
	// prepare arrays for constraints
	int nnz = 0;
	int*	index_out;	arr_malloc_s(index_out, g->nnodes - 1, int);
	double* value_out;	arr_malloc_s(value_out, g->nnodes - 1, double);
	int*	index_in;	arr_malloc_s(index_in, g->nnodes - 1, int);
	double* value_in;	arr_malloc_s(value_in, g->nnodes - 1, double);
	// add the degree constraints
	for (int h = 0; h < g->nnodes; h++)  // degree constraints
	{
		// reset nnz
		nnz = 0;
		// fill the added constraints
		for (int i = 0; i < g->nnodes; i++)
		{
			// no constraint with itself
			if (i == h) continue;
			// constraint with node i
			index_out[nnz] = xxpos(h, i, g->nnodes);
			index_in[nnz]  = xxpos(i, h, g->nnodes);
			value_out[nnz] = 1.0;
			value_in[nnz] =  1.0;
			nnz++;
		}
		// OUTGOING CONSTRAINT
		// define name of constraints
		sprintf(cname, "degree_out(%d)", h + 1);
		mip_add_cut(env, lp, -1, nnz, rhs, sense, index_out, value_out, cname, CUT_STATIC);
		// INGOING CONSTRAINT
		// define name of constraint
		sprintf(cname, "degree_in(%d)", h + 1);
		mip_add_cut(env, lp, -1, nnz, rhs, sense, index_in, value_in, cname, CUT_STATIC);
	}

	// CLEANUP
	free(index_out);	free(index_in);
	free(value_out);	free(value_in);
}


/* **************************************************************************************************
*			MILLER, TUCKER AND ZEMLIN COMPACT MODEL (MTZ)
************************************************************************************************** */
void build_model_mtz(instance* inst, CPXENVptr env, CPXLPptr lp)
{
	// construct the base model for directed graphs
	build_model_base_directed(inst, env, lp);

	// extract values
	graph* g = &inst->inst_graph;
	modeltype mt = inst->inst_params.model_type;
	// define name of constraints
	char cname[100];

	// define bounds of integer values
	double lb = 0.0;
	double ub = g->nnodes - 2;
	char integer = 'I';
	char* ptr_cname = cname;

	// ************************ ADD COLUMNS (VARIABLES) ************************
	// add integer var.s u(i) for all i but the first node
	for (int i = 1; i < g->nnodes; i++)
	{
		// define name of variable
		sprintf(cname, "u(%d)", i + 1);
		// add variable in CPX
		if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &integer, &ptr_cname))
			print_error(ERR_CPLEX, "wrong CPXnewcols on u var.s");
		// check correctness of upos
		if (CPXgetnumcols(env, lp) - 1 != upos(i, g->nnodes))
			print_error(ERR_INCORRECT_FUNCTION_IMPL, "wrong position for u var.s using function \"upos\"");
	}

	// ************************ ADD ROWS (CONSTRAINTS) ************************
	// add static/lazy constraints  1.0 * u_i - 1.0 * u_j + M * x_ij <= M - 1, for each arc (i,j) not touching node 0
	// get type to constraint to add to the model
	int constr_type = variant2constr(mt);
	// define arrays
	int index[3];
	double value[3];
	// define constants
	double big_M = g->nnodes - 1.0;
	double rhs = big_M - 1.0;
	char sense = 'L';
	int nnz = 3;
	for (int i = 1; i < g->nnodes; i++) // excluding node 0
	{
		for (int j = 1; j < g->nnodes; j++) // excluding node 0 
		{
			if (i == j) continue;
			sprintf(cname, "u-consistency for arc (%d,%d)", i + 1, j + 1);
			index[0] = upos(i, g->nnodes);
			value[0] = 1.0;
			index[1] = upos(j, g->nnodes);
			value[1] = -1.0;
			index[2] = xxpos(i, j, g->nnodes);
			value[2] = big_M;

			// add constraint
			mip_add_cut(env, lp, -1, nnz, rhs, sense, index, value, cname, constr_type);
		}
	}
}

/* **************************************************************************************************
*			GAVISH AND GRAVES SINGLE COMMODITY FLOW MODEL (GG)
************************************************************************************************** */
void build_model_gg(instance* inst, CPXENVptr env, CPXLPptr lp)
{
	// construct the base model for directed graphs
	build_model_base_directed(inst, env, lp);

	// extract values
	graph* g = &inst->inst_graph;
	modeltype mt = inst->inst_params.model_type;

	char cname[100];

	char integer = 'I';
	char* ptr_cname = cname;

	// ************************ ADD COLUMNS (VARIABLES) ************************
	// add integer var.s 0<=y(i, j)<=n-2 for all i=/=j (0<=y(i, 1)<=0)
	for (int i = 0; i < g->nnodes; i++)
	{
		for (int j = i + 1; j < g->nnodes; j++)
		{
			// define bounds of integer values
			double lb = 0.0;
			double ub = g->nnodes - 1.0;

			// 1 - define y(i,j)
			// define name of variable
			sprintf(cname, "y(%d,%d)", i + 1, j + 1);
			// add variable in CPX
			if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &integer, &ptr_cname))
				print_error(ERR_CPLEX, "wrong CPXnewcols on x var.s");
			// check correctness of xxpos
			if (CPXgetnumcols(env, lp) - 1 != ypos(i, j, g->nnodes))
				print_error(ERR_INCORRECT_FUNCTION_IMPL, "wrong position for x var.s using function \"ypos\"");

			// 2 - define y(j,i)
			if (i == 0) ub = 0;
			// define name of variable
			sprintf(cname, "y(%d,%d)", j + 1, i + 1);
			// add variable in CPX
			if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &integer, &ptr_cname))
				print_error(ERR_CPLEX, "wrong CPXnewcols on x var.s");
			// check correctness of xxpos
			if (CPXgetnumcols(env, lp) - 1 != ypos(j, i, g->nnodes))
				print_error(ERR_INCORRECT_FUNCTION_IMPL, "wrong position for x var.s using function \"ypos\"");
		}
	}

	// ************************ ADD ROWS (CONSTRAINTS) ************************
	int constr_type = variant2constr(mt);
	double rhs = 0.0;
	char sense = 'E';
	int index[2];
	double value[2];
	int nnz = 2;
	// add static/lazy constraints (n-1) * x_1j - 1.0 * y_1j == 0 for each j but the first node
	for (int j = 1; j < g->nnodes; j++)
	{
		sprintf(cname, "commodity quantity on arc (%d,%d)", 1, j + 1);
		index[0] = xxpos(0, j, g->nnodes);
		value[0] = g->nnodes - 1.0;
		index[1] = ypos(0, j, g->nnodes);
		value[1] = -1.0;
		mip_add_cut(env, lp, -1, nnz, rhs, sense, index, value, cname, constr_type);
	}

	sense = 'G';
	// add static/lazy constrains (n-2) * x_ij - 1.0 * y_ij >= 0.0, for each arc (i,j) with i,j=/=0
	for (int i = 1; i < g->nnodes; i++)
	{
		for (int j = 1; j < g->nnodes; j++)
		{
			if (i == j) continue;

			sprintf(cname, "commodity quantity on arc (%d,%d)", i + 1, j + 1);
			index[0] = xxpos(i, j, g->nnodes);
			value[0] = g->nnodes - 2.0;
			index[1] = ypos(i, j, g->nnodes);
			value[1] = -1.0;
			mip_add_cut(env, lp, -1, nnz, rhs, sense, index, value, cname, constr_type);
		}
	}
	
	// allocate arrays for holding indices and values for flow constraints
	int* index_flow;		arr_malloc_s(index_flow, 2 * (g->nnodes - 1), int);
	double* value_flow;		arr_malloc_s(value_flow, 2 * (g->nnodes - 1), double);

	sense = 'E';
	rhs = 1.0;
	nnz = 2 * (g->nnodes - 1);

	// add static/lazy constrains on flow through nodes
	for (int h = 1; h < g->nnodes; h++)
	{
		// define name of constraint
		sprintf(cname, "flow(%d)", h + 1);

		int arr_idx = 0;
		// fill the added constraints
		for (int i = 0; i < g->nnodes; i++)
		{
			// no constraint with itself
			if (i == h) continue;

			// add node i to the constraint
			index_flow[arr_idx] = ypos(i, h, g->nnodes);
			value_flow[arr_idx] = 1.0;
			index_flow[arr_idx+g->nnodes-1] = ypos(h, i, g->nnodes);
			value_flow[arr_idx+g->nnodes-1] = -1.0;
			arr_idx++;
		}
		mip_add_cut(env, lp, -1, nnz, rhs, sense, index_flow, value_flow, cname, constr_type);
	}

	free(index_flow);
	free(value_flow);

}

/* **************************************************************************************************
*				ADD SUBTOUR ELIMINATION CONSTRAINT FOR PAIRS IN ASYMMETRIC TSP
************************************************************************************************** */

void add_sec2_asymmetric(instance* inst, CPXENVptr env, CPXLPptr lp)
{
	graph* g = &inst->inst_graph;

	// add lazy constraints 1.0 * x_ij + 1.0 * x_ji <= 1, for each arc (i,j) with i < j
	char cname[100];
	double rhs = 1.0;
	char sense = 'L';
	int nnz = 2;
	int index[3];
	double value[3];
	for (int i = 0; i < g->nnodes; i++)
	{
		for (int j = i + 1; j < g->nnodes; j++)
		{
			sprintf(cname, "SEC on node pair (%d,%d)", i + 1, j + 1);
			index[0] = xxpos(i, j, g->nnodes);
			value[0] = 1.0;
			index[1] = xxpos(j, i, g->nnodes);
			value[1] = 1.0;
			mip_add_cut(env, lp, -1, nnz, rhs, sense, index, value, cname, CUT_LAZY);
		}
	}
}

/* **************************************************************************************************
*						SOLUTION USING DIRECTED GRAPHS MODELS
************************************************************************************************** */
void solve_asymmetric_tsp(instance* inst, CPXENVptr env, CPXLPptr lp)
{
	modeltype mt = inst->inst_params.model_type;
	// if the tsptype of the model is not asymmetric, throw error
	if (model_tsptype(mt) != MODEL_TSP_ASYMM) print_error(ERR_WRONG_TSP_PROCEDURE, NULL);

	// build asymmetric TSP model
	switch (model_archetype(mt))
	{
	case MODEL_AS_MTZ:
		build_model_mtz(inst, env, lp);
		break;
	case MODEL_AS_GG:
		build_model_gg(inst, env, lp);
		break;
	default:
		print_error(ERR_MODEL_NOT_IMPL, "asymmetric variant");
	}

	// add SEC on pairs if needed
	if (need_sec(mt)) add_sec2_asymmetric(inst, env, lp);

	// solve the model
	int error = 0;
	if (error = CPXmipopt(env, lp))
		print_error_ext(ERR_CPLEX, "CPXmipopt() Asymmetric, CPX error: %d", error);

}
