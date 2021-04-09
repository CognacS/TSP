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
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, &ptr_cname)) print_error("wrong CPXnewcols on x var.s", ERR_CPLEX);
			// check correctness of xpos
			if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, g->nnodes)) print_error("wrong position for x var.s using function \"xpos\"", ERR_INCORRECT_FUNCTION);
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
	instance* inst = (instance*)userhandle;
	model* m = &inst->inst_model;

	// get current solution and objective value
	double* xstar;	arr_malloc_s(xstar, m->ncols, double);
	double objval = CPX_INFBOUND;
	if (CPXcallbackgetcandidatepoint(context, xstar, 0, m->ncols - 1, &objval)) print_error("CPXcallbackgetcandidatepoint error", ERR_CPLEX);

	// add SEC on callback
	add_sec_on_subtours(context, NULL, inst, xstar, -1, CUT_CALLBACK_REJECT);

	// CLEANUP
	free(xstar);
	return 0;
	
}

/* **************************************************************************************************
*				ADD SUBTOUR ELIMINATION CONSTRAINTS ON AN INFEASABLE SOLUTION
************************************************************************************************** */
int add_sec_on_subtours(void* env, void* cbdata, instance* inst, double* xstar, int wherefrom, int purgeable)
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
			if (h >= g->nnodes) print_error("needed more components but found no node", ERR_GENERIC_INCONSIST);

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
		// solve model and get solution xstar
		if (error = CPXmipopt(env, lp))
		{
			printf("CPX error %d\n", error);
			print_error("CPXmipopt() Benders", ERR_CPLEX);
		}
		if (error = CPXgetx(env, lp, xstar, 0, ncols - 1))
		{
			printf("CPX error %d\n", error);
			print_error("CPXgetx() Benders", ERR_CPLEX);
		}

		// produce new cuts on violated constraints
		newcuts = add_sec_on_subtours(env, lp, inst, xstar, -1, CUT_STATIC);
		log_line_ext(VERBOSITY, LOGLVL_INFO, "Added %d new SEC's", newcuts);

	} while (newcuts);

	// CLEANUP
	free(xstar);

}

/* **************************************************************************************************
*						SOLUTION USING CALLBACK'S METHOD
************************************************************************************************** */
void solve_callback(instance* inst, CPXENVptr env, CPXLPptr lp)
{
	int error;

	// build naive model
	build_model_base_undirected(inst, env, lp);

	// add lazy constraints when there is a new candidate
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
	if (error = CPXcallbacksetfunc(env, lp, contextid, sec_callback, inst))
	{
		printf("CPX error %d\n", error);
		print_error("CPXcallbacksetfunc() error", ERR_CPLEX);
	}

	// solve the problem with the callback
	if (error = CPXmipopt(env, lp))
	{
		printf("CPX error %d\n", error);
		print_error("CPXmipopt() Callback", ERR_CPLEX);
	}

}

/* **************************************************************************************************
*						SOLUTION USING UNDIRECTED GRAPHS MODELS
************************************************************************************************** */
void solve_symmetric_tsp(instance* inst, CPXENVptr env, CPXLPptr lp)
{
	modeltype mt = inst->inst_params.model_type;
	// if the tsptype of the model is not asymmetric, throw error
	if (model_tsptype(mt) != MODEL_TSP_SYMM) print_error("", ERR_WRONG_TSP_PROCEDURE);

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
		print_error("symmetric variant", ERR_MODEL_NOT_IMPL);
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
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, &ptr_cname)) print_error("wrong CPXnewcols on x var.s", ERR_CPLEX);
			// check correctness of xxpos
			if (CPXgetnumcols(env, lp) - 1 != xxpos(i, j, g->nnodes)) print_error("wrong position for x var.s using function \"xxpos\"", ERR_INCORRECT_FUNCTION);

			// 2 - define x(j,i)
			// define name of variable
			sprintf(cname, "x(%d,%d)", j + 1, i + 1);
			// add variable in CPX
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, &ptr_cname)) print_error("wrong CPXnewcols on x var.s", ERR_CPLEX);
			// check correctness of xxpos
			if (CPXgetnumcols(env, lp) - 1 != xxpos(j, i, g->nnodes)) print_error("wrong position for x var.s using function \"xxpos\"", ERR_INCORRECT_FUNCTION);
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
		if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &integer, &ptr_cname)) print_error("wrong CPXnewcols on u var.s", ERR_CPLEX);
		// check correctness of upos
		if (CPXgetnumcols(env, lp) - 1 != upos(i, g->nnodes)) print_error("wrong position for u var.s using function \"upos\"", ERR_INCORRECT_FUNCTION);
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
			if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &integer, &ptr_cname)) print_error("wrong CPXnewcols on x var.s", ERR_CPLEX);
			// check correctness of xxpos
			if (CPXgetnumcols(env, lp) - 1 != ypos(i, j, g->nnodes)) print_error("wrong position for x var.s using function \"ypos\"", ERR_INCORRECT_FUNCTION);

			// 2 - define y(j,i)
			if (i == 0) ub = 0;
			// define name of variable
			sprintf(cname, "y(%d,%d)", j + 1, i + 1);
			// add variable in CPX
			if (CPXnewcols(env, lp, 1, NULL, &lb, &ub, &integer, &ptr_cname)) print_error("wrong CPXnewcols on x var.s", ERR_CPLEX);
			// check correctness of xxpos
			if (CPXgetnumcols(env, lp) - 1 != ypos(j, i, g->nnodes)) print_error("wrong position for x var.s using function \"ypos\"", ERR_INCORRECT_FUNCTION);
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
	if (model_tsptype(mt) != MODEL_TSP_ASYMM) print_error("", ERR_WRONG_TSP_PROCEDURE);

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
		print_error("asymmetric variant", ERR_MODEL_NOT_IMPL);
	}

	// add SEC on pairs if needed
	if (need_sec(mt)) add_sec2_asymmetric(inst, env, lp);

	// solve the model
	int error = 0;
	if (error = CPXmipopt(env, lp))
	{
		printf("CPX error %d\n", error);
		print_error("CPXmipopt() Asymmetric", ERR_CPLEX);
	}

}
