#include "../cpx_models.h"

/* **************************************************************************************************
*						BASE MODEL FOR UNDIRECTED GRAPHS
************************************************************************************************** */
void build_model_base_undirected(graph* g, modeltype mt, CPXENVptr env, CPXLPptr lp)
{

	// ********************************* SETUP *********************************
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

	// ************************ ADD ROWS (CONSTRAINTS) ************************
	// add the degree constraints
	for (int h = 0; h < g->nnodes; h++)  // degree constraints
	{
		// get row number
		int lastrow = CPXgetnumrows(env, lp);
		// set degree of contraint
		double rhs = 2.0;
		// 'E' for equality constraint 
		char sense = 'E';
		// define name of constraint
		sprintf(cname, "degree(%d)", h + 1);
		// add empty constraint in CPX
		// TODO: setup array of constraints: (coefficients, positions)
		// pass it inside CPXnewrows!
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, &ptr_cname)) print_error("wrong CPXnewrows [degree]", ERR_CPLEX);
		// fill the added constraint
		for (int i = 0; i < g->nnodes; i++)
		{
			// no constraint with itself
			if (i == h) continue;
			// constraint with node i
			if (CPXchgcoef(env, lp, lastrow, xpos(i, h, g->nnodes), 1.0)) print_error("wrong CPXchgcoef [degree]", ERR_CPLEX);
		}
	}
}

/* **************************************************************************************************
*						SOLUTION USING BENDERS' METHOD
************************************************************************************************** */
void solve_benders(graph* g, modeltype mt, CPXENVptr env, CPXLPptr lp)
{
	int error;

	// build naive model
	build_model_base_undirected(g, mt, env, lp);
	int ncols = CPXgetnumcols(env, lp);

	// build solution array
	double* xstar = NULL;
	calloc_s(xstar, ncols, double);
	// allocate auxiliary data structures
	int* succ = NULL, *comp = NULL;
	calloc_s(succ, g->nnodes, int);
	calloc_s(comp, g->nnodes, int);
	int ncomp = 0;

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

	// arrays for rows
	double* rhs = NULL;
	char* sense = NULL;
	char** rowname = NULL;

	// arrays for coefficients
	int* rowlist = NULL;
	int* collist = NULL;
	double* vallist = NULL;

	// auxiliary arrays
	int* compsizes = NULL;
	char* visited = NULL;

	// while the solution has subtours
	while (find_conncomps_dfs(g, xstar, succ, comp, &ncomp) > 1)
	{
		// ******************************* PREPARE ROWS ******************************
		int nrows = CPXgetnumrows(env, lp);
		// get component sizes
		int tot_coeff_num = 0;
		calloc_s(compsizes, ncomp, int);
		for (int i = 0; i < g->nnodes; i++) compsizes[comp[i]]++;
		// get total number of coefficients as the number of all the edges in all subtours
		for (int row = 0; row < ncomp; row++) tot_coeff_num += compsizes[row] * (compsizes[row] - 1) / 2;
		// allocate current rhs and sense
		calloc_s(rhs, ncomp, double);
		calloc_s(sense, ncomp, char);
		calloc_s(rowname, ncomp, char*);
		// setup rhs, sense and rowname
		for (int row = 0; row < ncomp; row++)
		{
			rhs[row] = compsizes[row]-1;
			sense[row] = 'L';
			calloc_s(rowname[row], 30, char);
			sprintf(rowname[row], "SEC(%d)", nrows + row+1);
		}

		// ****************************** PREPARE COEFFS *****************************
		// allocate the coefficients arrays
		calloc_s(rowlist, tot_coeff_num, int);
		calloc_s(collist, tot_coeff_num, int);
		calloc_s(vallist, tot_coeff_num, double);
		// helper variables
		int arr_idx = 0;
		calloc_s(visited, g->nnodes, char);
		// for node, prepare its values
		for (int i = 0; i < g->nnodes; i++)
		{
			// iterate over the subtour to add constraints
			for (int j = succ[i]; j != i; j = succ[j])
			{
				if (!visited[j])
				{
					if (arr_idx >= tot_coeff_num) print_error("arr_idx out of bound", ERR_GENERIC_INCONSIST);
					// set values for element x_ij
					rowlist[arr_idx] = nrows + comp[i];
					collist[arr_idx] = xpos(i, j, g->nnodes);
					vallist[arr_idx] = 1.0;
					arr_idx++;
				}
			}
			visited[i] = 1;
		}

		// **************************** CREATE CONSTRAINTS ***************************
		// create empty rows (1 for each conn component)
		if (CPXnewrows(env, lp, ncomp, rhs, sense, NULL, rowname)) print_error("wrong CPXnewrows [SEC]", ERR_CPLEX);
		// fill coefficients (1 for each node)
		if (error = CPXchgcoeflist(env, lp, tot_coeff_num, rowlist, collist, vallist))
		{
			printf("CPX error %d\n", error);
			print_error("wrong CPXchgcoeflist [SEC]", ERR_CPLEX);
		}
		// log the addition of constraints
		log_line_ext(VERBOSITY, LOGLVL_INFO, "Added %d SE constraints with %d coefficients", ncomp, tot_coeff_num);

		// ********************************* CLEANUP *********************************
		// free names
		for (int row = 0; row < ncomp; row++) free(rowname[row]);
		// free arrays for rows
		free(rhs);
		free(sense);
		free(rowname);
		// free arrays for coeffs
		free(rowlist);
		free(collist);
		free(vallist);
		// free auxiliary arrays
		free(compsizes);
		free(visited);

		// ************************* SOLVE MODEL WITH NEW SEC ************************
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
	}


	// ********** FINAL CLEANUP **********
	// free solution
	free(xstar);
	// free auxiliary data structures
	free(succ);
	free(comp);
	// ***********************************

}

/* **************************************************************************************************
*						SOLUTION USING UNDIRECTED GRAPHS MODELS
************************************************************************************************** */
void solve_symmetric_tsp(graph* g, modeltype mt, CPXENVptr env, CPXLPptr lp)
{
	// if the tsptype of the model is not asymmetric, throw error
	if (model_tsptype(mt) != TSP_SYMM) print_error("", ERR_WRONG_TSP_PROCEDURE);

	// build symmetric TSP model OR use symmetric TSP procedure
	switch (model_archetype(mt))
	{
	case SY_BEND:
		solve_benders(g, mt, env, lp);
		break;
	default:
		print_error("symmetric variant", ERR_MODEL_NOT_IMPL);
	}
}


/* **************************************************************************************************
*						BASE MODEL FOR DIRECTED GRAPHS
************************************************************************************************** */
void build_model_base_directed(graph* g, modeltype mt, CPXENVptr env, CPXLPptr lp)
{

	// ********************************* SETUP *********************************
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

	// ************************ ADD ROWS (CONSTRAINTS) ************************
	// add the degree constraints
	for (int h = 0; h < g->nnodes; h++)  // degree constraints
	{
		// get row number
		int lastrow = CPXgetnumrows(env, lp);
		// set degree of contraint
		double rhs = 1.0;
		// 'E' for equality constraint 
		char sense = 'E';
		// OUTGOING CONSTRAINT
		// define name of constraint
		sprintf(cname, "degree_out(%d)", h + 1);
		// add empty constraint in CPX
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, &ptr_cname)) print_error("wrong CPXnewrows [degree]", ERR_CPLEX);
		// INGOING CONSTRAINT
		// define name of constraint
		sprintf(cname, "degree_in(%d)", h + 1);
		// add empty constraint in CPX
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, &ptr_cname)) print_error("wrong CPXnewrows [degree]", ERR_CPLEX);
		// fill the added constraints
		for (int i = 0; i < g->nnodes; i++)
		{
			// no constraint with itself
			if (i == h) continue;
			// constraint with node i
			if (CPXchgcoef(env, lp, lastrow, xxpos(h, i, g->nnodes), 1.0)) print_error("wrong CPXchgcoef [degree]", ERR_CPLEX);
			if (CPXchgcoef(env, lp, lastrow+1, xxpos(i, h, g->nnodes), 1.0)) print_error("wrong CPXchgcoef [degree]", ERR_CPLEX);
		}
	}
}


/* **************************************************************************************************
*			MILLER, TUCKER AND ZEMLIN COMPACT MODEL (MTZ)
************************************************************************************************** */
void build_model_mtz(graph* g, modeltype mt, CPXENVptr env, CPXLPptr lp)
{
	// construct the base model for directed graphs
	build_model_base_directed(g, mt, env, lp);
	
	int model_v = model_variant(mt);
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
	int izero = 0;
	int index[3];
	double value[3];
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
			int numrows;
			switch (model_v)
			{
			case 0:
				numrows = CPXgetnumrows(env, lp);
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, &ptr_cname)) print_error("wrong CPXlazyconstraints() for u-consistency", ERR_CPLEX);
				for (int k = 0; k < nnz; k++)
				{
					if (CPXchgcoef(env, lp, numrows, index[k], value[k])) print_error("wrong CPXchgcoef [degree]", ERR_CPLEX);
				}
				break;
			case 1:
			case 2:
				if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, &ptr_cname)) print_error("wrong CPXlazyconstraints() for u-consistency", ERR_CPLEX);
				break;
			}
			
		}
	}
	if (model_v == 2)
	{
		// add lazy constraints 1.0 * x_ij + 1.0 * x_ji <= 1, for each arc (i,j) with i < j
		rhs = 1.0;
		char sense = 'L';
		nnz = 2;
		for (int i = 0; i < g->nnodes; i++)
		{
			for (int j = i + 1; j < g->nnodes; j++)
			{
				sprintf(cname, "SEC on node pair (%d,%d)", i + 1, j + 1);
				index[0] = xxpos(i, j, g->nnodes);
				value[0] = 1.0;
				index[1] = xxpos(j, i, g->nnodes);
				value[1] = 1.0;
				if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, &ptr_cname)) print_error("wrong CPXlazyconstraints on 2-node SECs", ERR_CPLEX);
			}
		}
	}
}

/* **************************************************************************************************
*			GAVISH AND GRAVES SINGLE COMMODITY FLOW MODEL (GG)
************************************************************************************************** */
void build_model_gg(graph* g, modeltype mt, CPXENVptr env, CPXLPptr lp)
{
	// construct the base model for directed graphs
	build_model_base_directed(g, mt, env, lp);

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
	int izero = 0;
	double rhs = 0.0;
	char sense = 'E';
	int index[2];
	double value[2];
	int nnz = 2;
	// add lazy constraints (n-1) * x_1j - 1.0 * y_1j == 0 for each j but the first node
	for (int j = 1; j < g->nnodes; j++)
	{
		sprintf(cname, "commodity quantity on arc (%d,%d)", 1, j + 1);
		index[0] = xxpos(0, j, g->nnodes);
		value[0] = g->nnodes - 1.0;
		index[1] = ypos(0, j, g->nnodes);
		value[1] = -1.0;
		if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, &ptr_cname)) print_error("wrong CPXlazyconstraints on comm qty", ERR_CPLEX);
	}

	sense = 'G';
	// add lazy constrains (n-2) * x_ij - 1.0 * y_ij >= 0.0, for each arc (i,j) with i,j=/=0
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
			if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, &ptr_cname)) print_error("wrong CPXlazyconstraints on comm qty", ERR_CPLEX);
		}
	}

	int* index_flow;
	double* value_flow;
	// allocate arrays for holding indices and values for flow constraints
	if (!(index_flow = (int*)calloc(2 * (g->nnodes-1), sizeof(int)))) print_error("index_flow", ERR_NO_MEM_FOR_ALLOC);
	if (!(value_flow = (double*)calloc(2 * (g->nnodes-1), sizeof(double)))) print_error("value_flow", ERR_NO_MEM_FOR_ALLOC);

	sense = 'E';
	rhs = 1.0;
	nnz = 2 * (g->nnodes - 1);

	// add lazy constrains on flow through nodes
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
		if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index_flow, value_flow, &ptr_cname)) print_error("wrong CPXlazyconstraints on flow constraints", ERR_CPLEX);
	}

	free(index_flow);
	free(value_flow);

}

/* **************************************************************************************************
*						SOLUTION USING DIRECTED GRAPHS MODELS
************************************************************************************************** */
void solve_asymmetric_tsp(graph* g, modeltype mt, CPXENVptr env, CPXLPptr lp)
{
	// if the tsptype of the model is not asymmetric, throw error
	if (model_tsptype(mt) != TSP_ASYMM) print_error("", ERR_WRONG_TSP_PROCEDURE);

	// build asymmetric TSP model
	switch (model_archetype(mt))
	{
	case AS_MTZ:
		build_model_mtz(g, mt, env, lp);
		break;
	case AS_GG:
		build_model_gg(g, mt, env, lp);
		break;
	default:
		print_error("asymmetric variant", ERR_MODEL_NOT_IMPL);
	}

	// solve the model
	int error = 0;
	if (error = CPXmipopt(env, lp))
	{
		printf("CPX error %d\n", error);
		print_error("CPXmipopt() Asymmetric", ERR_CPLEX);
	}

}
