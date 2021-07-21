#include "../include/cpx_utility.h"

int upos(int i, int nnodes)
{
	// if index is of first node throw error
	if (i == 0) print_error(ERR_INVALID_FUNC_ARGS, "got node 0 in upos");
	// else return upos
	return i - 1 + (nnodes - 1) * nnodes;
}

int ypos(int i, int j, int nnodes)
{
	// if indices are equal throw error
	if (i == j) print_error(ERR_INVALID_FUNC_ARGS, "got equal indices in ypos");
	// else return ypos
	return xxpos(i, j, nnodes) + (nnodes - 1) * nnodes;
}

/*********************************************************************************************************************************************************/
void cpx_add_cut(void* env, void* lp,
	int nnz, double rhs, char sense,
	int* index, double* value, int purgeable, char* name, int local)
	/*********************************************************************************************************************************************************/
{

	/*
		unified driver to add constraints to the LP; usage:
			cpx_add_cut(env,	 lp,	nnz, rhs, sense, index, value, name, -1,	CUT_STATIC);   			// add static constr to LP
			cpx_add_cut(env,	 lp,	nnz, rhs, sense, index, value, name, -1,	CUT_LAZY);   			// add lazy constr to LP
			cpx_add_cut(context, NULL,	nnz, rhs, sense, index, value, NULL, local, purgeable); 			// from callbacks, add user cut
			cpx_add_cut(context, NULL,	nnz, rhs, sense, index, value, NULL, -1,	CUT_CALLBACK_REJECT);	// from callbacks, reject candidate
	*/

	int izero = 0;

	// REFERENCE: https://www.ibm.com/docs/en/icos/20.1.0?topic=c-cpxxcallbackaddusercuts-cpxcallbackaddusercuts
	if (purgeable >= 0) 	// add the cut from a callback
	{
		if (CPXcallbackaddusercuts((CPXCALLBACKCONTEXTptr)env, 1, nnz, &rhs, &sense, &izero, index, value, &purgeable, &local))
			print_error(ERR_ADD_CUT, "error on CPXcallbackaddusercuts (>=0)");
		return;
	}

	// REFERENCE: https://www.ibm.com/docs/en/icos/20.1.0?topic=cpxxaddrows-cpxaddrows
	if (purgeable == CUT_STATIC) 	// statically add the cut to the original model (with casting (CPXLPptr) cbdata)
	{
		if (CPXaddrows((CPXENVptr)env, (CPXLPptr)lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &name))
			print_error(ERR_ADD_CUT, "error on CPXaddrows (CUT_STATIC)");
		return;
	}

	// REFERENCE: https://www.ibm.com/docs/en/icos/20.1.0?topic=cpxxaddlazyconstraints-cpxaddlazyconstraints
	if (purgeable == CUT_LAZY) 	// add the lazycut cut to the original model (with casting (CPXLPptr) cbdata)
	{
		if (CPXaddlazyconstraints((CPXENVptr)env, (CPXLPptr)lp, 1, nnz, &rhs, &sense, &izero, index, value, &name))
			print_error(ERR_ADD_CUT, "error on CPXaddlazyconstraints (CUT_LAZY)");
		return;
	}

	// REFERENCE: https://www.ibm.com/docs/es/cofz/12.10.0?topic=c-cpxxcallbackrejectcandidate-cpxcallbackrejectcandidate
	if (purgeable == CUT_CALLBACK_REJECT) 	// add the usercut *local* to the original model (with casting (CPXLPptr) cbdata)
	{
		if (CPXcallbackrejectcandidate((CPXCALLBACKCONTEXTptr)env, 1, nnz, &rhs, &sense, &izero, index, value))
			print_error(ERR_ADD_CUT, "error on CPXcallbackrejectcandidate (CUT_CALLBACK_REJECT)");
		return;
	}

	print_error(ERR_ADD_CUT, "purgeable flag unknown");
}

void cpx_setup_cplex(OptData* optdata)
{
	Instance* inst = optdata->inst;
	Params* p = &(inst->inst_params);

	// allocate datastruct
	if (optdata->cpx != NULL) print_error(ERR_OPT_CPLEX_REDEFINITION, NULL);
	malloc_s(optdata->cpx, CplexData);

	// allocate env and lp
	int error;
	CplexData* cpx = optdata->cpx;
	cpx->env = CPXopenCPLEX(&error);
	cpx->lp = CPXcreateprob(optdata->cpx->env, &error, "TSP");

	// setup parameters
	CPXsetintparam(cpx->env, CPX_PARAM_RANDOMSEED, p->randomseed);
	CPXsetintparam(cpx->env, CPX_PARAM_NODELIM, p->max_nodes);
	CPXsetdblparam(cpx->env, CPX_PARAM_EPGAP, p->cutoff);
	CPXsetdblparam(cpx->env, CPX_PARAM_EPINT, 0.0);
	if (VERBOSITY >= LOGLVL_CPLEXLOG)
	{
		CPXsetintparam(cpx->env, CPX_PARAM_SCRIND, CPX_ON);
		CPXsetintparam(cpx->env, CPXPARAM_MIP_Display, 5);
	}

	cpx_timelimit(optdata, p->timelimit);
}

void cpx_close_cplex(OptData* optdata)
{
	CplexData* cpx = optdata->cpx;
	if (cpx == NULL) return;
	// cleanup cplex
	CPXfreeprob(cpx->env, &(cpx->lp));
	CPXcloseCPLEX(&(cpx->env));
	free_s(optdata->cpx);
}

int cpx_solution_available(OptData* optdata)
{
	double zz;
	if (CPXgetobjval(optdata->cpx->env, optdata->cpx->lp, &zz)) return 0;
	return 1;
}

void cpx_warmstart(OptData* optdata, double* xstar)
{
	int ncols = optdata->inst->inst_model.ncols;

	int izero = 0;
	int effortlevel = CPX_MIPSTART_AUTO;
	int* varindices; arr_malloc_s(varindices, ncols, int);
	for (int i = 0; i < ncols; i++) varindices[i] = i;
	int curr_mipstarts = CPXgetnummipstarts(optdata->cpx->env, optdata->cpx->lp);
	log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Number of mip starts: %d", curr_mipstarts);
	if (curr_mipstarts > 0) CPXdelmipstarts(optdata->cpx->env, optdata->cpx->lp, 0, curr_mipstarts - 1);
	CPXaddmipstarts(optdata->cpx->env, optdata->cpx->lp,
		1, ncols, &izero, varindices, xstar, &effortlevel, NULL);
	free(varindices);
}

void cpx_extract_sol_obj_lb(OptData* optdata, char* zone_name)
{
	int ncols = optdata->inst->inst_model.ncols;
	GlobalData* gd = &(optdata->inst->inst_global_data);

	int error;
	// get the optimal solution
	if (error = CPXgetx(optdata->cpx->env, optdata->cpx->lp, gd->xstar, 0, ncols - 1))
		print_error_ext(ERR_CPLEX, "CPXgetx() %s, CPX error: %d", zone_name, error);

	// get the obj value
	if (error = CPXgetobjval(optdata->cpx->env, optdata->cpx->lp, &gd->zbest))
		print_error_ext(ERR_CPLEX, "CPXgetobjval() %s, CPX error: %d", zone_name, error);

	// get the lower bound
	if (error = CPXgetbestobjval(optdata->cpx->env, optdata->cpx->lp, &gd->lbbest))
		print_error_ext(ERR_CPLEX, "CPXgetbestobjval() %s, CPX error: %d", zone_name, error);
}

void cpx_extract_sol_obj(OptData* optdata, Solution* sol, char* zone_name)
{
	int ncols = optdata->inst->inst_model.ncols;
	int error;

	// check for pointer
	if (sol->xstar == NULL) print_error(ERR_OPT_SOLPTR_INCONSISTENT, "xstar");
	// get the optimal solution
	if (error = CPXgetx(optdata->cpx->env, optdata->cpx->lp, sol->xstar, 0, ncols - 1))
		print_error_ext(ERR_CPLEX, "CPXgetx() %s, CPX error: %d", zone_name, error);

	// get the obj value
	if (error = CPXgetobjval(optdata->cpx->env, optdata->cpx->lp, &sol->cost))
		print_error_ext(ERR_CPLEX, "CPXgetobjval() %s, CPX error: %d", zone_name, error);
}

void cpx_timelimit(OptData* optdata, double timelimit)
{
	double residual_time = optdata->inst->inst_global_data.tstart + optdata->inst->inst_params.timelimit - second();
	if (residual_time < 0.0) residual_time = 0.0;
	double computed_timelimit = min(residual_time, timelimit);
	CPXsetintparam(optdata->cpx->env, CPX_PARAM_CLOCKTYPE, 2);
	CPXsetdblparam(optdata->cpx->env, CPX_PARAM_TILIM, computed_timelimit);	// real time
}