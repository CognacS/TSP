#include "../mip_utility.h"

/*********************************************************************************************************************************************************/
void mip_add_cut(void* env, void* lp,
	int nnz, double rhs, char sense,
	int* index, double* value, int purgeable, char* name, int local)
/*********************************************************************************************************************************************************/
{

	/*
		unified driver to add constraints to the LP; usage:
			mip_add_cut(env,	 lp,	nnz, rhs, sense, index, value, name, -1,	CUT_STATIC);   			// add static constr to LP
			mip_add_cut(env,	 lp,	nnz, rhs, sense, index, value, name, -1,	CUT_LAZY);   			// add lazy constr to LP
			mip_add_cut(context, NULL,	nnz, rhs, sense, index, value, NULL, local, purgeable); 			// from callbacks, add user cut
			mip_add_cut(context, NULL,	nnz, rhs, sense, index, value, NULL, -1,	CUT_CALLBACK_REJECT);	// from callbacks, reject candidate
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

void mip_timelimit(CPXENVptr env, double timelimit, instance* inst)
{
	double residual_time = inst->inst_global_data.tstart + inst->inst_params.timelimit - second();
	if (residual_time < 0.0) residual_time = 0.0;
	CPXsetintparam(env, CPX_PARAM_CLOCKTYPE, 2);
	CPXsetdblparam(env, CPX_PARAM_TILIM, residual_time);					// real time
	CPXsetdblparam(env, CPX_PARAM_DETTILIM, TICKS_PER_SECOND * timelimit);	// ticks
}

int time_limit_expired(instance* inst)
{
	global_data* gd = &inst->inst_global_data;
	params* p = &inst->inst_params;

	double tspan = second() - gd->tstart;
	if (tspan > p->timelimit)
	{
		print_warn_ext(WARN_EXPIRED_TIMELIMIT, "limit of %10.1lf sec.s expired after %10.1lf sec.s", p->timelimit, tspan);
		return 1;
	}
	return 0;
}