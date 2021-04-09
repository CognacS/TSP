#include "../mip_utility.h"

/*********************************************************************************************************************************************************/
void mip_add_cut(void* env, void* cbdata,
	int wherefrom, int nnz, double rhs, char sense,
	int* index, double* value, char* cname, int purgeable)
/*********************************************************************************************************************************************************/
{

	/*
		unified driver to add constraints to the LP; usage:
			mip_add_cut(env,	lp,		-1,			nnz, rhs, sense, index, value, cname,	CUT_STATIC);   			// -1 means static cuts
			mip_add_cut(env,	lp,		-1,			nnz, rhs, sense, index, value, cname,	CUT_LAZY);   			// -2 means lazy cuts
			mip_add_cut(env,	lp,		-1,			nnz, rhs, sense, index, value, cname,	CUT_USER);   			// -3 means user cuts
			mip_add_cut(env,	cbdata, wherefrom,	nnz, rhs, sense, index, value, NULL,	purgeable); 			// from callbacks, user cut
			mip_add_cut(env,	cbdata, wherefrom,	nnz, rhs, sense, index, value, NULL,	CUT_CALLBACK_LOCAL); 	// from callbacks, local cut
			mip_add_cut(context,NULL,	-1,			nnz, rhs, sense, index, value, NULL,	CUT_CALLBACK_REJECT);	// from callbacks, reject candidate
	*/

	// REFERENCE: https://www.ibm.com/docs/en/icos/12.7.1.0?topic=cpxxcutcallbackadd-cpxcutcallbackadd
	if (purgeable >= 0) 	// add the cut from a callback
	{
		if (CPXcutcallbackadd((CPXCENVptr)env, cbdata, wherefrom, nnz, rhs, sense, index, value, purgeable)) print_error("1", ERR_ADD_CUT);
		return;
	}

	int izero = 0;

	// REFERENCE: https://www.ibm.com/docs/en/icos/20.1.0?topic=cpxxaddrows-cpxaddrows
	if (purgeable == CUT_STATIC) 	// statically add the cut to the original model (with casting (CPXLPptr) cbdata)
	{
		if (CPXaddrows((CPXENVptr)env, (CPXLPptr)cbdata, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname)) print_error("2", ERR_ADD_CUT);
		return;
	}

	// REFERENCE: https://www.ibm.com/docs/en/icos/20.1.0?topic=cpxxaddlazyconstraints-cpxaddlazyconstraints
	if (purgeable == CUT_LAZY) 	// add the lazycut cut to the original model (with casting (CPXLPptr) cbdata)
	{
		if (CPXaddlazyconstraints((CPXENVptr)env, (CPXLPptr)cbdata, 1, nnz, &rhs, &sense, &izero, index, value, &cname)) print_error("3", ERR_ADD_CUT);
		return;
	}

	// REFERENCE: https://www.ibm.com/docs/en/icos/20.1.0?topic=cpxxaddusercuts-cpxaddusercuts
	if (purgeable == CUT_USER) 	// add the usercut cut to the original model (with casting (CPXLPptr) cbdata)
	{
		if (CPXaddusercuts((CPXENVptr)env, (CPXLPptr)cbdata, 1, nnz, &rhs, &sense, &izero, index, value, &cname)) print_error("4", ERR_ADD_CUT);
		return;
	}

	// REFERENCE: https://www.ibm.com/docs/en/icos/12.10.0?topic=c-cpxxcutcallbackaddlocal-cpxcutcallbackaddlocal
	if (purgeable == CUT_CALLBACK_LOCAL) 	// add the usercut *local* to the original model (with casting (CPXLPptr) cbdata)
	{
		if (CPXcutcallbackaddlocal((CPXCENVptr)env, cbdata, wherefrom, nnz, rhs, sense, index, value)) print_error("local cuts", ERR_ADD_CUT);
		return;
	}

	// REFERENCE: https://www.ibm.com/docs/es/cofz/12.10.0?topic=c-cpxxcallbackrejectcandidate-cpxcallbackrejectcandidate
	if (purgeable == CUT_CALLBACK_REJECT) 	// add the usercut *local* to the original model (with casting (CPXLPptr) cbdata)
	{
		if (CPXcallbackrejectcandidate((CPXCALLBACKCONTEXTptr)env, 1, nnz, &rhs, &sense, &izero, index, value)) print_error("local cuts", ERR_ADD_CUT);
		return;
	}

	print_error("purgeable flag not known", ERR_ADD_CUT);
}

void mip_timelimit(CPXENVptr env, double timelimit, instance* inst)
{
	double residual_time = inst->inst_global_data.tstart + inst->inst_params.timelimit - second();
	if (residual_time < 0.0) residual_time = 0.0;
	CPXsetintparam(env, CPX_PARAM_CLOCKTYPE, 2);
	CPXsetdblparam(env, CPX_PARAM_TILIM, residual_time);					// real time
	CPXsetdblparam(env, CPX_PARAM_DETTILIM, TICKS_PER_SECOND * timelimit);	// ticks
}