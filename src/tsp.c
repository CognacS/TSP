#include "../include/tsp.h"

/* ***************************************************************************************************
*						DATASTRUCTURES PRINT FUNCTIONS
*************************************************************************************************** */

void print_params(Params* inst_p)
{
	printf("Parameters:\n"
		"\tmodel_type = %d\n"
		"\ttimelimit  = %f\n"
		"\tinput_file = %s\n"
		"\tbatch_file = %s\n"
		"\trandomseed = %d\n"
		"\tcutoff = %f\n"
		"\tmax_nodes = %d\n"
		"\theuristic_code = %s\n",
		inst_p->model_type,
		inst_p->timelimit,
		inst_p->input_file,
		inst_p->batch_file,
		inst_p->randomseed,
		inst_p->cutoff,
		inst_p->max_nodes,
		inst_p->heuristic_code
	);
}


void print_global_data(GlobalData* inst_global)
{
	printf("Global Data:\n"
		"\ttstart = %f\n"
		"\ttexec = %f\n"
		"\tzbest= %f\n"
		"\tlbbest = %f\n",
		inst_global->tstart,
		inst_global->texec,
		inst_global->zbest,
		inst_global->lbbest
	);
}
void print_model(Model* inst_model)
{
	printf("Model:\n"
		"\tncols = %d\n",
		inst_model->ncols
	);
}

void print_instance(Instance* inst)
{
	printf("********** INSTANCE REPORT **********\n");

	print_graph(&inst->inst_graph);
	print_params(&inst->inst_params);
	print_global_data(&inst->inst_global_data);
	print_model(&inst->inst_model);

	printf("*************************************\n");

}


void fill_inst_default(Instance* inst)
{
	// define the parameters for the instance
	Graph* g = &inst->inst_graph;
	Params* p = &inst->inst_params;
	GlobalData* gd = &inst->inst_global_data;
	Model* m = &inst->inst_model;

	// ************ DEFAULT PARAMETERS DEFINITION ************
	default_graph(g);

	p->model_type = DEF_MODEL_TYPE;
	strcpy(p->input_file, DEF_INPUT_FILE);
	strcpy(p->batch_file, DEF_BATCH_FILE);
	p->timelimit = DEF_TIMELIMIT;
	p->randomseed = DEF_RANDOMSEED;
	p->max_nodes = DEF_MAX_NODES;

	gd->xstar = NULL;

	m->ncols = DEF_NCOLS;
	// *******************************************************

}

/* ***************************************************************************************************
*						DATASTRUCTURES DESTROYERS
*************************************************************************************************** */

void free_global_data(GlobalData* inst_global)
{
	free_s(inst_global->xstar);
}

void free_instance(Instance* inst)
{
	free_graph(&inst->inst_graph);
	free_global_data(&inst->inst_global_data);
}

/* ***************************************************************************************************
*						DATASTRUCTURES LOGGING
*************************************************************************************************** */

void log_datastruct(void* object, DataType type, int runlvl, int loglvl)
{
	if (runlvl >= loglvl)
	{
		switch (type)
		{
		case TYPE_GRAPH:
			print_graph((Graph*)object);
			break;
		case TYPE_PARAM:
			print_params((Params*)object);
			break;
		case TYPE_GLOB:
			print_global_data((GlobalData*)object);
			break;
		case TYPE_MODEL:
			print_model((Model*)object);
			break;
		case TYPE_INST:
			print_instance((Instance*)object);
			break;
		default:
			print_warn(WARN_WRONG_DATASTRUCT, NULL);
		}
	}
}

/* ***************************************************************************************************
*						TME HANDLING
*************************************************************************************************** */

double residual_time(Instance* inst)
{
	return inst->inst_global_data.tstart + inst->inst_params.timelimit - second();
}

int time_limit_expired(Instance* inst)
{
	GlobalData* gd = &inst->inst_global_data;
	Params* p = &inst->inst_params;

	double tspan = second() - gd->tstart;
	if (tspan > p->timelimit)
	{
		print_warn_ext(WARN_EXPIRED_TIMELIMIT, "limit of %10.1lf sec.s expired after %10.1lf sec.s", p->timelimit, tspan);
		return 1;
	}
	return 0;
}