#include "../tsp_data.h"

/* ***************************************************************************************************
*						DATASTRUCTURES PRINT FUNCTIONS
*************************************************************************************************** */


void print_graph(graph* inst_graph)
{
	printf("Graph of %d nodes, distances are int=%d, coords:\n", inst_graph->nnodes, inst_graph->integer_costs);

	for (int i = 0; i < inst_graph->nnodes; i++)
	{
		printf("(%lf, %lf), ", inst_graph->xcoord[i], inst_graph->ycoord[i]);

	}
	printf("\n");
}

void print_params(params* inst_p)
{
	printf("Parameters:\n"
		"\tmodel_type = %d\n"
		"\ttimelimit  = %f\n"
		"\tinput_file = %s\n"
		"\tbatch_file = %s\n"
		"\trandomseed = %d\n"
		"\tcutoff = %f\n"
		"\tmax_nodes = %d\n",
		inst_p->model_type,
		inst_p->timelimit,
		inst_p->input_file,
		inst_p->batch_file,
		inst_p->randomseed,
		inst_p->cutoff,
		inst_p->max_nodes
	);
}


void print_global_data(global_data* inst_global)
{
	printf("Global Data:\n"
		"\tz_best = %f\n",
		inst_global->z_best
	);
}
void print_model(model* inst_model)
{
	printf("Model:\n"
		"\tplc_holder = %d\n",
		inst_model->plc_holder
	);
}

void print_instance(instance* inst)
{
	printf("********** INSTANCE REPORT **********\n");

	print_graph(&inst->inst_graph);
	print_params(&inst->inst_params);
	print_global_data(&inst->inst_global_data);
	print_model(&inst->inst_model);

	printf("*************************************\n");

}

/* ***************************************************************************************************
*						DATASTRUCTURES DESTROYERS
*************************************************************************************************** */

void free_graph(graph* inst_graph)
{
	free(inst_graph->xcoord);
	free(inst_graph->ycoord);
	free(inst_graph->tr_xcoord);
	free(inst_graph->tr_ycoord);
}
void free_instance(instance* inst)
{
	free_graph(&inst->inst_graph);
}