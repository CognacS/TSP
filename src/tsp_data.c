#include "../tsp_data.h"

/* ***************************************************************************************************
*						DATASTRUCTURES PRINT FUNCTIONS
*************************************************************************************************** */


void print_graph(graph* inst_graph)
{
	printf("Graph of %d nodes, distances are int=%d, coords:\n", inst_graph->nnodes, inst_graph->integer_costs);
	/*
	for (int i = 0; i < inst_graph->nnodes; i++)
	{
		printf("(%lf, %lf), ", inst_graph->xcoord[i], inst_graph->ycoord[i]);

	}
	printf("\n");
	*/
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
void print_model(model* inst_model)
{
	printf("Model:\n"
		"\tncols = %d\n",
		inst_model->ncols
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


void fill_inst_default(instance* inst)
{
	// define the parameters for the instance
	graph* g = &inst->inst_graph;
	params* p = &inst->inst_params;
	global_data* gd = &inst->inst_global_data;
	model* m = &inst->inst_model;

	// ************ DEFAULT PARAMETERS DEFINITION ************
	g->integer_costs = DEF_INTEGER_COSTS;
	g->xcoord = NULL;
	g->ycoord = NULL;
	g->tr_xcoord = NULL;
	g->tr_ycoord = NULL;

	p->model_type = DEF_MODEL_TYPE;
	strcpy(p->input_file, DEF_INPUT_FILE);
	strcpy(p->batch_file, DEF_BATCH_FILE);
	p->timelimit = DEF_TIMELIMIT;
	p->randomseed = DEF_RANDOMSEED;
	p->max_nodes = DEF_MAX_NODES;

	gd->xstar = DEF_XSTAR;

	m->ncols = DEF_NCOLS;
	// *******************************************************

}

/* ***************************************************************************************************
*						DATASTRUCTURES DESTROYERS
*************************************************************************************************** */

void free_graph(graph* inst_graph)
{
	free_s(inst_graph->xcoord);
	free_s(inst_graph->ycoord);
	free_s(inst_graph->tr_xcoord);
	free_s(inst_graph->tr_ycoord);
}

void free_global_data(global_data* inst_global)
{
	free_s(inst_global->xstar);
}

void free_solution(Solution* sol)
{
	free(sol->xstar);
	free(sol->succ);
}

void free_instance(instance* inst)
{
	free_graph(&inst->inst_graph);
	free_global_data(&inst->inst_global_data);
}


/* ***********************************************************************************
*						AUXILIARY STRUCTURES AND FUNCTIONS
*********************************************************************************** */

LinkedList* newLinkedList()
{
	LinkedList* ll;	malloc_s(ll, LinkedList);
	ll->start = NULL; ll->end = NULL;
	return ll;
}

int LL_is_empty(LinkedList* list)
{
	return list->start == NULL;
}

void LL_add_value(LinkedList* list, double x, double y)
{
	Cell* cell = NULL;	malloc_s(cell, Cell);
	cell->x = x;
	cell->y = y;
	cell->next = NULL;
	// if list is not initialized
	if (list->end == NULL) list->end = list->start = cell;
	// else add at the end
	else {
		list->end->next = cell;
		list->end = cell;
	}
}

void LL_free(LinkedList* list)
{
	if (list->start != NULL)
	{
		// iterate through all the list and free
		for (Cell* cell = list->start, *next_cell = list->start->next;
			next_cell != NULL;
			cell = next_cell, next_cell = next_cell->next)
		{
			free(cell);
		}
		free(list->end);
	}
	free(list);
}