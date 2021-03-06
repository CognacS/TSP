#include "../batch_tool.h"

int next_args_config(Grid* p_grid, char** new_argv)
{
	if (p_grid->end_reached) return 0;

	// ********************** WRITE VALUES ********************** 
	int curr_index;
	GridParam* curr_param;
	// for each parameter in the grid
	for (int i = 0; i < p_grid->params_num; i++)
	{
		// get the current value of the paramer
		curr_param = &p_grid->grid_params[i];
		// get the current index
		curr_index = p_grid->indices[i];

		// add the parameter and its value to the arguments array
		new_argv[2 * i] = curr_param->param_name;
		new_argv[2 * i + 1] = curr_param->values[curr_index];
	}

	// **************** INCREASE COORDINATES **************** 
	int increase = 1;
	int incr_idx = 0;
	// increase the current index if there is a carry
	for (incr_idx = 0; increase && incr_idx < p_grid->params_num; incr_idx++)
	{
		p_grid->indices[incr_idx] = (p_grid->indices[incr_idx] + 1) % p_grid->grid_params[incr_idx].values_num;
		increase = p_grid->indices[incr_idx] == 0;
	}
	p_grid->last_incridx = incr_idx - 1;

	// stop iterator if the increase carry reached the end
	if (increase) p_grid->end_reached = 1;
	
	return 1;
}

int next_inst_config(Grid* p_grid, Instance* inst)
{
	// allocate new arg arrays
	int new_argc = p_grid->params_num * 2;
	char** new_argv;
	calloc_s(new_argv, new_argc, char*);

	// get next configuration of argument if it exists
	int hasnext = next_args_config(p_grid, new_argv);
	// set parameters like a command line
	if (hasnext) parse_command_line(new_argc, new_argv, inst);
	
	// free and return
	free(new_argv);
	return hasnext;
}

void restart_grid(Grid* p_grid)
{
	p_grid->end_reached = 0;
	p_grid->started = 0;
	for (int i = 0; i < p_grid->params_num; i++)
	{
		p_grid->indices[i] = 0;
	}
}

void print_grid(Grid* p_grid)
{
	// allocate new arg arrays
	int new_argc = p_grid->params_num * 2;
	char** new_argv;
	calloc_s(new_argv, new_argc, char*);

	// iterate over all grid argument command lines
	restart_grid(p_grid);
	while (next_args_config(p_grid, new_argv))
	{
		for (int i = 0; i < new_argc; i++)
		{
			printf("%s ", new_argv[i]);
		}
		printf("\n");
	}
	restart_grid(p_grid);
	free(new_argv);

}

void print_BatchTool(BatchTool* bt)
{
	printf("BATCH TOOL\n");
	printf("\t - input file: %s\n", bt->input_file);
	print_grid(&bt->p_grid);
}


void free_GridParam(GridParam* params)
{
	for (int i = 0; i < params->values_num; i++)
	{
		free(params->values[i]);
		free(params->labels[i]);
	}
	free(params->values);
	free(params->labels);
}
void free_grid(Grid* p_grid)
{
	for (int i = 0; i < p_grid->params_num; i++)
	{
		free_GridParam(&p_grid->grid_params[i]);
	}
	free(p_grid->grid_params);
	free(p_grid->indices);
}
void free_BatchTool(BatchTool* bt)
{
	free_grid(&bt->p_grid);
}