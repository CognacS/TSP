#include "../batch_tool.h"

int next(grid* p_grid, instance* inst)
{
	int curr_coord;
	gridparam curr_param;
	gridparamtype curr_type;
	int increase = 1;

	for (int i = 0; i < p_grid->params_num; i++)
	{
		curr_coord = p_grid->coords[i];
		curr_param = p_grid->params[i];
		curr_type = curr_param.type;

		switch (curr_type)
		{
		case MODEL_TYPE:
			inst->inst_params.model_type = *(modeltype*)(curr_param.values[curr_coord]);
			break;
		case INPUT_FILE:
			strcpy(inst->inst_params.input_file, *(char**)(curr_param.values[curr_coord]));
			break;
		}

		if (increase)
		{
			p_grid->coords[i] = (curr_coord + 1) % curr_param.values_num;
			increase = p_grid->coords[i] == 0;
		}
	}

	return !increase;

}


void free_gridparam(gridparam* params)
{
	for (int i = 0; i < params->values_num; i++)
	{
		free(params->values[i]);
		free(params->labels[i]);
	}
	free(params->values);
	free(params->labels);
}
void free_grid(grid* p_grid)
{
	for (int i = 0; i < p_grid->params_num; i++)
	{
		free_gridparam(&p_grid->params[i]);
	}
	free(p_grid->params);
	free(p_grid->coords);
	free(p_grid->out_file_is_row);
}
void free_batchtool(batchtool* bt)
{
	free_grid(&bt->p_grid);
}