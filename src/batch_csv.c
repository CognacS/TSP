#include "../batch_csv.h"

void reorder_grid_csv(grid* p_grid)
{
	csvaxis order[] = { CELL, COLS, ROWS };
	csvaxis curr_axis;
	int axis_num = 3;

	gridparam aux_param;
	gridparam* param_list = p_grid->grid_params;
	int params_num = p_grid->params_num;
	
	for (int axis_idx = 0; axis_idx < axis_num; axis_idx++)
	{
		// search for the current axis
		curr_axis = order[axis_idx];

		// scan for a parameter with different axis from the current one
		for (int i = 0; i < params_num; i++)
		{
			if (param_list[i].axis != curr_axis)
			{
				int j;
				// look for a parameter with the current axis to switch
				for (j = i + 1; param_list[j].axis == curr_axis && j < params_num; j++);
				// if none was found then continue to the next axis
				if (j == params_num) break;
				// if one was found, then switch the two elements and go on
				printf("SWITCH %d<->%d\n", i, j);
				aux_param = param_list[i];
				param_list[i] = param_list[j];
				param_list[j] = aux_param;
			}
			
		}
	}
}