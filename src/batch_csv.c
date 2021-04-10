#include "../batch_csv.h"

void reorder_grid_csv(grid* p_grid)
{
	csvaxis order[] = { CELL, COLS, ROWS };
	csvaxis curr_axis;
	int axis_num = 3;

	gridparam aux_param;
	gridparam* param_list = p_grid->grid_params;
	int params_num = p_grid->params_num;
	
	int last_switch = -1;

	for (int axis_idx = 0; axis_idx < axis_num; axis_idx++)
	{
		// search for the current axis
		curr_axis = order[axis_idx];

		// scan for a parameter with different axis from the current one
		for (int i = last_switch+1; i < params_num; i++)
		{
			if (param_list[i].axis != curr_axis)
			{
				int j;
				// look for a parameter with the current axis to switch
				for (j = i + 1; param_list[j].axis != curr_axis && j < params_num; j++);
				// if none was found then continue to the next axis
				if (j == params_num) break;
				// if one was found, then switch the two elements and go on
				aux_param = param_list[i];
				param_list[i] = param_list[j];
				param_list[j] = aux_param;
				last_switch = i;
			}
			
		}
	}
}

void open_file_csv(csv_batchtool* bt)
{
	// open the csv file "write-only"
	bt->csv_fp = fopen(bt->csv_file, "w");
	// check if the file exists
	if (bt->csv_fp == NULL) print_error(ERR_INPUT_NOT_EXISTS, bt->csv_file);
}

void register_time_csv(csv_batchtool* bt, double seconds)
{

}