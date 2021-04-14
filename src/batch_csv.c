#include "../batch_csv.h"

void reorder_grid_csv(csv_batchtool* bt)
{
	grid* p_grid = &bt->bt.p_grid;
	csvaxis order[] = { CSV_CELL, CSV_COLS, CSV_ROWS };
	csvaxis curr_axis;
	int axis_num = 3;

	gridparam aux_param;
	gridparam* param_list = p_grid->grid_params;
	int params_num = p_grid->params_num;
	
	// elements are correct up to this index
	int corr_idx = -1;

	for (int axis_idx = 0; axis_idx < axis_num; axis_idx++)
	{
		// search for the current axis
		curr_axis = order[axis_idx];

		// scan for a parameter with different axis from the current one
		for (int i = corr_idx +1; i < params_num; i++)
		{
			if (param_list[i].axis != curr_axis)
			{
				int j;
				// look for a parameter with the current axis to switch
				for (j = i + 1; j < params_num && param_list[j].axis != curr_axis; j++);
				// if none was found then continue to the next axis
				if (j == params_num) break;
				// if one was found, then switch the two elements and go on
				aux_param = param_list[i];
				param_list[i] = param_list[j];
				param_list[j] = aux_param;
			}
			// update correct elements index
			corr_idx = i;
		}
	}

	bt->reordered = 1;
}

void open_file_csv(csv_batchtool* bt)
{
	// if the grid was not reordered, throw error
	if (!bt->reordered) print_error(ERR_CSV_NOT_REORDERED, NULL);

	// open the csv file "write-only"
	bt->csv_fp = fopen(bt->csv_file, "w");
	// check if the file was opened
	if (bt->csv_fp == NULL) print_error(ERR_CSV_CANNOT_OPEN, bt->csv_file);

	// ********* print first line *********
	grid* p_grid = &bt->bt.p_grid;
	// iterate through all cols
	int cols_start;
	int cols_end;
	// find start of cols
	for (cols_start = 0;
		p_grid->grid_params[cols_start].axis != CSV_COLS &&
		cols_start < p_grid->params_num; cols_start++);

	// find end of cols
	for (cols_end = cols_start;
		p_grid->grid_params[cols_end].axis == CSV_COLS &&
		cols_end < p_grid->params_num; cols_end++);
	// get size of indices to consider
	int size = cols_end - cols_start;

	// get number of columns
	int col_num = 0;
	for (int i = cols_start; i < cols_end; i++)
	{
		col_num += p_grid->grid_params[i].values_num;
	}
	// print number of columns
	fprintf(bt->csv_fp, "%d", col_num);

	int* indices;	calloc_s(indices, size, int);

	// iterate through each label and write them to the csv
	int increase;
	do
	{
		// initialize new cell with a comma
		fprintf(bt->csv_fp, ",");
		// for each col parameter
		for (int i = 0; i < size; i++)
		{
			// get param at position cols_start + i
			// get its label at index i
			int idx = indices[i];
			fprintf(bt->csv_fp, "%s+", p_grid->grid_params[cols_start + i].labels[idx]);

		}
		// delete '+' sign at the end with a backspace
		fseek(bt->csv_fp, -1, SEEK_CUR);

		// increase indices
		increase = 1;
		for (int i = 0; increase && i < size; i++)
		{
			// update index
			indices[i] = (indices[i] + 1) % p_grid->grid_params[cols_start + i].values_num;
			// increase if there was an overflow
			increase = indices[i] == 0;
		}

	} while (!increase);

	// make first line
	newline_csv(bt);

	free(indices);
}

void close_file_csv(csv_batchtool* bt)
{
	fclose(bt->csv_fp);
}

void register_time_csv(csv_batchtool* bt, double seconds)
{
	grid* p_grid = &bt->bt.p_grid;
	// define cumulative time and number
	static double accum_time;
	static int	  samples_num;
	// update value of current cell
	accum_time += seconds;
	samples_num++;

	// get last update index
	int last_incridx = p_grid->last_incridx;

	// get type of axis to update
	csvaxis axis_update = p_grid->grid_params[last_incridx].axis;

	// if update is out of cell, print the value of the cell
	if (axis_update != CSV_CELL && samples_num > 0)
	{
		fprintf(bt->csv_fp, ",%f", accum_time / samples_num);
		accum_time = 0;
		samples_num = 0;
	}

	// when updating to a new row, print a newline and the next name (if not end reached)
	if (axis_update == CSV_ROWS)
	{
		newline_csv(bt);
	}

}

void newline_csv(csv_batchtool* bt)
{
	grid* p_grid = &bt->bt.p_grid;

	fprintf(bt->csv_fp, "\n");

	if (!p_grid->end_reached)
	{
		// construct the new row label
		for (int i = 0; i < p_grid->params_num; i++)
		{
			if (p_grid->grid_params[i].axis == CSV_ROWS)
			{
				int idx = p_grid->indices[i];
				fprintf(bt->csv_fp, "%s+", p_grid->grid_params[i].labels[idx]);
			}
		}
		// delete '+' sign at the end with a backspace
		fseek(bt->csv_fp, -1, SEEK_CUR);
	}
}