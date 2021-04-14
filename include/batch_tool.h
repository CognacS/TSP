#ifndef BATCH_TOOL_H_

#define BATCH_TOOL_H_

#include "tsp_data.h"
#include "comm_parser.h"

#define MAX_SIZE_FILE_NAME 100
#define MAX_SIZE_PARAM_NAME 100

/* ***********************************************************************************
*						GRID SECTION
*********************************************************************************** */

// *********** definitions for csv file ***********
typedef enum
{
	CSV_ROWS,
	CSV_COLS,
	CSV_CELL
} csvaxis;

typedef struct
{
	char param_name[MAX_SIZE_PARAM_NAME];
	csvaxis axis;
	int values_num;
	char** labels;
	char** values;

} gridparam;

typedef struct
{
	int started;
	int end_reached;
	int params_num;
	int* indices;
	int last_incridx;
	gridparam* grid_params;
} grid;

/* ***********************************************************************************
*						BATCH TOOL SECTION
*********************************************************************************** */

typedef struct batch_tool_struct
{

	char input_file[MAX_SIZE_FILE_NAME];
	char log_file[MAX_SIZE_FILE_NAME];

	grid p_grid;

} batchtool;

void print_grid(grid* p_grid);
void print_batchtool(batchtool* bt);

void restart_grid(grid* p_grid);
int next_args_config(grid* p_grid, char** new_argv);
int next_inst_config(grid* p_grid, instance* inst);

void free_gridparam(gridparam* params);
void free_grid(grid* p_grid);
void free_batchtool(batchtool* bt);


#endif