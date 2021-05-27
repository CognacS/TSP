#ifndef BATCH_TOOL_H_

#define BATCH_TOOL_H_

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
} CsvAxis;

typedef struct
{
	char param_name[MAX_SIZE_PARAM_NAME];
	CsvAxis axis;
	int values_num;
	char** labels;
	char** values;

} GridParam;

typedef struct
{
	int started;
	int end_reached;
	int params_num;
	int* indices;
	int last_incridx;
	GridParam* grid_params;
} Grid;

/* ***********************************************************************************
*						BATCH TOOL SECTION
*********************************************************************************** */

typedef struct batch_tool_struct
{

	char input_file[MAX_SIZE_FILE_NAME];
	char log_file[MAX_SIZE_FILE_NAME];

	Grid p_grid;

} BatchTool;

void print_grid(Grid* p_grid);
void print_BatchTool(BatchTool* bt);

void restart_grid(Grid* p_grid);
int next_args_config(Grid* p_grid, char** new_argv);
int next_inst_config(Grid* p_grid, Instance* inst);

void free_GridParam(GridParam* params);
void free_grid(Grid* p_grid);
void free_BatchTool(BatchTool* bt);


#endif