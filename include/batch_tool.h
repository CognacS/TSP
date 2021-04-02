#ifndef BATCH_TOOL_H_

#define BATCH_TOOL_H_

#include "tsp.h"

typedef enum
{
	MODEL_TYPE,
	INPUT_FILE
} gridparamtype;

typedef struct
{
	char param_name[100];
	gridparamtype type;
	int values_num;
	char** labels;
	void** values;

} gridparam;

typedef struct
{
	int params_num;
	int* coords;
	gridparam* params;
	int* out_file_is_row;
} grid;

typedef struct
{

	char* input_file;
	char* log_file;

	grid p_grid;

} batchtool;


void tostring_batchtool(batchtool* bt);
int next(grid* p_grid, instance* inst);

void free_gridparam(gridparam* params);
void free_grid(grid* p_grid);
void free_batchtool(batchtool* bt);

#endif