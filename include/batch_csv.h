#ifndef BATCH_CSV_H_

#define BATCH_CSV_H_

#include "batch_tool.h"

typedef struct
{
	char csv_file[MAX_SIZE_FILE_NAME];
	batchtool bt;
	FILE* csv_fp;

}csv_batchtool;

void open_file_csv(csv_batchtool* bt);
void reorder_grid_csv(grid* p_grid);
void register_time_csv(csv_batchtool* bt, double seconds);

#endif