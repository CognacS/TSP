#ifndef BATCH_CSV_H_

#define BATCH_CSV_H_

#include "log.h"
#include "error.h"
#include "batch_tool.h"

typedef struct
{
	char csv_file[MAX_SIZE_FILE_NAME];
	batchtool bt;
	FILE* csv_fp;
	int reordered;

}csv_batchtool;

void open_file_csv(csv_batchtool* bt);
void close_file_csv(csv_batchtool* bt);
void reorder_grid_csv(csv_batchtool* bt);
void register_measure_csv(csv_batchtool* bt, double measure);
void newline_csv(csv_batchtool* bt);

#endif