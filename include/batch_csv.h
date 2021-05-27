#ifndef BATCH_CSV_H_

#define BATCH_CSV_H_

#include "batch_tool.h"

typedef struct
{
	char csv_file[MAX_SIZE_FILE_NAME];
	BatchTool bt;
	FILE* csv_fp;
	int reordered;

}CsvBatchTool;

void open_file_csv(CsvBatchTool* bt);
void close_file_csv(CsvBatchTool* bt);
void reorder_grid_csv(CsvBatchTool* bt);
void register_measure_csv(CsvBatchTool* bt, double measure);
void newline_csv(CsvBatchTool* bt);

#endif