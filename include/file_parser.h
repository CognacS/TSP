#ifndef FILE_PARSER_H_  

#define FILE_PARSER_H_

#include "tsp.h"
#include "parser.h"
#include "batch_tool.h"

// SECTIONS FOR TSP FILES
#define TOKEN_SECTION 0
#define COORD_SECTION 1

// SECTIONS FOR BATCH FILES
#define BATCH_GENERAL_SECTION 2
#define BATCH_PARAM_SECTION 3
#define BATCH_VALUES_SECTION 4

// random instances parameters
#define RI_X_MIN 0
#define RI_Y_MIN 0
#define RI_X_SIZE 1000
#define RI_Y_SIZE 1000
#define RI_DIST_TYPE EUC_2D

// *********** definitions for tsp-like file ***********
typedef struct
{
    FILE* fp;
    char line[180];
    char* pos;
} datasource_tsplike;

int next_arg_tsplike(void* ds, char* out_buffer, int section);
int next_value_tsplike(void* ds, char* out_buffer, int section);
// constructor
strings_iterator* build_tsplike_iter(FILE* fin);

void random_instance(Graph* g, char* line);
void read_input(Instance* inst);
void read_batchfile(BatchTool* bt);

#endif
