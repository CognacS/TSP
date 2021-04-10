#ifndef FILE_PARSER_H_  

#define FILE_PARSER_H_

#include "tsp_data.h"
#include "tsp_utility.h"
#include "log.h"
#include "error.h"
#include "parser.h"
#include "batch_tool.h"

// SECTIONS FOR TSP FILES
#define TOKEN_SECTION 0
#define COORD_SECTION 1

// SECTIONS FOR BATCH FILES
#define BATCH_GENERAL_SECTION 2
#define BATCH_PARAM_SECTION 3
#define BATCH_VALUES_SECTION 4

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


void read_input(instance* inst);
void read_batchfile(batchtool* bt);

#endif
