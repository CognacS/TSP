#ifndef ERROR_H_

#define ERROR_H_

#include <stdio.h>
#include <stdlib.h>
#include "log.h"


#define SUPPRESS_WARNINGS 0

// *************** define errors ***************
// command line errors
#define ERR_ARG_UNDEF			0
#define ERR_NOT_ENOUGH_VALUES	1
// input file errors
#define ERR_INPUT_NOT_EXISTS	2
#define ERR_INPUT_FORMAT		3
#define ERR_PARAM_REPEATED		4
#define ERR_PARAM_DISORDERED	5
#define ERR_UNKNOWN_NODE		6
// main errors
#define ERR_OPT_PROCEDURE		7
// cplex errors
#define ERR_CPLEX				8
// general purpose 
#define ERR_INCORRECT_FUNCTION	9
#define ERR_INPUT_VALUES		10
#define ERR_NO_MEM_FOR_ALLOC	11
// models errors
#define ERR_WRONG_TSP_PROCEDURE	12
#define ERR_MODEL_NOT_IMPL		13
#define ERR_GENERIC_INCONSIST	14
#define ERR_NO_SOLUTION			15
#define ERR_ADD_CUT				16

// *************** define warnings ***************
#define WARN_WRONG_DATASTRUCT			0	// unknown data structure (inst, graph, etc...)
#define WARN_IGNORED_INPUT_FILE_PARAM	1	// unknown parameter of input file to be skipped

static const char* const error_msgs[] =
{
	"argument type not correct: ",
	"not enough values on the command line for the argument, given/required = ",
	"input file does not exists: ",
	"format error: ",
	"repeated section in input file: ",
	"should be ordered: ",
	"unknown node number",
	"error in optimization procedure",
	"CPX error: ",
	"function error: ",
	"input values error: ",
	"NULL pointer returned when allocating array: ",
	"wrong TSP solver procedure used for TSP variant",
	"the requested model was not implemented: ",
	"generic inconsistency: ",
	"no solution available",
	"error encountered in mip_add_cut(): "
};

static const char* const warn_msgs[] =
{
	"unknown data structure for logging",
	"unknown parameter from input file, skipped: "
};

void print_warn(char* addline, int warn_type);
void print_error(char* addline, int error_type);

#define calloc_s(ARRAY, SIZE, TYPE) \
if (!(ARRAY = (TYPE*)calloc((SIZE), sizeof(TYPE)))) print_error(#ARRAY, ERR_NO_MEM_FOR_ALLOC)

#define malloc_s(STRUCT, TYPE) \
if (!(STRUCT = (TYPE*)malloc(sizeof(TYPE)))) print_error(#STRUCT, ERR_NO_MEM_FOR_ALLOC)

#define arr_malloc_s(ARRAY, SIZE, TYPE) \
if (!(ARRAY = (TYPE*)malloc((SIZE) * sizeof(TYPE)))) print_error(#ARRAY, ERR_NO_MEM_FOR_ALLOC)

#endif