#ifndef ERROR_H_

#define ERROR_H_

#include <stdio.h>
#include <stdlib.h>
#include "log.h"

#define SUPPRESS_WARNINGS 0

// *************** define errors ***************
typedef enum
{
	// command line errors
	ERR_CLINE_ARG_UNDEF =			0000,
	ERR_CLINE_NOT_ENOUGH_VALUES =	0001,
	// TSP input file errors
	ERR_INPUT_NOT_EXISTS =			0100,
	ERR_INPUT_FORMAT =				0101,
	ERR_INPUT_PARAM_REPEATED = 		0102,
	ERR_INPUT_PARAM_DISORDERED =	0103,
	ERR_INPUT_UNKNOWN_NODE = 		0104,
	// optimization errors
	ERR_OPT_PROCEDURE = 			0200,
	// cplex errors
	ERR_CPLEX = 					0300,
	// general purpose 
	ERR_INCORRECT_FUNCTION_IMPL =	0400,
	ERR_INVALID_FUNC_ARGS = 		0401,
	ERR_NO_MEM_FOR_ALLOC = 			0402,
	ERR_UNKNOWN_ERR =				0403,
	// models errors
	ERR_WRONG_TSP_PROCEDURE = 		0500,
	ERR_MODEL_NOT_IMPL = 			0501,
	ERR_MODEL_INCONSISTENT = 		0502,
	ERR_NO_SOLUTION = 				0503,
	ERR_ADD_CUT = 					0504,
	// batchtool errors
	ERR_CSV_CANNOT_OPEN =			0600,
	ERR_CSV_NOT_REORDERED =			0601
} err_code;


// *************** define warnings ***************
typedef enum
{
	WARN_WRONG_DATASTRUCT =			0000,	// unknown data structure (inst, graph, etc...)
	WARN_IGNORED_INPUT_FILE_PARAM = 0100,	// unknown parameter of input file to be skipped
	// optimization warnings
	WARN_EXPIRED_TIMELIMIT =		0200,
	// general purpose
	WARN_UNKNOWN_WARN =				0403
} warn_code;

typedef struct
{
	const int num;
	const char* str;

} numbered_str;

static const numbered_str error_msgs[] =
{
	{ERR_CLINE_ARG_UNDEF,			"argument type not correct:"},
	{ERR_CLINE_NOT_ENOUGH_VALUES,	"not enough values on the command line for the argument, given/required ="},
	{ERR_INPUT_NOT_EXISTS,			"input file does not exists:"},
	{ERR_INPUT_FORMAT,				"format error:"},
	{ERR_INPUT_PARAM_REPEATED,		"repeated section in input file:"},
	{ERR_INPUT_PARAM_DISORDERED,	"should be ordered:"},
	{ERR_INPUT_UNKNOWN_NODE,		"unknown node number"},
	{ERR_OPT_PROCEDURE,				"error in optimization procedure"},
	{ERR_CPLEX,						"CPX error:"},
	{ERR_INCORRECT_FUNCTION_IMPL,	"function error:"},
	{ERR_INVALID_FUNC_ARGS,			"invalid input arguments in function: "},
	{ERR_NO_MEM_FOR_ALLOC,			"NULL pointer returned when allocating array:"},
	{ERR_WRONG_TSP_PROCEDURE,		"wrong TSP solver procedure used for TSP variant"},
	{ERR_MODEL_NOT_IMPL,			"the requested model was not implemented:"},
	{ERR_MODEL_INCONSISTENT,		"generic inconsistency: "},
	{ERR_NO_SOLUTION,				"no solution available"},
	{ERR_ADD_CUT,					"error encountered in mip_add_cut():"},
	{ERR_CSV_CANNOT_OPEN,			"could not open output CSV file:"},
	{ERR_CSV_NOT_REORDERED,			"grid not reordered before opening csv file"},
	{ERR_UNKNOWN_ERR,				"unknown error:"}
};

static const numbered_str warn_msgs[] =
{
	{WARN_WRONG_DATASTRUCT,			"unknown data structure for logging"},
	{WARN_IGNORED_INPUT_FILE_PARAM,	"unknown parameter from input file, skipped:"},
	{WARN_EXPIRED_TIMELIMIT,		"time limit expired:"},
	{WARN_UNKNOWN_WARN,				"unknown warn:"}
};


void print_warn(warn_code warn_type, char* addline);
void print_error(err_code error_type, char* addline);

#define print_warn_ext(warn_type, format, ...) \
{\
char line[MAX_LOG_LINE_SIZE];\
sprintf(line, format, ##__VA_ARGS__);\
print_warn(warn_type, line);\
}

#define print_error_ext(error_type, format, ...) \
{\
char line[MAX_LOG_LINE_SIZE];\
sprintf(line, format, ##__VA_ARGS__);\
print_error(error_type, line);\
}

#define calloc_s(ARRAY, SIZE, TYPE) \
if (!(ARRAY = (TYPE*)calloc((SIZE), sizeof(TYPE)))) print_error(ERR_NO_MEM_FOR_ALLOC, #ARRAY)

#define malloc_s(STRUCT, TYPE) \
if (!(STRUCT = (TYPE*)malloc(sizeof(TYPE)))) print_error(ERR_NO_MEM_FOR_ALLOC, #STRUCT)

#define arr_malloc_s(ARRAY, SIZE, TYPE) \
if (!(ARRAY = (TYPE*)malloc((SIZE) * sizeof(TYPE)))) print_error(ERR_NO_MEM_FOR_ALLOC, #ARRAY)

#endif