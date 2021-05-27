#ifndef ERROR_H_

#define ERROR_H_

#include <stdlib.h>
#include "log.h"

#define SUPPRESS_WARNINGS 0

// *************** define errors ***************
typedef enum
{
	// command line errors
	ERR_CLINE_ARG_UNDEF =			000,
	ERR_CLINE_NOT_ENOUGH_VALUES =	001,
	// TSP input file errors
	ERR_INPUT_NOT_EXISTS =			100,
	ERR_INPUT_FORMAT =				101,
	ERR_INPUT_PARAM_REPEATED = 		102,
	ERR_INPUT_PARAM_DISORDERED =	103,
	ERR_INPUT_UNKNOWN_NODE = 		104,
	// optimization errors
	ERR_OPT_PROCEDURE = 			200,
	ERR_OPT_CPLEX_REDEFINITION =	201,
	ERR_OPT_SOLPTR_INCONSISTENT =	202,
	ERR_OPT_SOLFORMAT_NOTSET	=	203,
	// cplex errors
	ERR_CPLEX = 					300,
	// general purpose 
	ERR_INCORRECT_FUNCTION_IMPL =	400,
	ERR_INVALID_FUNC_ARGS = 		401,
	ERR_NO_MEM_FOR_ALLOC = 			402,
	ERR_UNKNOWN_ERR =				403,
	ERR_RNG_ERR =					404,
	ERR_GENERIC_NOT_IMPL =			405,
	// models errors
	ERR_WRONG_TSP_PROCEDURE = 		500,
	ERR_MODEL_NOT_IMPL = 			501,
	ERR_MODEL_INCONSISTENT = 		502,
	ERR_NO_SOLUTION = 				503,
	ERR_ADD_CUT = 					504,
	ERR_CB_UNDEF_PROCEDURE =		505,
	// BatchTool errors
	ERR_CSV_CANNOT_OPEN =			600,
	ERR_CSV_NOT_REORDERED =			601,
	// concorde errors
	ERR_CONCORDE =					700,
	// heuristic errors
	ERR_HEUR_DECODE =				800,
	ERR_HEUR_CONSTR_PARAM	=		801,
	ERR_HEUR_REFINE_PARAM	=		802,
	ERR_HEUR_INFEASABLE_SOL	=		803

} err_code;


// *************** define warnings ***************
typedef enum
{
	WARN_WRONG_DATASTRUCT =			000,	// unknown data structure (inst, graph, etc...)
	WARN_IGNORED_INPUT_FILE_PARAM = 100,	// unknown parameter of input file to be skipped
	// optimization warnings
	WARN_EXPIRED_TIMELIMIT =		200,
	// general purpose
	WARN_UNKNOWN_WARN =				403
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
	{ERR_OPT_CPLEX_REDEFINITION,	"trying to redefine cplex environment, close it first"},
	{ERR_OPT_SOLPTR_INCONSISTENT,	"pointer is invalid for solution format:"},
	{ERR_OPT_SOLFORMAT_NOTSET,		"solution format not set (no solution available?)"},
	{ERR_CPLEX,						"CPX error:"},
	{ERR_INCORRECT_FUNCTION_IMPL,	"function error:"},
	{ERR_INVALID_FUNC_ARGS,			"invalid input arguments in function: "},
	{ERR_NO_MEM_FOR_ALLOC,			"NULL pointer returned when allocating array:"},
	{ERR_RNG_ERR,					"the Random Number Generator procedure returned an error"},
	{ERR_GENERIC_NOT_IMPL,			"feature not implemented:"},
	{ERR_WRONG_TSP_PROCEDURE,		"wrong TSP solver procedure used for TSP variant"},
	{ERR_MODEL_NOT_IMPL,			"the requested model was not implemented:"},
	{ERR_MODEL_INCONSISTENT,		"generic inconsistency: "},
	{ERR_NO_SOLUTION,				"no solution available"},
	{ERR_ADD_CUT,					"error encountered in cpx_add_cut():"},
	{ERR_CB_UNDEF_PROCEDURE,		"undefined procedure in callback, check contextid for:"},
	{ERR_CSV_CANNOT_OPEN,			"could not open output CSV file:"},
	{ERR_CSV_NOT_REORDERED,			"grid not reordered before opening csv file"},
	{ERR_CONCORDE,					"Concorde error:"},
	{ERR_HEUR_DECODE,				"Heuristic decoding error:"},
	{ERR_HEUR_CONSTR_PARAM,			"Constructive heuristic got wrong parameter"},
	{ERR_HEUR_REFINE_PARAM,			"Refining heuristic got wrong parameter"},
	{ERR_HEUR_INFEASABLE_SOL,		"Got an infeasable solution during the procedure"},
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
do{\
if (!(ARRAY = (TYPE*)calloc((SIZE), sizeof(TYPE)))) print_error(ERR_NO_MEM_FOR_ALLOC, #ARRAY);\
}while(0)

#define malloc_s(STRUCT, TYPE) \
do{\
if (!(STRUCT = (TYPE*)malloc(sizeof(TYPE)))) print_error(ERR_NO_MEM_FOR_ALLOC, #STRUCT);\
}while (0)

#define arr_malloc_s(ARRAY, SIZE, TYPE) \
do{\
if (!(ARRAY = (TYPE*)malloc((SIZE) * sizeof(TYPE)))) print_error(ERR_NO_MEM_FOR_ALLOC, #ARRAY);\
}while(0)

#define free_s(POINTER) \
do{\
free(POINTER); POINTER = NULL;\
}while(0)

#endif