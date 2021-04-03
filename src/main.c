#include "../main.h"

int main(int argc, char **argv) 
{ 

	if ( argc < 2 ) { printf("Usage: \"%s -help\" for help\n", argv[0]); exit(1); }

	// print arguments if needed
	log_multilines(VERBOSITY, LOGLVL_INFO, argc, argv);

	// initialize start time
	double t1 = second();
	// initialize instance
	instance inst;

	// fill instance with default values
	fill_inst_default(&inst);
	// parse the arguments from the command line and set params of instance
	parse_command_line(argc-1, argv+1, &inst);

	// ***************************** BATCH AUTOMATIZATION ****************************
	// if a batch file was given, automatize runs
	if (strlen(inst.inst_params.batch_file) > 0)
	{
		// initialize batch tool
		batchtool bt;
		strcpy(bt.input_file, inst.inst_params.batch_file);

		read_batchfile(&bt);


		free_batchtool(&bt);
	}
	// *******************************************************************************
	else
		// ***************************** SINGLE RUN MODE ****************************
	{
		// read the input file
		read_input(&inst);

		// print instance if needed
		log_datastruct(&inst, TYPE_INST, VERBOSITY, LOGLVL_DEBUG);

		double* xstar = NULL;

		if (TSPopt(&inst, &xstar)) print_error(NULL, ERR_OPT_PROCEDURE);
		double t2 = second();


		log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] TSP problem solved in %lf sec.s\n", t2 - t1);

		plot_tsp_solution_directed(&inst.inst_graph, xstar);

		// free solution
		free(xstar);
	}
	// free the data structure
	free_instance(&inst);
	
	return 0; 
} 