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
		// read batch file
		read_batchfile(&bt);
		// reorder grid to be used for printing on csv
		reorder_grid_csv(&bt.p_grid);

		// iterate through each instance
		restart_grid(&bt.p_grid);
		while (next_inst_config(&bt.p_grid, &inst))
		{
			log_datastruct(&inst.inst_params, TYPE_INST, VERBOSITY, LOGLVL_DEBUG);
		}

		free_batchtool(&bt);
	}
	// *******************************************************************************
	else
	// ******************************* SINGLE RUN MODE *******************************
	{
		// read the input file
		read_input(&inst);

		// print instance if needed
		log_datastruct(&inst, TYPE_INST, VERBOSITY, LOGLVL_DEBUG);

		if (TSPopt(&inst)) print_error(NULL, ERR_OPT_PROCEDURE);
		double t2 = second();


		log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] TSP problem solved in %lf sec.s\n", t2 - t1);
		log_datastruct(&inst.inst_global_data, TYPE_GLOB, VERBOSITY, LOGLVL_MSG);

		double* xstar = inst.inst_global_data.xstar;

		switch (model_tsptype(inst.inst_params.model_type))
		{
		case TSP_ASYMM:
			print_directed_sol(&inst.inst_graph, xstar);
			plot_tsp_solution_directed(&inst.inst_graph, xstar);
			break;
		case TSP_SYMM:
			print_undirected_sol(&inst.inst_graph, xstar);
			plot_tsp_solution_undirected(&inst.inst_graph, xstar);
			break;
		}

	}
	// *******************************************************************************

	// free the data structure
	free_instance(&inst);
	
	return 0; 
} 