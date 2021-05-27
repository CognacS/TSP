#include "../include/main.h"


int main(int argc, char **argv) 
{

	if ( argc < 2 ) { printf("Usage: \"%s -help\" for help\n", argv[0]); exit(1); }

	// print arguments if needed
	log_multilines(VERBOSITY, LOGLVL_INFO, argc, argv);

	// initialize start time
	double t1 = second();
	// initialize instance
	Instance inst;

	// fill instance with default values
	fill_inst_default(&inst);
	// parse the arguments from the command line and set params of instance
	parse_command_line(argc-1, argv+1, &inst);

	// ***************************** BATCH AUTOMATIZATION ****************************
	// if a batch file was given, automatize runs
	if (strlen(inst.inst_params.batch_file) > 0)
	{

		// setup csv_batchfile
		CsvBatchTool csv_bt;
		csv_bt.reordered = 0;
		strcpy(csv_bt.csv_file, "output.csv");

		// initialize batch tool
		BatchTool* bt = &csv_bt.bt;
		strcpy(bt->input_file, inst.inst_params.batch_file);
		// read batch file
		read_batchfile(bt);
		log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] Parsed batch file: %s", bt->input_file);

		// reorder grid to be used for printing on csv
		reorder_grid_csv(&csv_bt);

		// open csv file
		open_file_csv(&csv_bt);
		log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] Opened output file: %s", csv_bt.csv_file);

		// iterate through each instance
		restart_grid(&bt->p_grid);
		log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] ############ STARTING BATCH PROCEDURE ############ ");
		while (next_inst_config(&bt->p_grid, &inst))
		{
			// log current instance
			log_line(VERBOSITY, LOGLVL_MSG, "[MESSAGE] Now solving:");
			log_datastruct(&inst.inst_params, TYPE_PARAM, VERBOSITY, LOGLVL_MSG);

			// read the input file
			read_input(&inst);

			// solve current instance
			OptResult out_code = TSPopt(&inst);
			double measure = inst.inst_global_data.perf_measure;
			double exec_time = inst.inst_global_data.texec;

			// log execution time
			if (out_code == OPT_OK)
				log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] Solved in %f sec.s", exec_time);
			// log timelimit expired
			if (out_code == OPT_TL_EXPIRED)
				log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] Timelimit of %f expired: ran %f sec.s",
					inst.inst_params.timelimit,  exec_time);
			// log heuristic method's final obj
			if (out_code == OPT_HEUR_OK)
				log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] Heuristic reached an obj of %f in %f sec.s",
					measure, exec_time);

			// register measure (exec time, final obj)
			register_measure_csv(&csv_bt, measure);

			// cleanup current graph
			free_graph(&inst.inst_graph);
		}

		// close csv file
		close_file_csv(&csv_bt);
		// cleanup
		free_BatchTool(bt);
	}
	// *******************************************************************************
	else
	// ******************************* SINGLE RUN MODE *******************************
	{
		// read the input file
		read_input(&inst);

		// print instance if needed
		log_datastruct(&inst, TYPE_INST, VERBOSITY, LOGLVL_DEBUG);

		OptResult out_code = TSPopt(&inst);

		if (out_code == OPT_OK)
		{

			log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] TSP problem solved in %lf sec.s\n", inst.inst_global_data.texec);
			log_datastruct(&inst.inst_global_data, TYPE_GLOB, VERBOSITY, LOGLVL_MSG);

			double* xstar = inst.inst_global_data.xstar;

			switch (model_tsptype(inst.inst_params.model_type))
			{
			case MODEL_TSP_ASYMM:
				if (VERBOSITY >= LOGLVL_INFO) print_directed_sol(&inst.inst_graph, xstar);
				if (VERBOSITY >= LOGLVL_PLOTSOL) plot_tsp_xstar_directed(&inst.inst_graph, xstar);
				break;
			case MODEL_TSP_SYMM:
				if (VERBOSITY >= LOGLVL_INFO) print_undirected_sol(&inst.inst_graph, xstar);
				if (VERBOSITY >= LOGLVL_PLOTSOL) plot_tsp_xstar_undirected(&inst.inst_graph, xstar);
				break;
			}
		}
		if (out_code == OPT_HEUR_OK)
		{
			log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] TSP problem solved heuristically in %lf sec.s\n", inst.inst_global_data.texec);
			log_datastruct(&inst.inst_global_data, TYPE_GLOB, VERBOSITY, LOGLVL_MSG);

			double* xstar = inst.inst_global_data.xstar;
			if (VERBOSITY >= LOGLVL_PLOTSOL) plot_tsp_xstar_undirected(&inst.inst_graph, xstar);
		}

	}
	// *******************************************************************************

	// free the data structure
	free_instance(&inst);

	return 0; 
} 