#include "../comm_parser.h"


void parse_command_line(int argc, char** argv, Instance* inst)
{
	// print command when debugging
	log_line_ext(VERBOSITY, LOGLVL_INFO, "Running command parser with %d parameters", argc);

	// define the parameters for the instance
	Graph* g = &inst->inst_graph;
	Params* p = &inst->inst_params;
    
	int help = 0; if (argc < 1) help = 1;
	strings_iterator* iter = build_cmdline_iter(argc, argv);

	tokenswitch(iter)
	{
		tokencase("-help|--help")
		{
			help = 1;
			break;
		}
		tokencase_1("-model_type|-model", mt_str)
		{
			p->model_type = atoi(mt_str);
		}
		tokencase_1("-time_limit", tl_str)
		{
			p->timelimit = atof(tl_str);
		}
		tokencase_1("-file|-input|-input_file|-f", file_name)
		{
			strcpy(p->input_file, file_name);
		}
		tokencase("-int")
		{
			g->integer_costs = 1;
		}
		tokencase_1("-seed|-random_seed", seed_str)
		{
			p->randomseed = atoi(seed_str);
		}
		tokencase_1("-max_nodes", mn_str)
		{
			p->max_nodes = atoi(mn_str);
		}
		tokencase_1("-cutoff", co_str)
		{
			p->cutoff = atoi(co_str);
		}
		tokencase_1("-batch|-batch_file|-bf", batchfile)
		{
			strcpy(p->batch_file, batchfile);
		}
		tokencase_1("-heur|-heuristic|-heuristics", heur_code)
		{
			strcpy(p->heuristic_code, heur_code);
		}
		tokencase("-ilpheur|-heurilp")
		{
			p->heuristic_ilp = 1;
		}
		tokenfinally(arg)
		{
			print_error(ERR_CLINE_ARG_UNDEF, arg);
		}
	}

	free_iter(iter);

	if (help)
	{
		printf("## HELP TSP SOLVER 06/08/2021 ##\n");
		printf(" -help          \tCall help utility                                 \tSynonims: --help\n");
		printf(" -model <XYZ>   \tModel type to indentify the optimization procedure\tSynonims: -model_type\n");
		printf(" -time_limit <f>\tTime limit for running the optimization process   \n");
		printf(" -file <str>    \tInput instance file name                          \tSynonims: -input, -input_file -f\n");
		printf(" -int           \tSet integer costs for the model                   \n");
		printf(" -seed <i>      \tRandom Seed for CPX RNG                           \n");
		printf(" -max_nodes <i> \tMaximum number of branching nodes in the final run\n");
		printf(" -cutoff <f>    \tCutoff (upper bound) for master                   \n");
		printf(" -batch <str>   \tBatch file name for instance automatization       \tSynonims: -batch_file -bf\n");
		printf(" -heur <str>    \tDefine parameters for heuristics (model=300)      \tSynonims: -heuristic, -heuristics\n");
		printf(" -ilpheur       \tEnable heuristic methods in ILP solvers           \tSynonims: -heurilp\n");

		exit(1);

	}

	// print parameters for information
	log_datastruct(p, TYPE_PARAM, VERBOSITY, LOGLVL_INFO);

}