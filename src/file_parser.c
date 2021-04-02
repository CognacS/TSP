#include "../file_parser.h"

void read_input(instance* inst)
{
	// open the input file "read-only"
	char input_file[100];
	sprintf(input_file, "./data/%s", inst->inst_params.input_file);
	FILE* fin = fopen(input_file, "r");
	// check if the file exists
	if (fin == NULL) print_error(input_file, ERR_INPUT_NOT_EXISTS);

	// get the instance graph
	graph* g = &inst->inst_graph;
	g->nnodes = -1;
	g->tr_xcoord = NULL;
	g->tr_ycoord = NULL;
	// get the instance params
	params* p = &inst->inst_params;


	strings_iterator* iter = build_tsplike_iter(fin);

	tokenswitch_sect(iter, TOKEN_SECTION)
	{
		section(TOKEN_SECTION)
		{
			tokencase_1("NAME", name)
			{
				log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] instance name: \"%s\"", name);

			}
			tokencase_1("COMMENT", comment)
			{
				log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] solving instance \"%s\" with model %d", comment, p->model_type);
			}
			tokencase_1("TYPE", type)
			{
				if (strncmp(type, "TSP", 4) != 0) print_error("only TYPE == TSP implemented so far", ERR_INPUT_FORMAT);
			}
			tokencase_1("DIMENSION", dim)
			{
				g->nnodes = atoi(dim);
				log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] nnodes = %d", g->nnodes);
				calloc_s(g->xcoord, g->nnodes, double);
				calloc_s(g->ycoord, g->nnodes, double);

			}
			tokencase_1("EDGE_WEIGHT_TYPE", ewt_str)
			{
				if (!strncmp(ewt_str, "ATT", 3))	{ g->distance_type = ATT; continue; }
				if (!strncmp(ewt_str, "EUC_2D", 6)) { g->distance_type = EUC_2D; continue; }
				if (!strncmp(ewt_str, "GEO", 3))	{ g->distance_type = GEO; continue; }

				print_error("only EDGE_WEIGHT_TYPE == ATT, EUC_2D, GEO implemented so far", ERR_INPUT_FORMAT);

			}
			tokencase("NODE_COORD_SECTION")
			{
				set_section(COORD_SECTION);
				if (g->nnodes <= 0) print_error("DIMENSION -> NODE_COORD_SECTION", ERR_PARAM_DISORDERED);
			}
		}
		section(COORD_SECTION)
		{
			tokencase("EOF")
			{
				set_section(COORD_SECTION);
			}
			tokencase_n(ANY_STR, 2, coords_str)
			{
				int i = atoi(_buffer) - 1;
				if (i < 0 || i >= g->nnodes) print_error(NULL, ERR_UNKNOWN_NODE);
				g->xcoord[i] = atof(coords_str[0]);
				g->ycoord[i] = atof(coords_str[1]);
				log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] node %4d at coordinates ( %15.7lf , %15.7lf )", i + 1, g->xcoord[i], g->ycoord[i]);
			}
		}
		tokenfinally(param)
		{
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Final active section %d", _section);
			print_warn(param, WARN_IGNORED_INPUT_FILE_PARAM);
		}
	}
	
	free_iter(iter);

	// print graph for information
	log_datastruct(g, TYPE_GRAPH, VERBOSITY, LOGLVL_DEBUG);

	fclose(fin);

}