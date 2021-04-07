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

void read_batchfile(batchtool* bt)
{
	// open the input file "read-only"
	char input_file[100];
	sprintf(input_file, "./batching/%s", bt->input_file);
	FILE* fin = fopen(input_file, "r");
	// check if the file exists
	if (fin == NULL) print_error(input_file, ERR_INPUT_NOT_EXISTS);

	// set default values
	grid* p_grid = &bt->p_grid;
	p_grid->end_reached = 0;
	p_grid->params_num = -1;

	// general section variables
	int params_counter = 0;
	gridparam* active_param = NULL;

	// param section variables
	int values_counter = 0;
	

	strings_iterator* iter = build_tsplike_iter(fin);

	tokenswitch_sect(iter, BATCH_GENERAL_SECTION)
	{
		section(BATCH_GENERAL_SECTION)
		{
			tokencase_1("PARAMS_NUM", pnum_str)
			{
				if (p_grid->params_num >= 0) print_error("duplicate PARAMS_NUM", ERR_INPUT_FORMAT);
				// set params_num
				p_grid->params_num = atoi(pnum_str);
				// allocate parameters grid
				calloc_s(p_grid->indices, p_grid->params_num, int);
				calloc_s(p_grid->grid_params, p_grid->params_num, gridparam);

				log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] number of batch parameters: %d", p_grid->params_num);

			}
			tokencase_1("PARAM", param_name)
			{
				if (params_counter >= p_grid->params_num) print_error("out of bound for number of batching parameters", ERR_INPUT_FORMAT);
				// setup new parameter
				active_param = &p_grid->grid_params[params_counter++];
				strcpy(active_param->param_name, param_name);
				active_param->values_num = -1;
				values_counter = 0;
				// change section to the specific parameter section
				set_section(BATCH_PARAM_SECTION);

				log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] found parameter \"%s\"", param_name);
			}
		}
		section(BATCH_PARAM_SECTION)
		{
			tokencase_1("CSV_AXIS", axis_str)
			{
				if (!strncmp(axis_str, "ROWS", 4) != 0) { active_param->axis = ROWS; continue; }
				if (!strncmp(axis_str, "COLS", 4) != 0) { active_param->axis = COLS; continue; }
				if (!strncmp(axis_str, "CELL", 4) != 0) { active_param->axis = CELL; continue; }

				log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] csv type axis is: %s", axis_str);
				print_error("only CSV_AXIS == ROWS, COLS possible", ERR_INPUT_FORMAT);
			}
			tokencase_1("VALUES_NUM", vnum_str)
			{
				if (active_param->values_num >= 0) print_error("duplicate VALUES_NUM", ERR_INPUT_FORMAT);
				// set values_num of current parameter
				active_param->values_num = atoi(vnum_str);
				// allocate labels and values of current parameter
				calloc_s(active_param->labels, active_param->values_num, char*);
				calloc_s(active_param->values, active_param->values_num, char*);

				log_line_ext(VERBOSITY, LOGLVL_MSG, "[MESSAGE] number of values for %s : %d", active_param->param_name, active_param->values_num);

			}
			tokencase("VALUES_SECTION")
			{
				set_section(BATCH_VALUES_SECTION);
				if (active_param->values_num <= 0) print_error("VALUES_NUM -> VALUES_SECTION", ERR_PARAM_DISORDERED);
			}
			tokencase("END_PARAM")
			{
				set_section(BATCH_GENERAL_SECTION);

			}
		}
		section(BATCH_VALUES_SECTION)
		{
			tokencase("END_VALUES")
			{
				set_section(BATCH_PARAM_SECTION);
			}
			tokencase_1(ANY_STR, value)
			{
				if (values_counter >= active_param->values_num) print_error("out of bound for number of values in parameter", ERR_INPUT_FORMAT);

				// allocate strings
				calloc_s(active_param->labels[values_counter], strlen(_buffer)+1, char);
				calloc_s(active_param->values[values_counter], strlen(value)+1, char);

				// copy content
				strcpy(active_param->labels[values_counter], _buffer);
				strcpy(active_param->values[values_counter], value);

				log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] value %s : %s", active_param->labels[values_counter], active_param->values[values_counter]);
				values_counter++;
			}
		}
		tokenfinally(param)
		{
			log_line_ext(VERBOSITY, LOGLVL_DEBUG, "[DEBUG] Final active section %d", _section);
			print_warn(param, WARN_IGNORED_INPUT_FILE_PARAM);
		}
	}

	free_iter(iter);

	fclose(fin);

	print_batchtool(bt);

}


int next_arg_tsplike(void* ds, char* out_buffer, int section)
{
	datasource_tsplike* tsplike = (datasource_tsplike*)ds;
	if (fgets(tsplike->line, sizeof(tsplike->line), tsplike->fp) != NULL)
	{
		tsplike->pos = tsplike->line;
		if (strlen(tsplike->line) <= 1) return next_arg_tsplike(ds, out_buffer, section);

		char* token = strsep(&tsplike->pos, " :\t\n");
		strcpy(out_buffer, token);

		return 1;

	}
	return 0;

}

int next_value_tsplike(void* ds, char* out_buffer, int section)
{
	datasource_tsplike* tsplike = (datasource_tsplike*)ds;
	char* token;
	char* sep = NULL;
	int offset = 0;

	switch (section)
	{
	case TOKEN_SECTION:
	case BATCH_GENERAL_SECTION:
	case BATCH_PARAM_SECTION:
	case BATCH_VALUES_SECTION:
		sep = ":\n";
		offset = 1;
		break;
	case COORD_SECTION:
		sep = " \t\n";
		break;
	default:
		return 0;
	}

	if ((token = (strsep(&tsplike->pos, sep) + offset)) != NULL)
	{
		strcpy(out_buffer, token);
		return 1;
	}
	return 0;

}

strings_iterator* build_tsplike_iter(FILE* fin)
{
	datasource_tsplike* ds = (datasource_tsplike*)malloc(sizeof(datasource_tsplike));
	ds->fp = fin;

	strings_iterator* iter = (strings_iterator*)malloc(sizeof(strings_iterator));
	iter->datasource = (void*)ds;
	iter->next_arg = next_arg_tsplike;
	iter->next_value = next_value_tsplike;

	return iter;
}