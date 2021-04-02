#include "../parser.h"

int next_arg_cmdline(void* ds, char* out_buffer, int section)
{
	datasource_cmdline* cmdline = (datasource_cmdline*)ds;
	if (cmdline->pos < cmdline->argc)
	{
		strcpy(out_buffer, cmdline->argv[cmdline->pos++]);
		return 1;
	}
	else {
		return 0;
	}
}

int next_value_cmdline(void* ds, char* out_buffer, int section)
{
	return next_arg_cmdline(ds, out_buffer, section);
}


strings_iterator* build_cmdline_iter(int argc, char** argv)
{
	datasource_cmdline* ds = (datasource_cmdline*)malloc(sizeof(datasource_cmdline));
	ds->argc = argc;
	ds->argv = argv;
	ds->pos = 0;

	strings_iterator* iter = (strings_iterator*)malloc(sizeof(strings_iterator));
	iter->datasource = (void*)ds;
	iter->next_arg = next_arg_cmdline;
	iter->next_value = next_value_cmdline;

	return iter;
}

char* strsep(char** stringp, const char* delim)
{
	if (*stringp == NULL) return NULL;

	size_t len = strlen(*stringp);
	if (len == 0)
	{
		*stringp = NULL;
		return NULL;
	}

	char* token = strtok(*stringp, delim);
	size_t toklen = strlen(token);

	if (len == toklen)
	{
		*stringp = NULL;
	}
	else
	{
		*stringp = token + toklen + 1;
	}
	return token;
}

int next_arg_tsplike(void* ds, char* out_buffer, int section)
{
	datasource_tsplike* tsplike = (datasource_tsplike*)ds;
	if (fgets(tsplike->line, sizeof(tsplike->line), tsplike->fp) != NULL)
	{
		tsplike->pos = tsplike->line;
		if (strlen(tsplike->line) <= 1) return next_arg_tsplike(ds, out_buffer, section);


		strcpy(out_buffer, strsep(&tsplike->pos, " :\n"));
		
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


void free_iter(strings_iterator* iter)
{
	free(iter->datasource);
	free(iter);
}

int parsestr(char* in_str, char* format)
{
	if (!format[0]) return 1;

	// left trim
	while (*in_str == ' ') in_str++;
	// right trim
	char* last = in_str + strlen(in_str) - 1;
	while (*last == ' ') last--;
	*(last + 1) = 0;

	char* token = strtok(format, "|");

	// while there are new strings to parse
	while (token != NULL)
	{
		// if it matches, return true
		if (!strcmp(in_str, token))
			return 1;

		token = strtok(NULL, "|");
	}
	// if in_str did not match with anything return false
	return 0;
}