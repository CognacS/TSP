#include "../error.h"

void print_warn(warn_code warn_type, char* addline)
{
	if (!SUPPRESS_WARNINGS)
	{
		// search for the right warn OR wait for warn unknown
		int i = 0;
		while (	warn_msgs[i].num != warn_type &&
				warn_msgs[i].num != WARN_UNKNOWN_WARN) i++;
		// when found, log it
		if (addline == NULL) addline = "";

		char line[1000];
		// concat error line and addline
		sprintf(line, "[WARN] %s %s\n", warn_msgs[i].str, addline);

		// log warn
		log_line(VERBOSITY, LOGLVL_WARN, line);
	}
}
void print_error(err_code error_type, char* addline)
{
	// search for the right error OR wait for error unknown
	int i = 0;
	while( error_msgs[i].num != error_type &&
		   error_msgs[i].num != ERR_UNKNOWN_ERR) i++;
	// when found, log it and exit
	if (addline == NULL) addline = "";

	char line[1000];
	// concat error line and addline
	sprintf(line, "[ERROR] %s %s\n", error_msgs[i].str, addline);
	
	// log error
	log_line(VERBOSITY, LOGLVL_ERROR, line);
	fflush(NULL);
	exit(1);
}