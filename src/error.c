#include "../error.h"

void print_warn(char* addline, int warn_type)
{
	if (!SUPPRESS_WARNINGS)
	{
		char line[1000];
		// concat warn line and addline
		sprintf(line, "WARN: %s %s\n", warn_msgs[warn_type], addline);
		// log warn
		log_line(VERBOSITY, LOGLVL_WARN, line);
	}
}
void print_error(char* addline, int error_type)
{
	char line[1000];
	// concat error line and addline
	sprintf(line, "ERROR: %s %s\n", error_msgs[error_type], addline);
	// log error
	log_line(VERBOSITY, LOGLVL_ERROR, line);
	fflush(NULL);
	exit(1);
}