#include "../log.h"

void log_line(int runlvl, int loglvl, char* line)
{
	// if the running log level is bigger than this log level then print the line
	if (runlvl >= loglvl)
		printf("%s\n", line);

}
void log_multilines(int runlvl, int loglvl, int nlines, char** lines)
{
	// if the running log level is bigger than this log level then print the multilines
	if (runlvl >= loglvl)
		for (int i = 0; i < nlines; i++)
			printf("%s\n", lines[i]);

}