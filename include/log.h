#ifndef LOG_H_  

#define LOG_H_

#include <stdio.h>

#define MAX_LOG_LINE_SIZE 100

// define log levels
#define LOGLVL_ERROR		10
#define LOGLVL_WARN			20
#define LOGLVL_MSG			40
#define LOGLVL_PLOTSOL		45
#define LOGLVL_INFO			50
#define LOGLVL_EXTINFO		60
#define LOGLVL_DEBUG		70
#define LOGLVL_DEBUGPLOT	75
#define LOGLVL_DEBUGPLOT_P	76
#define LOGLVL_DEBUGPLOT_PP	77
#define LOGLVL_PEDANTIC		90
#define LOGLVL_CPLEXLOG		100

#define VERBOSITY			50		// printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)

void log_line(int runlvl, int loglvl, char* line);
void log_multilines(int runlvl, int loglvl, int nlines, char** lines);

#define log_line_ext(runlvl, loglvl, format, ...) \
if (runlvl >= loglvl) printf(format "\n", __VA_ARGS__);



#endif