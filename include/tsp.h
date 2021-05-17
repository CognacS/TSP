#ifndef TSP_H_  

#define TSP_H_

#include "optimizer_data.h"
#include "log.h"
#include "error.h"
#include "tsp_utility.h"
#include "heuristics.h"
#include "chrono.h"

#define OPT_TL_PENALIZATION_MULT 10



// optimization functions
int TSPopt(instance* inst);


#endif   /* TSP_H_ */ 