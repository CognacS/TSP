#ifndef TSP_H_  

#define TSP_H_

#include <cplex.h>
#include "tsp_data.h"
#include "log.h"
#include "error.h"
#include "tsp_utility.h"
#include "cpx_models.h"

// optimization functions
int TSPopt(instance* inst, double** xstar);

//inline
inline int imax(int i1, int i2) { return (i1 > i2) ? i1 : i2; }
inline double dmin(double d1, double d2) { return (d1 < d2) ? d1 : d2; }
inline double dmax(double d1, double d2) { return (d1 > d2) ? d1 : d2; }

#endif   /* TSP_H_ */ 