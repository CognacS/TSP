#ifndef OPT_UTILITY_H_  

#define OPT_UTILITY_H_

#include "tsp.h"

// 2-opt move
double move_2opt(int* succ, Graph* g, char allow_unimproving);
double remove_crossings(int* succ, Graph* g);

#endif