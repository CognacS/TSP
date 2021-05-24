#ifndef RANDOM_H_

#define RANDOM_H_

#include <stdlib.h>

#define DEFAULT_SEED 123456

void safe_rand(double* rand_num, int a, int b, int c);

// return a uniformly distributed random value in the interval [0,1)
inline double random() { return (double)rand() / (RAND_MAX + 1); }


#endif
