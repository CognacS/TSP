#ifndef RANDOM_H_

#define RANDOM_H_

#include <stdlib.h>

#define DEFAULT_SEED 123456

void safe_rand(double* rand_num, int a, int b, int c);

// return a uniformly distributed random value in the interval [0,1)
inline double random() { return (double)rand() / (RAND_MAX + 1); }
// return a random integer uniformly distributed in the interval [0,max)
inline int rand_int(int max) { return (int)(random() * max); }

// return a random integer uniformly distributed in the interval [min, max]
inline int rand_int_range(int min, int max) { return rand_int(max - min + 1) + min; }


#endif
