#ifndef GLOBAL_H_

#define GLOBAL_H_

// ********************************** DEFINITIONS ************************************

typedef enum
{
	NO = 0,
	YES = 1
} answer;

// ******************************** VALUE INFERENCE **********************************
#define ROUNDING_EPSILON	1e-5	// tolerance for rounding integer values
// value rounding to nearest integer
int nearest_integer(double value);

// variable inference of 0/1 value functions
int is_one(double num);
int is_zero_strict(double num);
int is_one_strict(double num);

#endif