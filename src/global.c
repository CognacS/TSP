#include "../include/global.h"

int nearest_integer(double value)
{
	int int_dist = (int)(value + 0.5 - ROUNDING_EPSILON);
	return int_dist;
}

int is_one(double num)
{
	return num > 0.5;
}

int is_zero_strict(double num)
{
	return num < ROUNDING_EPSILON;
}

int is_one_strict(double num)
{
	return num > 1.0 - ROUNDING_EPSILON;
}