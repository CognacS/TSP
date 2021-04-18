#ifndef RANDOM_H_

#define RANDOM_H_


#ifdef __linux__

#include <stdlib.h>
#include "error.h"

#define safe_rand(rvalue_ptr) \
	do\
	{\
		unsigned int rand_state = time(NULL) ^ getpid() ^ pthread_self();\
		int rand_int = rand_r(&rand_state)); \
		*(rvalue_ptr) = (double)rand_int / ((double)RAND_MAX + 1); \
	}while(0)


#elif _WIN32

#define _CRT_RAND_S
#include <stdlib.h>
#include "error.h"

#define safe_rand(rvalue_ptr) \
	do\
	{\
		errno_t err;\
		unsigned int rand_int = 0;\
		if (err = rand_s(&rand_int)) \
			print_error_ext(ERR_RNG_ERR, "errno: %d", err); \
		*(rvalue_ptr) = (double)rand_int / ((double)UINT_MAX + 1); \
	}while(0)

#endif



#endif
