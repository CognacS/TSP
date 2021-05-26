#ifndef TOKENS_H_  

#define TOKENS_H_

#include <string.h>

#ifdef _WIN32
#define strtok_u strtok_s

#else
#define strotok_u strtok_r

#endif

#endif