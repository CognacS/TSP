#ifndef PARSER_H_

#define PARSER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#define MAX_SIZE 100

typedef struct
{
    void* datasource;
    int (*next_arg)(void*, char*, int);
    int (*next_value)(void*, char*, int);

} strings_iterator;

// utility functions to assign in strings_iterator
// *********** definitions for command line ***********
typedef struct
{
    int argc;
    char** argv;
    int pos;
} datasource_cmdline;

int next_arg_cmdline(void* ds, char* out_buffer, int section);
int next_value_cmdline(void* ds, char* out_buffer, int section);
// constructor
strings_iterator* build_cmdline_iter(int argc, char** argv);

// *********** definitions for tsp-like file ***********
#define TOKEN_SECTION 0
#define COORD_SECTION 1

typedef struct
{
    FILE* fp;
    char line[180];
    char* pos;
} datasource_tsplike;

int next_arg_tsplike(void* ds, char* out_buffer, int section);
int next_value_tsplike(void* ds, char* out_buffer, int section);
// constructor
strings_iterator* build_tsplike_iter(FILE* fin);


// common functions
void free_iter(strings_iterator* iter);

// parse an input string with a simple regex with OR options divided by "|"
int parsestr(char* in_str, char* format);


#define ANY_STR ""

#define tokenswitch_sect(ITERATOR, STARTING_SECTION) \
int _section = (STARTING_SECTION);\
strings_iterator* _iterator = ITERATOR;\
char _buffer[MAX_SIZE];\
for(int _match = 0; _iterator->next_arg(_iterator->datasource, _buffer, _section); _match = 0){

#define tokenswitch(ITERATOR) \
tokenswitch_sect(ITERATOR, 0)

#define set_section(SECTION) \
_section = (SECTION)

#define section(SECTION) \
    }if (_section == (SECTION)){

#define tokencase_n(COMM_REGEX, NUM, CONTAINER)\
}{\
    if (_match) continue;\
    char _format[] = COMM_REGEX;\
    _match = parsestr(_buffer, _format);\
    char CONTAINER[NUM][MAX_SIZE];\
    char _val_buffer[MAX_SIZE];\
    if (_match)\
    {\
        for (int _j = 0; _j < (NUM); _j++)\
        {\
            _iterator->next_value(_iterator->datasource, _val_buffer, _section);\
            strcpy(CONTAINER[_j], _val_buffer);\
        }\
    }\
    if (_match)

#define tokencase_1(COMM_REGEX, CONTAINER)\
}{\
    if (_match) continue;\
    char _format[] = COMM_REGEX;\
    _match = parsestr(_buffer, _format);\
    char CONTAINER[MAX_SIZE];\
    char _val_buffer[MAX_SIZE];\
    if (_match)\
    {\
        _iterator->next_value(_iterator->datasource, _val_buffer, _section);\
        strcpy(CONTAINER, _val_buffer);\
    }\
    if (_match)

#define tokencase(COMM_REGEX)\
}{\
    if (_match) continue;\
    char _format[] = COMM_REGEX;\
    _match = parsestr(_buffer, _format);\
    if (_match)


#define tokenfinally(ARG) \
}\
if (_match) continue;\
char ARG[MAX_SIZE];\
strcpy(ARG, _buffer);

/*
#define tokenswitch(SIZE, STRINGS) \
char** _str_to_parse = (char**)STRINGS;\
int _size = SIZE;\
for (int _i = 0; _i < _size; _i++)\
{ int match = 0;{

#define tokencase_n(COMM_REGEX, NUM, CONTAINER)\
}}{\
    if (match) continue;\
    char format[] = COMM_REGEX;\
    match = parsestr(_str_to_parse[_i], format);\
    if (match)\
    {\
        char CONTAINER[NUM][MAX_SIZE];\
        for (int _j = 0; _j < NUM; _j++)\
        {\
            strcpy(CONTAINER[_j], _str_to_parse[++_i]);\
        }

#define tokencase_1(COMM_REGEX, CONTAINER)\
}}{\
    if (match) continue;\
    char format[] = COMM_REGEX;\
    match = parsestr(_str_to_parse[_i], format);\
    if (match)\
    {\
        char CONTAINER[MAX_SIZE];\
        strcpy(CONTAINER, _str_to_parse[++_i]);\

#define tokencase(COMM_REGEX)\
}}{\
    if (match) continue;\
    char format[] = COMM_REGEX;\
    match = parsestr(_str_to_parse[_i], format);\
    if (match)\
    {
        

#define tokenfinally(ARG) \
}}\
if (match) continue;\
char ARG[MAX_SIZE];\
strcpy(ARG, _str_to_parse[_i]);\

*/

#endif // !PARSER_H_
