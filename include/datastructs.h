#ifndef DATASTRUCTS_H_  

#define DATASTRUCTS_H_

#include <math.h>
#include "error.h"
#include "random.h"

/* ***********************************************************************************
*						LINKED LIST STRUCTURE DEFINITION
*********************************************************************************** */

// ************* linked cell *************
typedef struct Cell
{
	double x;
	double y;
	struct Cell* next;
} Cell;

// ************* linked list *************
typedef struct
{
	Cell* start;
	Cell* end;
} LinkedList;

LinkedList* LL_new();
int LL_is_empty(LinkedList* list);
void LL_add_value(LinkedList* list, double x, double y);
void LL_free(LinkedList* list);

/* ***********************************************************************************
*						SET STRUCTURE DEFINITION
*********************************************************************************** */

// ************* set of nodes *************
/*
* Set structure with the property:
*	nodes[indices[node]] = node
*/
typedef struct
{
	int* nodes;
	int* indices;
	int max_size;
	int curr_size;
	int nnodes;
} SetOfNodes;

SetOfNodes* SETN_new(int max_size, int nnodes);
SetOfNodes* SETN_new_allnodes(int nnodes);

inline char SETN_exists(SetOfNodes* set, int node) { return set->indices[node] >= 0; }
inline char SETN_isfull(SetOfNodes* set) { return set->curr_size >= set->max_size; }
inline char SETN_isempty(SetOfNodes* set) { return set->curr_size == 0; }
inline int SETN_get(SetOfNodes* set, int idx) { return set->nodes[idx]; }	// O(1)
char SETN_add(SetOfNodes* set, int node);									// O(1)
char SETN_remove(SetOfNodes* set, int node);								// O(1)
int SETN_rand_node(SetOfNodes* set);										// O(1)
char SETN_reposition(SetOfNodes* set, int node, int repos_idx);				// O(1)
void SETN_free(SetOfNodes* set);
void SETN_deepcopy(SetOfNodes* dst, SetOfNodes* src);						// O(n)
	

/* ***********************************************************************************
*					ORDERED INDEXED ARRAY STRUCTURE DEFINITION
*********************************************************************************** */

// ************* ordered indexed array *************
#define MAX_INDEX_SIZE 3
typedef struct
{
	int arr[MAX_INDEX_SIZE];
} Index;

typedef struct
{
	Index index;
	double value;
} IndexedValue;

// Ordered Indexed array has elements ordered from bigger to smaller, i.e.: oiarr[0] has the biggest value
inline char OIA_eligible(IndexedValue* oiarr, double newvalue) { return oiarr[0].value > newvalue; }
char OIA_insert(IndexedValue* oiarr, int size, IndexedValue newelem);
void OIA_clear(IndexedValue* oiarr, int size);
IndexedValue OIA_pack1(int index, double value);
IndexedValue OIA_pack3(int a, int b, int c, double value);
inline IndexedValue OIA_best(IndexedValue* oiarr, int size) { return oiarr[size - 1]; }
IndexedValue OIA_choose(IndexedValue* oiarr, int size, char includebest);

void IV_deepcopy(IndexedValue* dst, IndexedValue* src);

/* ***********************************************************************************
*					TABU LIST STRUCTURE DEFINITION
*********************************************************************************** */
typedef struct
{
	int* tabu_nodes;
	int tenure;
	int now;
	int nnodes;

} TabuList;

TabuList* TABU_new(int nnodes, int tenure);
char TABU_istabu(TabuList* tabu, int node);
inline void TABU_update(TabuList* tabu, int node) { tabu->tabu_nodes[node] = tabu->now; }
inline void TABU_advance(TabuList* tabu) { tabu->now++; }
inline void TABU_set_tenure(TabuList* tabu, int newtenure) { tabu->tenure = newtenure; }
void TABU_free(TabuList* tabu);
void TABU_print(TabuList* tabu);

inline char using_tabu(TabuList* tabu) { return tabu != NULL; }

#endif