#include "../include/datastructs.h"

/* **************************************************************************************************
*							LINKED LIST
************************************************************************************************** */
LinkedList* LL_new()
{
	LinkedList* ll;	malloc_s(ll, LinkedList);
	ll->start = NULL; ll->end = NULL;
	return ll;
}

int LL_is_empty(LinkedList* list)
{
	return list->start == NULL;
}

void LL_add_value(LinkedList* list, double x, double y)
{
	Cell* cell = NULL;	malloc_s(cell, Cell);
	cell->x = x;
	cell->y = y;
	cell->next = NULL;
	// if list is not initialized
	if (list->end == NULL) list->end = list->start = cell;
	// else add at the end
	else {
		list->end->next = cell;
		list->end = cell;
	}
}

void LL_free(LinkedList* list)
{
	if (list->start != NULL)
	{
		// iterate through all the list and free
		for (Cell* cell = list->start, *next_cell = list->start->next;
			next_cell != NULL;
			cell = next_cell, next_cell = next_cell->next)
		{
			free(cell);
		}
		free(list->end);
	}
	free(list);
}

/* **************************************************************************************************
*							SET OF NODES
************************************************************************************************** */
SetOfNodes* SETN_new(int max_size, int nnodes)
{
	SetOfNodes* set = NULL;		malloc_s(set, SetOfNodes);
	arr_malloc_s(set->nodes, max_size, int);
	arr_malloc_s(set->indices, nnodes, int);
	for (int i = 0; i < nnodes; i++) set->indices[i] = -1;
	set->max_size = max_size;
	set->curr_size = 0;
	set->nnodes = nnodes;
	return set;
}

SetOfNodes* SETN_new_allnodes(int nnodes)
{
	SetOfNodes* set = SETN_new(nnodes, nnodes);
	// fill with all nodes
	for (int i = 0; i < nnodes; i++)
	{
		set->nodes[i] = i;
		set->indices[i] = i;
	}
	set->curr_size = nnodes;
	return set;
}

char SETN_add(SetOfNodes* set, int node)
{
	// if set is full or node already in set, return false
	if (SETN_isfull(set) || SETN_exists(set, node)) return 0;

	// else add node at end
	set->nodes[set->curr_size] = node;
	set->indices[node] = set->curr_size++;
	return 1;
}
char SETN_remove(SetOfNodes* set, int node)
{
	if (SETN_exists(set, node))
	{
		// get removed node index
		int node_idx = set->indices[node];
		// get last node
		int last_node = set->nodes[set->curr_size-1];
		// replace node with last node
		set->nodes[node_idx] = last_node;
		// set last node index as the previous removed node index
		set->indices[last_node] = node_idx;
		set->indices[node] = -1;
		set->curr_size--;
		return 1;
	}
	return 0;
}

int SETN_rand_node(SetOfNodes* set)
{
	int idx = rand_int(set->curr_size);
	return set->nodes[idx];
}

char SETN_reposition(SetOfNodes* set, int node, int repos_idx)
{
	if (SETN_exists(set, node))
	{
		// place node
		int repos_node = set->nodes[repos_idx];
		set->nodes[repos_idx] = node;
		// change node index
		int node_idx = set->indices[node];
		set->indices[node] = repos_idx;
		// swap previous node
		set->nodes[node_idx] = repos_node;
		set->indices[repos_node] = node_idx;

		return 1;
	}
	return 0;
}

void SETN_free(SetOfNodes* set)
{
	free(set->nodes);
	free(set->indices);
	free(set);
}

void SETN_deepcopy(SetOfNodes* dst, SetOfNodes* src)
{
	dst->curr_size = src->curr_size;
	for (int i = 0; i < dst->nnodes; i++)
	{
		dst->indices[i] = src->indices[i];
	}
	for (int i = 0; i < src->curr_size; i++)
	{
		dst->nodes[i] = src->nodes[i];
	}
}

/* **************************************************************************************************
*							ORDERED INDEXED ARRAY
************************************************************************************************** */
char OIA_insert(IndexedValue* oiarr, int size, IndexedValue newelem)
{
	IndexedValue swap_elem;
	char inserted = 0;

	// if eligible place at last place
	if (oiarr[0].value > newelem.value)
	{
		oiarr[0] = newelem;
		inserted = 1;
	}
	// use insertion (from insertion sort)
	for (int i = 1; i < size && oiarr[i].value > oiarr[i - 1].value; i++)
	{
		// swap
		swap_elem = oiarr[i];
		oiarr[i] = oiarr[i - 1];
		oiarr[i - 1] = swap_elem;
	}

	return inserted;
}

void OIA_clear(IndexedValue* oiarr, int size)
{
	for (int i = 0; i < size; i++)
	{
		oiarr[i].index.arr[0] = -1;	oiarr[i].value = INFINITY;
	}
}

IndexedValue OIA_pack1(int index, double value)
{
	return (IndexedValue) { (Index) { {index, -1, -1} }, value };
}
IndexedValue OIA_pack3(int a, int b, int c, double value)
{
	return (IndexedValue) { (Index) { {a, b, c} }, value };
}

IndexedValue OIA_choose(IndexedValue* oiarr, int size, char includebest)
{
	// compute max size (max index + 1)
	int max_size = includebest ? size : size - 1;
	// compute min index (avoid unasigned elements)
	int min_idx;
	for (min_idx = 0; min_idx < size && oiarr[min_idx].index.arr[0] == -1; min_idx++);
	// choose random index
	int rand_idx = (int)(random() * (max_size - min_idx)) + min_idx;

	return oiarr[rand_idx];
}

void IV_deepcopy(IndexedValue* dst, IndexedValue* src)
{
	dst->index = src->index;
	dst->value = dst->value;
}

/* **************************************************************************************************
*							TABU LIST
************************************************************************************************** */

TabuList* TABU_new(int nnodes, int tenure)
{
	// allocate tabu list
	TabuList* tabu = NULL;	malloc_s(tabu, TabuList);

	// assign tabu nodes and tenure
	arr_malloc_s(tabu->tabu_nodes, nnodes, int);
	// fill tabu with minimum integer time
	for (int i = 0; i < nnodes; i++) tabu->tabu_nodes[i] = -1;
	tabu->tenure = tenure;
	tabu->now = 0;
	tabu->nnodes = nnodes;

	return tabu;
}

char TABU_istabu(TabuList* tabu, int node)
{
	if (tabu->tabu_nodes[node] < 0)	return 0;
	else return (tabu->now - tabu->tabu_nodes[node]) <= tabu->tenure;
		
}

void TABU_free(TabuList* tabu)
{
	free(tabu->tabu_nodes);
	free(tabu);
}

void TABU_print(TabuList* tabu)
{
	printf("**** TABU LIST ****\n");
	for (int i = 0; i < tabu->nnodes; i++)
	{
		if (TABU_istabu(tabu, i)) printf("%d ", i);
	}
	printf("\n");
}