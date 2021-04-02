#ifndef CPX_MODELS_H_  

#define CPX_MODELS_H_

#include "tsp.h"

void build_model_base_undirected(graph* g, modeltype mt, CPXENVptr env, CPXLPptr lp);
void build_model_base_directed(graph* g, modeltype mt, CPXENVptr env, CPXLPptr lp);
void build_model_mtz(graph* g, modeltype mt, CPXENVptr env, CPXLPptr lp);
void build_model_gg(graph* g, modeltype mt, CPXENVptr env, CPXLPptr lp);

#endif