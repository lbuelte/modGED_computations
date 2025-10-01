#ifndef MODGED_H
#define MODGED_H

#include "graph_boost.h"

typedef tuple<Vertex_boost, Vertex_boost, Vertex_boost, double> DirectedPath2;

// The two differnt ILP formulations for modGED return the same solution up to numerical inaccuracies occuring in the ILP solver.
// However, the formulations differ in their performance on different instances.
// The function compute_modGED_F2Plus implements the ILP formulation F2+
// and the function compute_modGED_FORI implements the ILP formulation FORI, both as described in the paper 
// "Enhancing Graph Edit Distance Computation: Stronger and Orientation-based ILP Formulations"
// by Andrea D'Ascenzo, Julian Meffert, Petra Mutzel and Fabrizio Rossi (2025).
// When running both formulations on a large set of medial axis graphs, 
// we observed that FORI clearly dominates F2+ in total runtime, usually by a factor 4-5.
// Further, with a time limit of 60 seconds, FORI had significantly less timeouts than F2+, 
// though both formulations only struggled when both input graphs were structural outliers or had many vertices.
// We therefore recommend using FORI as the default formulation for modGED computation of medial axis graphs.

long double compute_modGED_F2Plus(const Graph_boost &G, const Graph_boost &H, double cGlobal, double cLocal, double cInsDel, bool &timeout);
long double compute_modGED_FORI(const Graph_boost &G, const Graph_boost &H, double cGlobal, double cLocal, double cInsDel, bool &timeout);


#endif // MODGED_H