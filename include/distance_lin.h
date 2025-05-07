#ifndef POLYGONCLUSTERING_DISTANCE_H
#define POLYGONCLUSTERING_DISTANCE_H

#include <gurobi_c++.h>
#include <iostream>
#include "../include/graph_boost.h"

void distance(Graph_boost &G, Graph_boost &H, double cGlobal, double cLocal, double cInsDel);
double distanceVal_lin(Graph_boost &G, Graph_boost &H, double cGlobal, double cLocal, double cInsDel);
double distanceVal_lin2(Graph_boost &G, Graph_boost &H, double cGlobal, double cLocal, double cInsDel);
double distanceVal_lin3(const Graph_boost &G, const Graph_boost &H, double cGlobal, double cLocal, double cInsDel);
double distVal_est(const Graph_boost &G, const Graph_boost &H, double E1, double E2, double longest1, double longest2);

class distance {

};


#endif //POLYGONCLUSTERING_DISTANCE_H
