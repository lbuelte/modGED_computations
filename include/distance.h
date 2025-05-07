#ifndef POLYGONCLUSTERING_DISTANCE_H
#define POLYGONCLUSTERING_DISTANCE_H

#include <gurobi_c++.h>
#include <iostream>
#include "../include/graph_boost.h"

void distance(Graph_boost &G, Graph_boost &H, double cGlobal, double cLocal, double cInsDel);
double distanceVal(const Graph_boost &G, const Graph_boost &H, double cGlobal, double cLocal, double cInsDel);

class distance {

};


#endif //POLYGONCLUSTERING_DISTANCE_H
