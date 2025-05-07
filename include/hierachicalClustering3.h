#ifndef POLYGONCLUSTERING_HIERACHICALCLUSTERING2_H
#define POLYGONCLUSTERING_HIERACHICALCLUSTERING2_H

#include "distance_lin.h"
#include <random>

std::vector<int>
get_hierarchicalClustering(const std::vector<Graph_boost> &input_graphs, int h, int k, double epsilon, double initCutOff, std::string filename);


class hierachicalClustering3 {

};


#endif //POLYGONCLUSTERING_HIERACHICALCLUSTERING2_H
