#ifndef POLYGON_CLUSTERING_BINARY_IO_H
#define POLYGON_CLUSTERING_BINARY_IO_H


#include <fstream>
#include <iostream>
#include <vector>
#include "graph_boost.h"

void writeGraphVectorToBin(const std::vector<Graph_boost>& graphVector, const std::string& filepath);
void readGraphVectorFromBin(std::vector<Graph_boost>& graphVector, const std::string& filePath);

#endif //POLYGON_CLUSTERING_BINARY_IO_H
