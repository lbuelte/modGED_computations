#ifndef _shapefile_io_operations_included_
#define _shapefile_io_operations_included_

#include "shapefil.h"
//#include "medial_axis_boost.h"

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include "graph_boost.h"

void ReadShapeFile(const std::string &filename, MultiPolygon_boost& multiPolygon);

void writeToShapeFile(std::vector<Graph_boost> mas, std::string path, double angle);

void writeToShapeFile(MultiPolygon_boost polys, std::vector<int> clusters, std::string path);


void writeToShapeFile(MultiPolygon_boost polys, std::string path);


#endif