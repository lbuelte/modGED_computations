#ifndef POLYGON_CLUSTERING_GRAPH_BOOST_H
#define POLYGON_CLUSTERING_GRAPH_BOOST_H

#include <cmath>
#include <cstdio>
#include <ctime>
#include <string>

// This will work properly only with GCC compiler.
#include <ieee754.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/polygon/point_concept.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/io/wkt/wkt.hpp>
#include <boost/polygon/voronoi.hpp>
#include <boost/geometry/geometries/linestring.hpp>

#include <boost/graph/adjacency_list.hpp>
//#include "shapefil.h"


#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/algorithms/convex_hull.hpp>
#include <boost/geometry/algorithms/transform.hpp>
#include <boost/geometry/strategies/transform.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/io/io.hpp>
#include <fstream>
#include <boost/geometry/geometries/segment.hpp>

#include <boost/functional/hash.hpp>

//#include "voronoi_boost_advanced.h"
#include <utility>


using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;


// Define the Boost Polygon type
//typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> Point;










namespace bg = boost::geometry;
namespace bgm = bg::model;
namespace bgs = bg::strategy;





using namespace boost;


struct Point_boost {
    long double a;
    long double b;

    Point_boost(long double x = 0, long double y = 0) : a(x), b(y) {}

    // Custom comparison for std::map
    bool operator<(const Point_boost &other) const {
        return std::tie(a, b) < std::tie(other.a, other.b);
    }

    bool operator==(const Point_boost &other) const {
        return std::tie(a, b) == std::tie(other.a, other.b);
    }

    bool operator!=(const Point_boost &other) const {
        return !(std::tie(a, b) == std::tie(other.a, other.b));
    }

    Point_boost operator+(const Point_boost &other) const {
        return Point_boost(a + other.a, b + other.b);
    }

    Point_boost operator-(const Point_boost &other) const {
        return Point_boost(a - other.a, b - other.b);
    }

    Point_boost operator/(const double &divisor) const {
        return Point_boost(a / divisor, b / divisor);
    }
};

//BOOST_GEOMETRY_REGISTER_POINT_2D(Point_boost, long double, boost::geometry::cs::cartesian, a, b)

namespace boost {
    namespace polygon {
        template<>
        struct geometry_concept<Point_boost> {
            typedef point_concept type;
        };

        template<>
        struct point_traits<Point_boost> {
            typedef long double coordinate_type;

            static inline coordinate_type get(const Point_boost &point, orientation_2d orient) {
                return (orient == HORIZONTAL) ? point.a : point.b;
            }
        };
    }
}


/*bool operator<(const Point_boost& u, const Point_boost& v)
{
    // Usage of std::tie :
    // compares a.x to b.x,
    // then a.y to b.y
    return std::tie(u.a, u.b) < std::tie(v.a, v.b);
}*/
struct VertexProps_boost {
    //radius specifying the medial axis in this point
    double minRadius;
    double maxRadius;
    //location of the corresponding point in the plane
    Point_boost p;
};

struct EdgeProps_boost {
    //edge weights
    long double weight;
    //integer label of the edge
    int label;
};

struct VAEdgeProps_boost {
    //edge weight
    double weight;
    //integer label of the edge
    int label;
    //Vanishing Angle (according to )
    double vanishingAngle;
};


//typedef adjacency_list<vecS, vecS, undirectedS, VertexProps, property<edge_weight_t, double>> Graph;
typedef boost::polygon::voronoi_diagram<double> VoronoiDiagram_boost;

BOOST_GEOMETRY_REGISTER_POINT_2D(Point_boost, long double, bg::cs::cartesian, a, b)
static constexpr bool closed_polygons = false;
using Polygon_boost = bgm::polygon<Point_boost, false, closed_polygons>;


//typedef boost::geometry::model::polygon<Point> Polygon;
typedef boost::geometry::model::multi_polygon<Polygon_boost> MultiPolygon_boost;
typedef adjacency_list<vecS, vecS, undirectedS, VertexProps_boost, VAEdgeProps_boost> Graph_boost;
typedef adjacency_list<vecS, vecS, bidirectionalS, VertexProps_boost, VAEdgeProps_boost> DirectedGraph_boost;
typedef std::pair<int, int> Edge_boost;
typedef typename graph_traits<Graph_boost>::vertex_descriptor Vertex_boost;
typedef boost::geometry::model::segment<Point_boost> Segment_boost;
typedef std::map<Point_boost, unsigned long> VH_Vertex_Map_boost;
typedef std::pair<Vertex_boost, Vertex_boost> SimpleEdge;
typedef graph_traits<DirectedGraph_boost>::vertex_iterator VertexIterator;
typedef graph_traits<DirectedGraph_boost>::edge_iterator EdgeIterator;
typedef graph_traits<DirectedGraph_boost>::out_edge_iterator OutEdgeIterator;
typedef graph_traits<DirectedGraph_boost>::in_edge_iterator InEdgeIterator;
typedef boost::geometry::model::linestring<Point_boost> Linestring;


struct EdgeComparator {
    bool
    operator()(const DirectedGraph_boost::edge_descriptor &e1, const DirectedGraph_boost::edge_descriptor &e2) const {
        return std::tie(e1.m_source, e1.m_target) > std::tie(e2.m_source, e2.m_target);
    }
};

struct SimpleEdgeHash {
    std::size_t operator()(const std::pair<unsigned long, unsigned long> &edge) const {
        std::size_t h1 = std::hash<unsigned long>()(edge.first);
        std::size_t h2 = std::hash<unsigned long>()(edge.second);
        // Asymmetrische Kombination von h1 und h2
        return h1 ^ (h2 * 0x9e3779b9 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    }
};

// Define equality operator for std::pair
struct SimpleEdgeEqual {
    bool operator()(const std::pair<unsigned long, unsigned long> &e1,
                    const std::pair<unsigned long, unsigned long> &e2) const {
        return e1.first == e2.first && e1.second == e2.second;
    }
};

Graph_boost simplifyGraphPaths(Graph_boost &graph, double epsilon);

Graph_boost simplifyGraphPathsWithEdgeSum(Graph_boost &graph, double epsilon, double &edgeSum, double &longest);
Graph_boost simplifyAndFilterGraphPathsWithEdgeSum(Graph_boost graph, double epsilon, double &edgeSum, double &longest, double minAngle);


#endif // !POLYGON_CLUSTERING_GRAPH_BOOST_H
