// #include "graph_boost.h"
#include "graph.hpp"
#include "modGED.hpp"

graph<Point, EdgeLength> testgraph1(){
    graph<Point, EdgeLength> g;
    Point p1 = {1, 1}, 
          p2 = {0, 0}, 
          p3 = {0, 2}, 
          p4 = {5, 1},
          p5 = {6, 2}, 
          p6 = {6, 0};
    g.add_node(0, p1);
    g.add_node(1, p2);
    g.add_node(2, p3);
    g.add_node(3, p4);
    g.add_node(4, p5);
    g.add_node(5, p6);

    auto add_edge = [] (graph<Point, EdgeLength> &g, node i, node j) {
        g.add_edge(i, j, compute_euclidean_distance(g.get_node_label(i), g.get_node_label(j)));
    };
    add_edge(g, 0, 1);
    add_edge(g, 0, 2);
    add_edge(g, 0, 3);
    add_edge(g, 3, 4);
    add_edge(g, 3, 5);
    
    return g;
}

graph<Point, EdgeLength> testgraph2(){
    graph<Point, EdgeLength> g;
    Point p1 = {1, 1}, 
          p2 = {0, 0}, 
          p3 = {0, 2}, 
          p4 = {4, 1},
          p5 = {5, 2}, 
          p6 = {4, -1},
          p7 = {3, -2},
          p8 = {5, -2};

    g.add_node(0, p1);
    g.add_node(1, p2);
    g.add_node(2, p3);
    g.add_node(3, p4);
    g.add_node(4, p5);
    g.add_node(5, p6);
    g.add_node(6, p7);
    g.add_node(7, p8);

    auto add_edge = [] (graph<Point, EdgeLength> &g, node i, node j) {
        g.add_edge(i, j, compute_euclidean_distance(g.get_node_label(i), g.get_node_label(j)));
    };
    add_edge(g, 0, 1);
    add_edge(g, 0, 2);
    add_edge(g, 0, 3);
    add_edge(g, 3, 4);
    add_edge(g, 3, 5);
    add_edge(g, 5, 6);
    add_edge(g, 5, 7);
    
    return g;
}

int main()
{
    auto g = testgraph1();
    auto h = testgraph2();

    auto print_graph = [] (graph<Point, EdgeLength> g){
        for (int i = 0; i < g.number_of_nodes(); i++){
            std::cout << "Node " << i << " : (" << g.get_node_label(i).first << "," << g.get_node_label(i).second << ")" << std::endl; 
        }
        for (int ij = 0; ij < g.number_of_edges(); ij++){
            std::cout << "Edge " << ij << " = (" << g.get_edge(ij).first << "," << g.get_edge(ij).second << ") : " << g.get_edge_label(g.get_edge(ij)) << std::endl;
        }
    };
    print_graph(g);
    std::cout << std::endl;
    print_graph(h);


    compute_modfied_GED(g, h, 1.0, 1.0);
}