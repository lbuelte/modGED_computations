#include <iostream>
#include <fstream>
#include "../include/distance_lin.h"
#include "../include/binary_io.h"
#include "../include/shapefile_io_operations.h"
#include "../include/hierachicalClustering3.h"

// TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or
// click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
Graph_boost testgraph1() {
    Graph_boost G;
    Point_boost p1 = {0, 0},
            p2 = {2, 0}, p3 = {1, 1}, p4 = {1, 2},
            p5 = {0, 3}, p6 = {2, 3}, p7 = {2, 2};
    Vertex_boost vertex1 = boost::add_vertex({0, 0, p1}, G);
    Vertex_boost vertex2 = boost::add_vertex({0, 0, p2}, G);
    Vertex_boost vertex3 = boost::add_vertex({0, 0, p3}, G);
    Vertex_boost vertex4 = boost::add_vertex({0, 0, p4}, G);
    Vertex_boost vertex5 = boost::add_vertex({0, 0, p5}, G);
    Vertex_boost vertex6 = boost::add_vertex({0, 0, p6}, G);
    Vertex_boost vertex7 = boost::add_vertex({0, 0, p7}, G);

    boost::add_edge(vertex1, vertex3, {static_cast<double>(boost::geometry::distance(p1, p3)), -1, 0.0}, G);
    boost::add_edge(vertex2, vertex3, {static_cast<double>(boost::geometry::distance(p2, p3)), -1, 0.0}, G);
    boost::add_edge(vertex3, vertex4, {static_cast<double>(boost::geometry::distance(p3, p4)), -1, 0.0}, G);
    boost::add_edge(vertex4, vertex5, {static_cast<double>(boost::geometry::distance(p4, p5)), -1, 0.0}, G);
    boost::add_edge(vertex4, vertex6, {static_cast<double>(boost::geometry::distance(p4, p6)), -1, 0.0}, G);
    boost::add_edge(vertex3, vertex7, {static_cast<double>(boost::geometry::distance(p3, p7)), -1, 0.0}, G);

    return G;
}

Graph_boost testgraph2() {
    Graph_boost G;
    Point_boost p1 = {3, 0},
            p2 = {3, 2}, p3 = {4, 1}, p4 = {5, 1},
            p5 = {6, 0}, p6 = {6, 2};
    Vertex_boost vertex1 = boost::add_vertex({0, 0, p1}, G);
    Vertex_boost vertex2 = boost::add_vertex({0, 0, p2}, G);
    Vertex_boost vertex3 = boost::add_vertex({0, 0, p3}, G);
    Vertex_boost vertex4 = boost::add_vertex({0, 0, p4}, G);
    Vertex_boost vertex5 = boost::add_vertex({0, 0, p5}, G);
    Vertex_boost vertex6 = boost::add_vertex({0, 0, p6}, G);

    boost::add_edge(vertex1, vertex3, {static_cast<double>(boost::geometry::distance(p1, p3)), -1, 0.0}, G);
    boost::add_edge(vertex2, vertex3, {static_cast<double>(boost::geometry::distance(p2, p3)), -1, 0.0}, G);
    boost::add_edge(vertex3, vertex4, {static_cast<double>(boost::geometry::distance(p3, p4)), -1, 0.0}, G);
    boost::add_edge(vertex4, vertex5, {static_cast<double>(boost::geometry::distance(p4, p5)), -1, 0.0}, G);
    boost::add_edge(vertex4, vertex6, {static_cast<double>(boost::geometry::distance(p4, p6)), -1, 0.0}, G);

    return G;
}


Graph_boost testgraph3() {
    Graph_boost G;
    Point_boost p1 = {0, 0}, p2 = {0.25, 0.25}, p3 = {0.5, 0.5},
            p4 = {2, 0}, p5 = {1.5, 0.5}, p6 = {1.25, 0.75},
            p7 = {1, 1},
            p8 = {1, 2},
            p11 = {0, 3}, p10 = {0.5, 2.5}, p9 = {0.75, 2.25},
            p14 = {2, 3}, p13 = {1.5, 2.5}, p12 = {1.25, 2.25};
    Vertex_boost vertex1 = boost::add_vertex({0, 0, p1}, G);
    Vertex_boost vertex2 = boost::add_vertex({0, 0, p2}, G);
    Vertex_boost vertex3 = boost::add_vertex({0, 0, p3}, G);
    Vertex_boost vertex4 = boost::add_vertex({0, 0, p4}, G);
    Vertex_boost vertex5 = boost::add_vertex({0, 0, p5}, G);
    Vertex_boost vertex6 = boost::add_vertex({0, 0, p6}, G);
    Vertex_boost vertex7 = boost::add_vertex({0, 0, p7}, G);
    Vertex_boost vertex8 = boost::add_vertex({0, 0, p8}, G);
    Vertex_boost vertex9 = boost::add_vertex({0, 0, p9}, G);
    Vertex_boost vertex10 = boost::add_vertex({0, 0, p10}, G);
    Vertex_boost vertex11 = boost::add_vertex({0, 0, p11}, G);
    Vertex_boost vertex12 = boost::add_vertex({0, 0, p12}, G);
    Vertex_boost vertex13 = boost::add_vertex({0, 0, p13}, G);
    Vertex_boost vertex14 = boost::add_vertex({0, 0, p14}, G);

    boost::add_edge(vertex1, vertex2, {static_cast<double>(boost::geometry::distance(p1, p3)), -1, 0.0}, G);
    boost::add_edge(vertex2, vertex3, {static_cast<double>(boost::geometry::distance(p1, p3)), -1, 0.0}, G);
    boost::add_edge(vertex3, vertex7, {static_cast<double>(boost::geometry::distance(p1, p3)), -1, 0.0}, G);
    boost::add_edge(vertex4, vertex5, {static_cast<double>(boost::geometry::distance(p1, p3)), -1, 0.0}, G);
    boost::add_edge(vertex5, vertex6, {static_cast<double>(boost::geometry::distance(p2, p3)), -1, 0.0}, G);
    boost::add_edge(vertex6, vertex7, {static_cast<double>(boost::geometry::distance(p3, p4)), -1, 0.0}, G);
    boost::add_edge(vertex7, vertex8, {static_cast<double>(boost::geometry::distance(p4, p5)), -1, 0.0}, G);
    boost::add_edge(vertex8, vertex9, {static_cast<double>(boost::geometry::distance(p4, p6)), -1, 0.0}, G);
    boost::add_edge(vertex9, vertex10, {static_cast<double>(boost::geometry::distance(p3, p7)), -1, 0.0}, G);
    boost::add_edge(vertex10, vertex11, {static_cast<double>(boost::geometry::distance(p3, p7)), -1, 0.0}, G);
    boost::add_edge(vertex8, vertex12, {static_cast<double>(boost::geometry::distance(p3, p7)), -1, 0.0}, G);
    boost::add_edge(vertex12, vertex13, {static_cast<double>(boost::geometry::distance(p3, p7)), -1, 0.0}, G);
    boost::add_edge(vertex13, vertex14, {static_cast<double>(boost::geometry::distance(p3, p7)), -1, 0.0}, G);

    return G;
}


int main() {
    std::vector<Graph_boost> graphVec, graphVec2, graphVec3;
    std::vector<Graph_boost> graphVec5, graphVec10, graphVec20, graphVec30;
    std::vector<int> resultVec;
    graphVec.push_back(testgraph1());
    graphVec.push_back(testgraph2());
    //std::cout<<distanceVal_lin3(graphVec[0], graphVec[1], 1, 1, 1)<<std::endl;
    //distance(graphVec[0], graphVec[1], 1, 1, 1);
    std::vector<Graph_boost> medialAxes;
    readGraphVectorFromBin(medialAxes, "../input/endenich.bin");
    //readGraphVectorFromBin(medialAxes, "../input/zentrum.bin");
    //writeToShapeFile(medialAxes, "../output/osmdings", 0.0);
    int i = 0;
    for (auto const &graph: medialAxes) {
        // if(num_vertices(graph)<100) std::cout<<i<< " :   "<< num_vertices(graph)<<std::endl;
        graphVec2.push_back(graph);
        i++;
        if (i == 50) break;
    }
    //graphVec.push_back(simplifyGraphPaths(medialAxes[817], 1));
    //graphVec.push_back(simplifyGraphPaths(medialAxes[818], 1));
    double edgeSum, longestCurr;
    std::vector<double> edgeSums;
    std::vector<double> longest;
    std::cout << "Checkpoint1\n";
    /*for (auto graph: medialAxes) {
        //auto nextrgraph = simplifyGraphPathsWithEdgeSum(graph, 1, edgeSum, longestCurr);
        auto nextrgraph5 = simplifyAndFilterGraphPathsWithEdgeSum(graph, 1, edgeSum, longestCurr, 0.00);
        auto nextrgraph10 = simplifyAndFilterGraphPathsWithEdgeSum(graph, 1, edgeSum, longestCurr, 0.18);
        auto nextrgraph20 = simplifyAndFilterGraphPathsWithEdgeSum(graph, 1, edgeSum, longestCurr, 0.35);

        auto nextrgraph30 = simplifyAndFilterGraphPathsWithEdgeSum(graph, 1, edgeSum, longestCurr, 0.7);

        //std::cout<< "Checkpoint2\n";
        //std::cout<< num_edges(graph)-num_edges(resultVec[i])<<std::endl;
        graphVec5.push_back(nextrgraph5);
        graphVec10.push_back(nextrgraph10);
        graphVec20.push_back(nextrgraph20);
        graphVec30.push_back(nextrgraph30);

        edgeSums.push_back(edgeSum);
        longest.push_back(longestCurr);

    }*/



    //writeToShapeFile(medialAxes, "../output/osmdings_org", 0.0);
    //graphVec3.push_back(graphVec2[5]);
    //graphVec3.push_back(graphVec2[19]);
    //std::cout<< num_edges(graphVec2[342])<<std::endl;
    //std::cout<< num_edges(graphVec2[1923])<<std::endl;
    /*resultVec.push_back(simplifyGraphPaths(medialAxes[0], 1));
    resultVec.push_back(simplifyGraphPaths(medialAxes[1], 1));

    std::cout<<num_vertices(resultVec[0])<<"   "<<num_vertices(resultVec[1])<<std::endl;
    std::cout<<num_edges(resultVec[0])<<"   "<<num_edges(resultVec[1])<<std::endl;
    std::cout<<num_vertices(medialAxes[0])<<"   "<<#include <fstream>num_vertices(medialAxes[1])<<std::endl;
    std::cout<<num_edges(medialAxes[0])<<"   "<<num_edges(medialAxes[1])<<std::endl;*/

    //std::cout<<distanceVal_lin3(graphVec[0], graphVec[1], 1, 1, 1)<<std::endl;
    //std::cout<<distanceVal_lin3(graphVec3[0], graphVec3[1], 1, 1, 1)<<std::endl;

    //std::cout<< num_edges(resultVec[0])<<std::endl;
    //resultVec=kMedianClustering(graphVec, 4);
    std::cout << "Berechnung geglückt\n";
    //MultiPolygon_boost source;
    //ReadShapeFile("../input/shp/osm_siegburg", source);
    //writeToShapeFile(source, graphVec, "../output/osmdings_smpl");
    //writeToShapeFile(medialAxes, "../output/osmdings_smpl2", 0.0);
    //writeToShapeFile(graphVec2, "../output/osmdings_smpl", 0.0);
    //writeToShapeFile(graphVec5, "../output/endenich_5", 0.0);
    //writeToShapeFile(graphVec10, "../output/endenich_10", 0.0);
    //writeToShapeFile(graphVec20, "../output/endenich_20", 0.0);
    //writeToShapeFile(graphVec30, "../output/endenich_30", 0.0);
    //writeToShapeFile(medialAxes, "../output/endenich_30_2", 0.54);

    //auto stuff = kMedoidsWithSampling(graphVec5, 3, 100, edgeSums, longest, 5, "../output/endenich00.txt");
    //auto stuff1 = kMedoidsWithSampling(graphVec10, 3, 100, edgeSums, longest, 5, "../output/endenich10.txt");
    //auto stuff2 = kMedoidsWithSampling(graphVec20, 3, 100, edgeSums, longest, 5, "../output/endenich20.txt");
    //auto stuff3 = kMedoidsWithSampling(graphVec30, 3, 100, edgeSums, longest, 5, "../output/endenich30.txt");
    //auto stuff = clarans(graphVec2, 3, edgeSums, longest, 5, 250);
    std::ofstream datei("../output/timelogTotalh1k3zentrum.txt", std::ios::app);
    if (!datei) {
        std::cerr << "Fehler beim Öffnen der Datei!\n";
        return 1;
    }
    auto start = std::chrono::high_resolution_clock::now();
    //auto stuff = get_hierarchicalClustering(medialAxes, 1, 3, 1, 0.87, "../output/endenich50h1.txt");
    auto stuff = get_hierarchicalClustering(medialAxes, 2, 2, 1, 0.7, "../output/zentrum40h2.txt");
    auto ende = std::chrono::high_resolution_clock::now();
    auto dauer = std::chrono::duration_cast<std::chrono::milliseconds>(ende - start);

    datei << "10    " << dauer.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    auto stuff1 = get_hierarchicalClustering(medialAxes, 1, 4, 1, 0.7, "../output/zentrum40k4.txt");
    ende = std::chrono::high_resolution_clock::now();
    dauer = std::chrono::duration_cast<std::chrono::milliseconds>(ende - start);
    datei << "0    " << dauer.count() << std::endl;
    /*start = std::chrono::high_resolution_clock::now();
    auto stuff2 = get_hierarchicalClustering(medialAxes, 1, 3, 1, 0.52, "../output/zentrum30h1.txt");
    ende = std::chrono::high_resolution_clock::now();
    dauer = std::chrono::duration_cast<std::chrono::milliseconds>(ende - start);
    datei << "30    " << dauer.count() << std::endl;
    start = std::chrono::high_resolution_clock::now();
    auto stuff3 = get_hierarchicalClustering(medialAxes, 1, 3, 1, 0.35, "../output/zentrum20h1.txt");
    ende = std::chrono::high_resolution_clock::now();
    dauer = std::chrono::duration_cast<std::chrono::milliseconds>(ende - start);
    datei << "20    " << dauer.count() << std::endl;*/
    datei.close();
    //auto stuff = get_hierarchicalClustering(graphVec2, 2, 3, 1);

    //std::vector<int> clsut = run_clustering_pipeline(graphVec2, 3, 5, 0.5);

    //writeToShapeFile(graphVec, "../output/example", 0.0);
    //distance(medialAxes[342], medialAxes[1923], 1, 1, 1);
    /*std::vector<Graph_boost> exmplGraph = {testgraph3()};
    writeToShapeFile(exmplGraph, "../output/osmdings_org", 0.0);
    simplifyGraphPaths(exmplGraph[0], 0.1);
    writeToShapeFile(exmplGraph, "../output/osmdings_smpl", 0.0);*/




    return EXIT_SUCCESS;
}
