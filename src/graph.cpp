#include "../include/graph_boost.h"


unsigned long
findOrAddVertex(VH_Vertex_Map_boost &c, Point_boost item, Graph_boost &ma, double minRadius, double maxRadius) {

    if (c.count(item) == 0) {

        auto vertex = boost::add_vertex({minRadius, maxRadius, item}, ma);


        c.emplace(item, vertex);
        return vertex;
    }

    return c.at(item);
}

struct CuteDFS : public boost::default_dfs_visitor {
    std::vector<Linestring> &paths;
    const Graph_boost &graph;
    Linestring currentPath;

    CuteDFS(std::vector<Linestring> &paths, const Graph_boost &g) : paths(paths), graph(g) {}

    void tree_edge(boost::graph_traits<Graph_boost>::edge_descriptor e, const Graph_boost &g) {
        if (g[e].vanishingAngle != 0) {
            int v = target(e, g); // Zielknoten der Kante
            int u = source(e, g); // Zielknoten der Kante
            if (currentPath.empty()) currentPath.push_back(g[u].p);
            currentPath.push_back(g[v].p);
            //std::cout << "Überquere Kante von: " << u << " (" << g[u].p.a << ", " << g[u].p.b << "zu: " << v << " (" << g[v].p.a << ", " << g[v].p.b << ")\n";
            if (degree(v, g) != 4) {
                paths.push_back(currentPath);
                currentPath.clear();
                //std::cout<<degree(v, g)<<std::endl;
            }
        } else std::cout << "well";
        //linestring.push_back({g[v].p, g[v].y});

    }
};

int find_start_vertex(const Graph_boost &g) {
    for (auto v: boost::make_iterator_range(vertices(g))) {
        if (boost::degree(v, g) != 2 && boost::degree(v, g) != 0) {
            return v;  // Erster Knoten mit Grad ≠ 2 gefunden
        }
    }
    return *vertices(g).first;  // Falls alle Knotengrade 2 sind, wähle den ersten Knoten
}

Graph_boost simplifyGraphPaths(Graph_boost &graph, double epsilon) {
    std::vector<Linestring> paths;
    std::unordered_set<Vertex_boost> visited;
    Linestring ls;
    VH_Vertex_Map_boost vertex_handles;

    Graph_boost simplifiedGraph;


    // DFS auf allen Knoten mit Grad != 2 starten
    /*for (auto vp = vertices(graph); vp.first != vp.second; ++vp.first) {
        Vertex_boost v = *vp.first;
        if (boost::degree(v, graph) != 2 && !visited.count(v)) {
            boost::depth_first_search(graph, boost::visitor(PathFinderDFS(paths, visited, graph)));
        }
    }*/

    //boost::depth_first_search(graph, boost::visitor(PathFinderDFS(paths, visited, graph)));
    boost::depth_first_search(graph, boost::visitor(CuteDFS(paths, graph)).root_vertex(find_start_vertex(graph)));

    /*std::cout << "Pfade:  " << paths.size() << std::endl;
    std::cout << "Pfade:  " << paths[0].size() << std::endl;
    std::cout << "Pfade:  " << paths[1].size() << std::endl;
    std::cout << "Pfade:  " << paths[2].size() << std::endl;*/

    for (const auto &path: paths) {
        //std::vector<Point_boost> points;
        Linestring simplified_lstr, points_lstr;
        /*for (Vertex_boost v : path) {
            //points.push_back(graph[v].p);
            points_lstr.push_back(graph[v].p);
        }*/

        std::vector<Point_boost> simplified;
        //bg::simplify(points, simplified, epsilon);


        bg::simplify(path, simplified_lstr, epsilon);




        /*std::cout<<"Länge des Pfades: "<<path.size()<<std::endl;
        for(auto vertex: path){
            std::cout<< "("<<vertex.a<<" , "<<vertex.b<<") ";
        }
        std::cout<<std::endl;
        std::cout<<"Länge des vereinfachten Pfades: "<<simplified_lstr.size()<<std::endl;
        for(auto vertex: simplified_lstr){
            std::cout<< "("<<vertex.a<<" , "<<vertex.b<<") ";
        }
        std::cout<<std::endl;
        // Entferne alte Knoten (außer Start/Ende)
        /*for (size_t i = 1; i < path.size() - 1; ++i) {
            clear_vertex(path[i], graph);
            remove_vertex(path[i], graph);
        }

        // Neue Kante zwischen Start- und Endknoten
        if (simplified.size() >= 2) {
            add_edge(path.front(), path.back(), graph);
        }*/
        Vertex_boost previousV, currentV;
        Point_boost prevPoint;
        for (auto point: simplified_lstr) {
            currentV = findOrAddVertex(vertex_handles, point, simplifiedGraph, 0, 0);
            if (point != path.front()) {
                double weight = boost::geometry::distance(prevPoint, point);
                add_edge(currentV, previousV,
                         {weight, -1, 0},
                         simplifiedGraph);
            }
            previousV = currentV;
            prevPoint = point;
        }

    }
    return simplifiedGraph;
}

Graph_boost simplifyGraphPathsWithEdgeSum(Graph_boost &graph, double epsilon, double &edgeSum, double &longest) {
    std::vector<Linestring> paths;
    std::unordered_set<Vertex_boost> visited;
    Linestring ls;
    VH_Vertex_Map_boost vertex_handles;

    Graph_boost simplifiedGraph;

    longest = 0;

    // DFS auf allen Knoten mit Grad != 2 starten
    /*for (auto vp = vertices(graph); vp.first != vp.second; ++vp.first) {
        Vertex_boost v = *vp.first;
        if (boost::degree(v, graph) != 2 && !visited.count(v)) {
            boost::depth_first_search(graph, boost::visitor(PathFinderDFS(paths, visited, graph)));
        }
    }*/

    //boost::depth_first_search(graph, boost::visitor(PathFinderDFS(paths, visited, graph)));
    boost::depth_first_search(graph, boost::visitor(CuteDFS(paths, graph)).root_vertex(find_start_vertex(graph)));

    /*std::cout << "Pfade:  " << paths.size() << std::endl;
    std::cout << "Pfade:  " << paths[0].size() << std::endl;
    std::cout << "Pfade:  " << paths[1].size() << std::endl;
    std::cout << "Pfade:  " << paths[2].size() << std::endl;*/

    for (const auto &path: paths) {
        //std::vector<Point_boost> points;
        Linestring simplified_lstr, points_lstr;
        /*for (Vertex_boost v : path) {
            //points.push_back(graph[v].p);
            points_lstr.push_back(graph[v].p);
        }*/

        std::vector<Point_boost> simplified;
        //bg::simplify(points, simplified, epsilon);


        bg::simplify(path, simplified_lstr, epsilon);




        /*std::cout<<"Länge des Pfades: "<<path.size()<<std::endl;
        for(auto vertex: path){
            std::cout<< "("<<vertex.a<<" , "<<vertex.b<<") ";
        }
        std::cout<<std::endl;
        std::cout<<"Länge des vereinfachten Pfades: "<<simplified_lstr.size()<<std::endl;
        for(auto vertex: simplified_lstr){
            std::cout<< "("<<vertex.a<<" , "<<vertex.b<<") ";
        }
        std::cout<<std::endl;
        // Entferne alte Knoten (außer Start/Ende)
        /*for (size_t i = 1; i < path.size() - 1; ++i) {
            clear_vertex(path[i], graph);
            remove_vertex(path[i], graph);
        }

        // Neue Kante zwischen Start- und Endknoten
        if (simplified.size() >= 2) {
            add_edge(path.front(), path.back(), graph);
        }*/
        Vertex_boost previousV, currentV;
        Point_boost prevPoint;
        edgeSum = 0;
        for (auto point: simplified_lstr) {
            currentV = findOrAddVertex(vertex_handles, point, simplifiedGraph, 0, 0);
            if (point != path.front()) {
                double weight = boost::geometry::distance(prevPoint, point);
                add_edge(currentV, previousV,
                         {weight, -1, 0},
                         simplifiedGraph);
                edgeSum += weight;
                if (weight > longest) longest = weight;
            }
            previousV = currentV;
            prevPoint = point;
        }

    }
    return simplifiedGraph;
}




Graph_boost simplifyAndFilterGraphPathsWithEdgeSum(Graph_boost graph, double epsilon, double &edgeSum, double &longest,
                                                   double minAngle) {
    std::vector<Linestring> paths;
    std::unordered_set<Vertex_boost> visited;
    Linestring ls;
    VH_Vertex_Map_boost vertex_handles;
    std::vector<std::pair<unsigned long, unsigned long>> toRemove;
    Graph_boost simplifiedGraph;
    edgeSum = 0;

    longest = 0;
    for (auto e : make_iterator_range(edges(graph))) {
        auto u = source(e, graph);
        auto v = target(e, graph);
        auto [e1, exists1] = edge(u, v, graph);
        auto [e2, exists2] = edge(v, u, graph);
        double va1 = exists1 ? graph[e1].vanishingAngle : 0.0;
        double va2 = exists2 ? graph[e2].vanishingAngle : 0.0;
        if(va1 < minAngle && va2<minAngle) toRemove.push_back({source(e, graph), target(e, graph)});
    }
    for (auto e: toRemove) {
        if(edge(e.first, e.second, graph).second) remove_edge(e.first, e.second, graph);
    }




    // DFS auf allen Knoten mit Grad != 2 starten
    /*for (auto vp = vertices(graph); vp.first != vp.second; ++vp.first) {
        Vertex_boost v = *vp.first;
        if (boost::degree(v, graph) != 2 && !visited.count(v)) {
            boost::depth_first_search(graph, boost::visitor(PathFinderDFS(paths, visited, graph)));
        }
    }*/

    //boost::depth_first_search(graph, boost::visitor(PathFinderDFS(paths, visited, graph)));
    boost::depth_first_search(graph, boost::visitor(CuteDFS(paths, graph)).root_vertex(
            find_start_vertex(graph)));


    /*std::cout << "Pfade:  " << paths.size() << std::endl;
    std::cout << "Pfade:  " << paths[0].size() << std::endl;
    std::cout << "Pfade:  " << paths[1].size() << std::endl;
    std::cout << "Pfade:  " << paths[2].size() << std::endl;*/

    for (const auto &path: paths) {
        //std::vector<Point_boost> points;
        Linestring simplified_lstr, points_lstr;
        /*for (Vertex_boost v : path) {
            //points.push_back(graph[v].p);
            points_lstr.push_back(graph[v].p);
        }*/

        std::vector<Point_boost> simplified;
        //bg::simplify(points, simplified, epsilon);


        bg::simplify(path, simplified_lstr, epsilon);




        /*std::cout<<"Länge des Pfades: "<<path.size()<<std::endl;
        for(auto vertex: path){
            std::cout<< "("<<vertex.a<<" , "<<vertex.b<<") ";
        }
        std::cout<<std::endl;
        std::cout<<"Länge des vereinfachten Pfades: "<<simplified_lstr.size()<<std::endl;
        for(auto vertex: simplified_lstr){
            std::cout<< "("<<vertex.a<<" , "<<vertex.b<<") ";
        }
        std::cout<<std::endl;
        // Entferne alte Knoten (außer Start/Ende)
        /*for (size_t i = 1; i < path.size() - 1; ++i) {
            clear_vertex(path[i], graph);
            remove_vertex(path[i], graph);
        }

        // Neue Kante zwischen Start- und Endknoten
        if (simplified.size() >= 2) {
            add_edge(path.front(), path.back(), graph);
        }*/
        Vertex_boost previousV, currentV;
        Point_boost prevPoint;
        for (auto point: simplified_lstr) {
            currentV = findOrAddVertex(vertex_handles, point, simplifiedGraph, 0, 0);
            if (point != path.front()) {
                double weight = boost::geometry::distance(prevPoint, point);
                add_edge(currentV, previousV,
                         {weight, -1, 0},
                         simplifiedGraph);
                edgeSum += weight;
                if (weight > longest) longest = weight;
            }
            previousV = currentV;
            prevPoint = point;
        }

    }
    return simplifiedGraph;
}
//example code
/*void simplifyGraph(Graph_boost& graph, double tolerance) {
    Graph_boost simplifiedGraph;

    // Iteriere durch alle Kanten
    for (auto [ei, ei_end] = boost::edges(graph); ei != ei_end; ++ei) {
        auto u = boost::source(*ei, graph);
        auto v = boost::target(*ei, graph);

        // Extrahiere die Kanten als Linienzug
        Linestring line{graph[u].p, graph[v].p};
        Linestring simplifiedLine;

        // Wende Douglas-Peucker-Simplifizierung an
        bg::simplify(line, simplifiedLine, tolerance);

        // Neuen Graphen mit vereinfachten Knoten erstellen
        auto prev = boost::add_vertex(simplifiedGraph);
        simplifiedGraph[prev] = simplifiedLine.front();

        for (size_t i = 1; i < simplifiedLine.size(); ++i) {
            auto next = boost::add_vertex(simplifiedGraph);
            simplifiedGraph[next] = simplifiedLine[i];
            boost::add_edge(prev, next, simplifiedGraph);
            prev = next;
        }
    }

    // Ersetze den ursprünglichen Graphen durch den neuen
    graph = std::move(simplifiedGraph);
}*/


void createGraph() {
    //initialize graph
    Graph_boost g;

    //add vertices to graph
    int v0 = add_vertex(g);
    int v1 = add_vertex(g);

    //add edge without weight information or label
    add_edge(v0, v1, g);

    //add edge with weight information and label
    add_edge(v0, v1, {1.0, 1}, g);
}