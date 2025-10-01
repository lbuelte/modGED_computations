#include "modGED.h"

#include <numeric>
#include <gurobi_c++.h>

namespace{

// Returns the signed angle in radians between the line segments ab and bc, in (-π, π].
// Positive angles are counter-clockwise rotations from ab to bc. Negative angles are clockwise rotations.
long double compute_signed_angle_between_lines(Point_boost const &a, Point_boost const &b, Point_boost const &c){

    // v1 = a - b (incoming edge)
    long double vx1 = a.a - b.a;
    long double vy1 = a.b - b.b;

    // v2 = c - b (outgoing edge)
    long double vx2 = c.a - b.a;
    long double vy2 = c.b - b.b;

    // angle of v1 relative to x-axis
    long double theta = std::atan2(vy1, vx1);

    // rotate v2 by -theta, so v1 aligns with +x axis
    long double cos_t = std::cos(-theta);
    long double sin_t = std::sin(-theta);
    long double vx2r = vx2 * cos_t - vy2 * sin_t;
    long double vy2r = vx2 * sin_t + vy2 * cos_t;

    // now compute angle of rotated v2
    long double angle = std::atan2(vy2r, vx2r); // in (-π, π]

    return angle;
}


// Finds all directed paths of length 2 in the graph g, i.e., all triples of vertices (i, j, h)
// such that (i, j) and (j, h) are edges in g. The paths ijh and hji are considered different.
// Further, the function computes the signed angle between the two edges (i, j) and (j, h) 
// for each path and returns it together with the paths.
std::vector<DirectedPath2> findAllPathsOfLength2(Graph_boost const &g) {
    std::vector<DirectedPath2> paths;
    for (auto j: make_iterator_range(vertices(g))) {
        for (auto i: make_iterator_range(adjacent_vertices(j, g))) {
            if (i == j) continue;
            for (auto h: make_iterator_range(adjacent_vertices(j, g))) {
                if (h == j || h == i) continue;
                auto const angle = compute_signed_angle_between_lines(g[i].p, g[j].p, g[h].p);
                paths.emplace_back(i, j, h, angle);
            }
        }
    }
    return paths;
}


// Returns the sum of the variables in column 'index' of the 2D vector 'vars'.
GRBLinExpr colSum(const std::vector<std::vector<GRBVar>> &vars, int index) {
    GRBLinExpr sum = 0;
    for (const auto &vec: vars) sum += vec[index];
    return sum;
}


// Returns the sum of the variables in the vector 'vec'.
GRBLinExpr vecSum(const std::vector<GRBVar> &vec) {
    GRBLinExpr sum = 0;
    for (const auto &val: vec) sum += val;
    return sum;
}


// Computes the constant cost of deleting all edges of G and inserting all edges of H.
long double compute_gamma(const Graph_boost &G, const Graph_boost &H, double cInsDel) {
    double gamma = 0;
    for (auto e: boost::make_iterator_range(edges(G))) {
        gamma += G[e].weight * cInsDel;
    }
    for (auto e: boost::make_iterator_range(edges(H))) {
        gamma += H[e].weight * cInsDel;
    }
    return gamma;
}


// Computes the linear expression representing the objective cost based on the edge mapping variables y
// and the path mapping variables p. Here, the edges are undirected 
// and the variables y[ij][kl] represent the mapping of edge ij in G to edge kl in H.
GRBLinExpr compute_objective_cost_y(const Graph_boost &G, const Graph_boost &H, 
        double cLocal, double cGlobal, double cInsDel, const std::vector<std::vector<GRBVar>> &y, 
        std::vector<DirectedPath2> paths2G, std::vector<DirectedPath2> paths2H, std::vector<std::vector<GRBVar>> p) {

    GRBLinExpr costsum = 0;
    // Add the cost of mapping edge ij in G to edge kl in H 
    // and subtract the cost of deleting edge ij and inserting edge kl for all edge pairs.
    int ij = 0, kl;
    for (const auto &eG: boost::make_iterator_range(edges(G))) {
        auto edgeG = G[eG];
        kl = 0;
        for (const auto &eH: boost::make_iterator_range(edges(H))) {
            auto edgeH = H[eH];
            // The cost of mapping edge ij in G to edge kl in H is given by 
            // the difference of their weights (=lengths), weighted by the constant cLocal.
            costsum += ((abs(edgeG.weight - edgeH.weight)) * cLocal - (edgeG.weight + edgeH.weight) * cInsDel) * y[ij][kl];
            kl++;
        }
        ij++;
    }
    // Add the cost of mapping path ijh in G to path klm in H for all path pairs.
    int ijh = 0, klm;
    for (const auto &pG: paths2G) {
        auto [i, j, h, alphaijh] = pG;
        klm = 0;
        for (const auto &pH: paths2H) {
            auto [k, l, m, alphaklm] = pH;
            // The cost of mapping path ijh in G to path klm in H is given by 
            // the difference of their signed angles, weighted by the constant cGlobal.
            costsum += p[ijh][klm] * abs(alphaijh - alphaklm) * cGlobal;
            klm++;
        }
        ijh++;
    }
    return costsum;
}


// Computes the linear expression representing the objective cost based on the edge mapping variables z
// and the path mapping variables p. Here, the edges of G are directed and the for the edges of H we consider both directions. 
// The variables z[ij][kl] represent the mapping of directed edge ij in G to the directed edge kl in H and the variables z[ij][kl+1]
// represent the mapping of directed edge ij in G to the directed edge lk in H.
GRBLinExpr compute_objective_cost_z(const Graph_boost &G, const Graph_boost &H, 
        double cLocal, double cGlobal, double cInsDel, const std::vector<std::vector<GRBVar>> &z, 
        std::vector<DirectedPath2> paths2G, std::vector<DirectedPath2> paths2H, std::vector<std::vector<GRBVar>> p) {

    GRBLinExpr costsum = 0;
    // Add the cost of mapping directed edge ij in G to directed edges kl and lk in H
    // and subtract the cost of deleting directed edge ij and inserting directed edges kl and lk for all edge pairs.
    int ij = 0, kl;
    for (const auto &eG: boost::make_iterator_range(edges(G))) {
        auto edgeG = G[eG];
        kl = 0;
        for (const auto &eH: boost::make_iterator_range(edges(H))) {
            auto edgeH = H[eH];
            costsum += ((abs(edgeG.weight - edgeH.weight)) * cLocal - (edgeG.weight + edgeH.weight) * cInsDel) * z[ij][kl];
            costsum += ((abs(edgeG.weight - edgeH.weight)) * cLocal - (edgeG.weight + edgeH.weight) * cInsDel) * z[ij][kl+1];
            kl += 2;
        }
        ij++;
    }
    // Add the cost of mapping path ijh in G to path klm in H for all path pairs.
    int ijh = 0, klm;
    for (const auto &pG: paths2G) {
        auto [i, j, h, alphaijh] = pG;
        klm = 0;
        for (const auto &pH: paths2H) {
            auto [k, l, m, alphaklm] = pH;
            // The cost of mapping path ijh in G to path klm in H is given by 
            // the difference of their signed angles, weighted by the constant cGlobal.
            costsum += p[ijh][klm] * abs(alphaijh - alphaklm) * cGlobal;
            klm++;
        }
        ijh++;
    }
    return costsum;
}


// Print the solution of the ILP model, i.e., all variables with non-zero value.
void print_solution(GRBModel const &model, std::vector<std::vector<GRBVar>> const &x, 
        std::vector<std::vector<GRBVar>> const &y, std::vector<std::vector<GRBVar>> const &p, double gamma){
    std::cout << "Optimum modGED value: " << model.get(GRB_DoubleAttr_ObjVal) + gamma <<" with objective value " 
                    << model.get(GRB_DoubleAttr_ObjVal) << " and constant " << gamma << std::endl;
    auto n_g = x.size();
    auto n_h = x.front().size();
    auto m_g = y.size();
    auto m_h = y.front().size();
    auto p_g = p.size();
    auto p_h = p.front().size();

    for (int i = 0; i < n_g; i++) {
        for (int k = 0; k < n_h; k++) {
            if (x[i][k].get(GRB_DoubleAttr_X) > 0){
                std::cout << x[i][k].get(GRB_StringAttr_VarName) << " = " << x[i][k].get(GRB_DoubleAttr_X) << std::endl;
            }
        }
    }
    for (int ij = 0; ij < m_g; ij++) {
        for (int kl = 0; kl < m_h; kl++) {
            if (y[ij][kl].get(GRB_DoubleAttr_X) > 0){
                std::cout << y[ij][kl].get(GRB_StringAttr_VarName) << " = " << y[ij][kl].get(GRB_DoubleAttr_X) << std::endl;
            }
        }
    }
    for (int ijh = 0; ijh < p_g; ijh++) {
        for (int klm = 0; klm < p_h; klm++) {
            if (p[ijh][klm].get(GRB_DoubleAttr_X) > 0){
                std::cout << p[ijh][klm].get(GRB_StringAttr_VarName) << " = " << p[ijh][klm].get(GRB_DoubleAttr_X) << std::endl;
            }
        }
    }
}

} // end namespace


// Computes modGED using the ILP formulation F2+ as described in the paper 
// "Enhancing Graph Edit Distance Computation: Stronger and Orientation-based ILP Formulations"
// by Andrea D'Ascenzo, Julian Meffert, Petra Mutzel and Fabrizio Rossi (2025).
long double compute_modGED_F2Plus(const Graph_boost &G, const Graph_boost &H, 
        double cGlobal, double cLocal, double cInsDel, bool &timeout) {

    int nh = num_vertices(H);
    const int ng = num_vertices(G);
    const int eh = num_edges(H);
    const int eg = num_edges(G);
    auto const gamma = compute_gamma(G, H, cInsDel);

    auto timelimit = 60;

    try {
        GRBEnv env;
        
        env.set("OutputFlag", "0");

        GRBModel model = GRBModel(env);
        model.set(GRB_StringAttr_ModelName, "modGED_F2Plus");
        model.set(GRB_IntParam_Threads, 1);
        if (timelimit != 0) {
            model.set(GRB_DoubleParam_TimeLimit, timelimit);
        }

        model.set(GRB_IntParam_Seed, 1);
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

        // Adding binary node mapping variables with no contribution to objective
        std::vector<std::vector<GRBVar>> x(ng, std::vector<GRBVar>(nh));
        for (int i = 0; i < ng; i++) {
            for (int k = 0; k < nh; k++) {
                x[i][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + std::to_string(i) + "_" + std::to_string(k));
            }
        }
        
        // Storing the mapping from node pairs to edge index in G and H. 
        // Since the edges are undirected, we store the edge ij at both (i,j) and (j,i).
        std::vector<std::vector<int>> node_pair_to_yidx_G(ng, std::vector<int>(ng, -1));
        std::vector<std::vector<int>> node_pair_to_yidx_H(nh, std::vector<int>(nh, -1));

        // Adding binary edge mapping variables, here without contribution to objective.
        std::vector<std::vector<GRBVar>> y(eg, std::vector<GRBVar>(eh));
        int ij = 0, kl;
        for (auto edgeG: boost::make_iterator_range(edges(G))) {
            node_pair_to_yidx_G[edgeG.m_source][edgeG.m_target] = ij;
            node_pair_to_yidx_G[edgeG.m_target][edgeG.m_source] = ij;
            kl = 0;
            for (auto edgeH: boost::make_iterator_range(edges(H))) {
                if (ij == 0) {
                    node_pair_to_yidx_H[edgeH.m_source][edgeH.m_target] = kl;
                    node_pair_to_yidx_H[edgeH.m_target][edgeH.m_source] = kl;
                }
                y[ij][kl] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + std::to_string(edgeG.m_source) + std::to_string(edgeG.m_target) 
                                                                    + "_" + std::to_string(edgeH.m_source) + std::to_string(edgeH.m_target));
                kl++;
            }
            ij++;
        }

        // Find all undirected paths of length 2 in G and H.
        auto const undirected_paths_G = findAllPathsOfLength2(G);
        auto const undirected_paths_H = findAllPathsOfLength2(H);

        // Adding binary path mapping variables, here without contribution to objective. 
        // Add constraints linking path mapping variables to edge mapping variables, i.e.,
        // p[ijh][klm] = y[ij][kl] * y[jh][lm] for all paths ijh in G and klm in H.
        std::vector<std::vector<GRBVar>> p(undirected_paths_G.size(), std::vector<GRBVar>(undirected_paths_H.size()));
        for (size_t a = 0; a < undirected_paths_G.size(); a++) {
            auto const [i, j, h, alphaijh] = undirected_paths_G[a];
            for (size_t b = 0; b < undirected_paths_H.size(); b++) {
                auto [k, l, m, alphaklm] = undirected_paths_H[b];
                p[a][b] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "p_" + std::to_string(i) + std::to_string(j) + std::to_string(h) 
                                                                  + "_" + std::to_string(k) + std::to_string(l) + std::to_string(m));
                model.addConstr(p[a][b] >=
                                y[node_pair_to_yidx_G[i][j]][node_pair_to_yidx_H[k][l]] + y[node_pair_to_yidx_G[j][h]][node_pair_to_yidx_H[l][m]] - 1,
                                "Path Edge Relationship"); 
                model.addConstr(p[a][b] <= y[node_pair_to_yidx_G[i][j]][node_pair_to_yidx_H[k][l]], "Path UB1");
                model.addConstr(p[a][b] <= y[node_pair_to_yidx_G[j][h]][node_pair_to_yidx_H[l][m]], "Path UB2");                
            }
        }

        // Add constraints ensuring a valid node mapping, i.e., 
        // each node in G is mapped to at most one node in H and vice versa.
        for (int i = 0; i < ng; i++){
            model.addConstr(vecSum(x[i]) <= 1, "Vertex Mapping G");
        }
        for (int k = 0; k < nh; k++){
            model.addConstr(colSum(x, k) <= 1, "Vertex Mapping H");
        }

        // Add topological constraints linking node mapping variables to edge mapping variables.
        for (auto const &edgeG : boost::make_iterator_range(edges(G))){    
            auto const i = edgeG.m_source;
            auto const j = edgeG.m_target;        
            for (int k = 0; k < nh; k++){                    
                GRBLinExpr sum = 0;                
                for (auto const l : boost::make_iterator_range(adjacent_vertices(k, H))){
                    sum += y[node_pair_to_yidx_G[i][j]][node_pair_to_yidx_H[k][l]];
                }
                model.addConstr(sum <= x[i][k] + x[j][k], "Topological Constraints H");              
            }
        }
        for (auto const &edgeH : boost::make_iterator_range(edges(H))){
            auto const k = edgeH.m_source;
            auto const l = edgeH.m_target;
            for (int i = 0; i < ng; i++){
                GRBLinExpr sum = 0;
                for (auto const j : boost::make_iterator_range(adjacent_vertices(i, G))){
                    sum += y[node_pair_to_yidx_G[i][j]][node_pair_to_yidx_H[k][l]];
                }
                model.addConstr(sum <= x[i][k] + x[i][l], "Topological Constraints G");
            }
        }

        // Set objective function as defined in the function compute_objective_cost_y plus the constant gamma.
        model.setObjective(compute_objective_cost_y(G, H, cGlobal, cLocal, cInsDel, y, undirected_paths_G, undirected_paths_H, p) + gamma,
                           GRB_MINIMIZE);

        // model.update();
        // model.write("check_lp_f2plus.lp");

        model.optimize();
        model.update();
        // print_solution(model, x, y, p, gamma, G, H);

        auto status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL) {
            return model.get(GRB_DoubleAttr_ObjVal);
        } 
        else if (status == GRB_TIME_LIMIT) {
            std::cout << "Time limit reached, no proven optimal solution\n";
            timeout = true;
        } else if (status == GRB_INFEASIBLE) {
            std::cout << "Model is infeasible\n";
        } else if (status == GRB_INF_OR_UNBD) {
            std::cout << "Model is infeasible or unbounded\n";
        } else if (status == GRB_UNBOUNDED) {
            std::cout << "Model is unbounded\n";
        } else {
            std::cout << "Optimization stopped with status = " << status << "\n";
        }
        return -1;

    } catch (GRBException &e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        return -1;
    }
}


// Computes modGED using the ILP formulation FORI as described in the paper 
// "Enhancing Graph Edit Distance Computation: Stronger and Orientation-based ILP Formulations"
// by Andrea D'Ascenzo, Julian Meffert, Petra Mutzel and Fabrizio Rossi (2025).
long double compute_modGED_FORI(const Graph_boost &G, const Graph_boost &H, 
        double cGlobal, double cLocal, double cInsDel, bool &timeout) {

    int nh = num_vertices(H);
    const int ng = num_vertices(G);
    const int eh = num_edges(H);
    const int eg = num_edges(G);
    const double gamma = compute_gamma(G, H, cInsDel);

    auto timelimit = 60;

    try {
        GRBEnv env;
        
        env.set("OutputFlag", "0");

        GRBModel model = GRBModel(env);
        model.set(GRB_StringAttr_ModelName, "modGED_FORI");
        model.set(GRB_IntParam_Threads, 1);
        if (timelimit != 0) {
            model.set(GRB_DoubleParam_TimeLimit, timelimit);
        }

        model.set(GRB_IntParam_Seed, 1);
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

        // Adding binary node mapping variables with no contribution to objective
        std::vector<std::vector<GRBVar>> x(ng, std::vector<GRBVar>(nh));
        for (int i = 0; i < ng; i++) {
            for (int k = 0; k < nh; k++) {
                x[i][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + std::to_string(i) + "_" + std::to_string(k));
            }
        }

        // Storing the mapping from node pairs to edge index in G and H.
        // We consider the edges in G as dierected and for each edge in H we consider both directions.
        std::vector<std::vector<int>> node_pair_to_zidx_G(ng, std::vector<int>(ng, -1));
        std::vector<std::vector<int>> node_pair_to_zidx_H(nh, std::vector<int>(nh, -1));

        // Adding binary edge mapping variables, here without contribution to objective.
        // We need 2*eh variables for each edge in G, 
        // and we consecutively add two variables for mapping the edge ij in G to the edges kl and lk in H.
        std::vector<std::vector<GRBVar>> z(eg, std::vector<GRBVar>(2*eh));
        int ij = 0, kl;
        for (auto edgeG: boost::make_iterator_range(edges(G))) {
            node_pair_to_zidx_G[edgeG.m_source][edgeG.m_target] = ij;            
            kl = 0;
            for (auto edgeH: boost::make_iterator_range(edges(H))) {
                if (ij == 0) {
                    node_pair_to_zidx_H[edgeH.m_source][edgeH.m_target] = kl;
                    node_pair_to_zidx_H[edgeH.m_target][edgeH.m_source] = kl + 1;
                }
                z[ij][kl] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z_" + std::to_string(edgeG.m_source) + std::to_string(edgeG.m_target) 
                                                                    + "_" + std::to_string(edgeH.m_source) + std::to_string(edgeH.m_target));
                z[ij][kl+1] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z_" + std::to_string(edgeG.m_source) + std::to_string(edgeG.m_target) 
                                                                      + "_" + std::to_string(edgeH.m_target) + std::to_string(edgeH.m_source));
                kl += 2;
            }
            ij++;
        }

        // Find all directed paths of length 2 in G and H. 
        // Here, we consider the edges to be undirected in both graphs.
        auto directed_paths_G = findAllPathsOfLength2(G);
        auto directed_paths_H = findAllPathsOfLength2(H);

        // Adding binary path mapping variables, here without contribution to objective.
        // Add constraints linking path mapping variables to edge mapping variables, i.e.,
        // p[ijh][klm] = z[ij][kl] * z[jh][lm] + z[ij][lk] * z[jh][ml] for all directed paths ijh in G and klm in H.
        std::vector<std::vector<GRBVar>> p(directed_paths_G.size(), std::vector<GRBVar>(directed_paths_H.size()));
        for (size_t a = 0; a < directed_paths_G.size(); a++) {
            auto [i, j, h, alphaijh] = directed_paths_G[a];
            for (size_t b = 0; b < directed_paths_H.size(); b++) {
                auto [k, l, m, alphaklm] = directed_paths_H[b];
                p[a][b] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "p_" + std::to_string(i) + std::to_string(j) + std::to_string(h) 
                                                                  + "_" + std::to_string(k) + std::to_string(l) + std::to_string(m));

                // Since the paths do not care about orientation, we have to make sure that we get the edge in G in the right direction.
                auto get_directed_edge_G = [&](auto const i, auto const j){
                    if (node_pair_to_zidx_G[i][j] != -1){
                        return node_pair_to_zidx_G[i][j];
                    }
                    else if (node_pair_to_zidx_G[j][i] != -1){
                        return node_pair_to_zidx_G[j][i];
                    }
                    else {
                        std::cerr << "Error: For j adjacent to i, we stored neither the edge ij nor the edge ji." << std::endl;
                        return -1;
                    }       
                };
                model.addConstr(p[a][b] 
                                >= z[get_directed_edge_G(i,j)][node_pair_to_zidx_H[k][l]] + z[get_directed_edge_G(j,h)][node_pair_to_zidx_H[l][m]]
                                + z[get_directed_edge_G(i,j)][node_pair_to_zidx_H[l][k]] + z[get_directed_edge_G(j,h)][node_pair_to_zidx_H[m][l]] - 1,
                            "Path Edge Relation");    
                model.addConstr(p[a][b] <= z[get_directed_edge_G(i,j)][node_pair_to_zidx_H[k][l]] + z[get_directed_edge_G(i,j)][node_pair_to_zidx_H[l][k]], "Path UB1");
                model.addConstr(p[a][b] <= z[get_directed_edge_G(j,h)][node_pair_to_zidx_H[l][m]] + z[get_directed_edge_G(j,h)][node_pair_to_zidx_H[m][l]], "Path UB2");            
            }
        }

        // Add constraints ensuring a valid node mapping, i.e., 
        // each node in G is mapped to at most one node in H and vice versa.
        for (int i = 0; i < ng; i++){
            model.addConstr(vecSum(x[i]) <= 1, "Vertex Mapping G");
        }
        for (int k = 0; k < nh; k++){
            model.addConstr(colSum(x, k) <= 1, "Vertex Mapping H");
        }

        // Add topological constraints linking node mapping variables to edge mapping variables.
        for (auto const &edgeG : boost::make_iterator_range(edges(G))){    
            auto const i = edgeG.m_source;
            auto const j = edgeG.m_target;        
            for (int k = 0; k < nh; k++){                    
                GRBLinExpr sum = 0;    
                // Outgoing edges            
                for (auto const l : boost::make_iterator_range(adjacent_vertices(k, H))){
                    sum += z[node_pair_to_zidx_G[i][j]][node_pair_to_zidx_H[k][l]];
                }
                model.addConstr(sum <= x[i][k], "Outgoing Edge Sum");

                sum = 0;
                // Ingoing edges
                for (auto const l : boost::make_iterator_range(adjacent_vertices(k, H))){
                    sum += z[node_pair_to_zidx_G[i][j]][node_pair_to_zidx_H[l][k]];
                }
                model.addConstr(sum <= x[j][k], "Ingoing Edge Sum");              
            }
        }
        for (auto const &edgeH : boost::make_iterator_range(edges(H))){
            auto const k = edgeH.m_source;
            auto const l = edgeH.m_target;
            for (int i = 0; i < ng; i++){
                GRBLinExpr sum = 0;
                for (auto const j : boost::make_iterator_range(adjacent_vertices(i, G))){
                    if (node_pair_to_zidx_G[i][j] != -1){
                        sum += z[node_pair_to_zidx_G[i][j]][node_pair_to_zidx_H[k][l]];
                    }
                    else if (node_pair_to_zidx_G[j][i] != -1){
                        sum += z[node_pair_to_zidx_G[j][i]][node_pair_to_zidx_H[l][k]];
                    }
                    else {
                        std::cerr << "Error: For j adjacent to i, we stored neither the edge ij nor the edge ji." << std::endl;
                    }                    
                }
                model.addConstr(sum <= x[i][k], "Edge Sum H");
            }
        }

        // Set objective function as defined in the function compute_objective_cost_z plus the constant gamma.
        model.setObjective(compute_objective_cost_z(G, H, cGlobal, cLocal, cInsDel, z, directed_paths_G, directed_paths_H, p) + gamma,
                           GRB_MINIMIZE);

        // model.update();
        // model.write("check_lp_fori.lp");

        model.optimize();
        // print_solution(model, x, z, p, gamma, G, H);

        auto status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL) {
            return model.get(GRB_DoubleAttr_ObjVal);
        } 
        else if (status == GRB_TIME_LIMIT) {
            std::cout << "Time limit reached, no proven optimal solution\n";
            timeout = true;
        } else if (status == GRB_INFEASIBLE) {
            std::cout << "Model is infeasible\n";
        } else if (status == GRB_INF_OR_UNBD) {
            std::cout << "Model is infeasible or unbounded\n";
        } else if (status == GRB_UNBOUNDED) {
            std::cout << "Model is unbounded\n";
        } else {
            std::cout << "Optimization stopped with status = " << status << "\n";
        }
        return -1;

    } catch (GRBException &e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        return -1;
    }
}
