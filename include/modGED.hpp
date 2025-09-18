#ifndef MODGED_HPP
#define MODGED_HPP

#include "graph.hpp"

#include "gurobi_c++.h"



// using GeometricGraph = std::pair<graph<Point, EdgeLength>, AngleData>;
using Coordinate = double;
using Point = std::pair<Coordinate, Coordinate>;
using EdgeLength = double;

EdgeLength compute_euclidean_distance(Point p, Point q){
    return sqrt(pow(abs(p.first - q.first), 2) + pow(abs(p.second - q.second), 2));
}

template <typename T, typename U>
double compute_modfied_GED(const graph<T, U> &g, const graph<T, U> &h, double crot, double cmod)
{
    int n_h = h.number_of_nodes();
    int n_g = g.number_of_nodes();
    int m_h = h.number_of_edges();
    int m_g = g.number_of_edges();

    try{
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "gurobi.log");
        env.set("OutputFlag", "0");
        env.start();

        GRBModel model = GRBModel(env);
        // double gammaVar = gamma(G, H, cInsDel);

        GRBVar node_sub[n_g][n_h];
        GRBVar edge_sub[m_g][m_h];

        int num_paths_of_length2_g = 0;
        int num_paths_of_length2_h = 0;

        for (int i = 0; i < n_g; i++){
            num_paths_of_length2_g += (g.get_degree(i) * (g.get_degree(i) - 1)) / 2;
            for (int k = 0; k < n_h; k++){
                node_sub[i][k] = model.addVar(0, 1, 0, GRB_BINARY, "x_" + std::to_string(i) + "_" + std::to_string(k));
                if (i == 0){
                    num_paths_of_length2_h += (h.get_degree(k) * (h.get_degree(k) - 1)) / 2;
                }
            }
        }

        for (int ij = 0; ij < m_g; ij++){
            for (int kl = 0; kl < m_h; kl++){
                auto cost = cmod * (std::abs(g.get_edge_label(ij) - h.get_edge_label(kl))
                            - g.get_edge_label(ij) - h.get_edge_label(kl));
                edge_sub[ij][kl] = model.addVar(0, 1, cost, GRB_BINARY, "y_" + 
                                                std::to_string(g.get_edge(ij).first) + std::to_string(g.get_edge(ij).second) 
                                                + "_" + std::to_string(h.get_edge(kl).first) + std::to_string(h.get_edge(kl).second));
            }
        }

        auto get_angle = [&] (const graph<T, U> &g, node i, node j, node n) {

            assert(g.has_edge(i, j));
            assert(g.has_edge(j, n));

            Point a = g.get_node_label(i);
            Point v = g.get_node_label(j);
            Point b = g.get_node_label(n);

            double d1x = v.first - a.first;
            double d1y = v.second - a.second;

            double d2x = v.first - b.first;
            double d2y = v.second - b.second;

            // Scalarproduct
            double dotProduct = d1x * d2x + d1y * d2y;

            // vector length
            double norm1 = std::sqrt(d1x * d1x + d1y * d1y);
            double norm2 = std::sqrt(d2x * d2x + d2y * d2y);

            // Angle calcualtion
            double cosTheta = dotProduct / (norm1 * norm2);
            double theta = std::acos(cosTheta);

            // Umwandlung in Grad
            //return theta * (180.0 / M_PI);

            // Kreuzprodukt berechnen
            double cross = d1x * d2y - d1y * d2x;

            // Choose the right andle (from (a,v)'s point of view)
            if (cross < 0) {
                theta = 2 * M_PI - theta;
            }
            return theta;
        };

        GRBVar path_sub[num_paths_of_length2_g][num_paths_of_length2_h];

        idx path_counter_g = 0;
        for (node j = 0; j < n_g; j++){
            auto neighbors_j = g.get_neighbors(j);
            for (idx idx1 = 0; idx1 + 1 < neighbors_j.size(); idx1++){
                for (idx idx2 = idx1 + 1; idx2 < neighbors_j.size(); idx2++){
                    node i = neighbors_j[idx1];
                    node n = neighbors_j[idx2];

                    auto angle_g = get_angle(g, i, j, n);

                    idx path_counter_h = 0;

                    for (node l = 0; l < n_h; l++){
                        auto neighbors_l = h.get_neighbors(l);
                        for (idx idx3 = 0; idx3 < neighbors_l.size(); idx3++){
                            for (idx idx4 = idx3 + 1; idx4 < neighbors_l.size(); idx4++){
                                node k = neighbors_l[idx3];
                                node m = neighbors_l[idx4];

                                auto angle_h = get_angle(h, k, l, m);
                                auto cost = crot * std::abs(angle_h - angle_g);

                                path_sub[path_counter_g][path_counter_h] = model.addVar(0, 1, cost, GRB_BINARY, "p_" + std::to_string(i)
                                                                + std::to_string(j) + std::to_string(n) + "_" + std::to_string(k)
                                                                + std::to_string(l) + std::to_string(m));
                                
                                model.addConstr(path_sub[path_counter_g][path_counter_h] <= edge_sub[g.get_edge_id(i, j)][h.get_edge_id(k, l)]);
                                model.addConstr(path_sub[path_counter_g][path_counter_h] <= edge_sub[g.get_edge_id(j, n)][h.get_edge_id(l, m)]);
                                model.addConstr(path_sub[path_counter_g][path_counter_h] >= edge_sub[g.get_edge_id(i, j)][h.get_edge_id(k, l)] + edge_sub[g.get_edge_id(j, n)][h.get_edge_id(l, m)] - 1);
                                model.addConstr(path_sub[path_counter_g][path_counter_h] >= edge_sub[g.get_edge_id(i, j)][h.get_edge_id(l, m)] + edge_sub[g.get_edge_id(j, n)][h.get_edge_id(k, l)] - 1);

                                path_counter_h++;
                            }
                        }
                    }
                    path_counter_g++;
                }
            }
        }

        GRBLinExpr le;
        for (int i = 0; i < n_g; ++i){
            le = 0;
            for (int k = 0; k < n_h; ++k)
                le += node_sub[i][k];            
            model.addConstr(le <= 1, "Node_Sub");
        }

        GRBLinExpr le1;
        for (int k = 0; k < n_h; k++) {
            le1 = 0;
            for (int i = 0; i < n_g; i++)
                le1 += node_sub[i][k];
            model.addConstr(le1, GRB_LESS_EQUAL, 1, "Node_Sub");
        }

        GRBLinExpr le2;
        GRBLinExpr le3;
        for (int ij = 0; ij < m_g; ++ij){
            for (int k = 0; k < n_h; ++k){
                le2 = 0;
                le3 = node_sub[g.get_edge(ij).first][k] + node_sub[g.get_edge(ij).second][k];
                for (int l : h.get_neighbors(k)){
                    int kl = h.get_edge_id(k,l);
                    le2 += edge_sub[ij][kl];
                }
                model.addConstr(le2, GRB_LESS_EQUAL, le3, "Topological_(" + std::to_string(g.get_edge(ij).first) + "," +
                        std::to_string(g.get_edge(ij).second) + ")" + "_" + std::to_string(k));
            }
        }

        GRBLinExpr le4;
        GRBLinExpr le5;
        for (int i = 0; i < n_g; ++i){
            for (int kl = 0; kl < m_h; ++kl){
                le5 = 0;
                le4 = node_sub[i][h.get_edge(kl).first] + node_sub[i][h.get_edge(kl).second];
                for (int j : g.get_neighbors(i)){
                    int ij = g.get_edge_id(i,j);
                    le5 += edge_sub[ij][kl];
                }
                model.addConstr(le5, GRB_LESS_EQUAL, le4, "Topological_" + std::to_string(i) + "_(" + std::to_string(h.get_edge(kl).first)
                        + "," + std::to_string(h.get_edge(kl).second) + ")");
            }
        }

        model.update();

        model.optimize();

        double constant = 0;
        for (int ij = 0; ij < m_g; ij++){
            constant += g.get_edge_label(ij);
        }
        for (int kl = 0; kl < m_h; kl++){
            constant += h.get_edge_label(kl);
        }
        constant = constant * cmod;

        double objective_value = model.get(GRB_DoubleAttr_ObjVal);

        std::cout << "ged val: " << model.get(GRB_DoubleAttr_ObjVal) + constant <<" obj val: " 
                    << model.get(GRB_DoubleAttr_ObjVal) << " constant: " << constant << std::endl;

        std::vector<bool> inserted_k = std::vector<bool>(n_h, true);
        for (int i = 0; i < n_g; i++) {
            bool deleted_i = true;
            for (int k = 0; k < n_h; k++) {
                if (node_sub[i][k].get(GRB_DoubleAttr_X) > 0){
                    deleted_i = false;
                    inserted_k[k] = false;
                    std::cout << "node " << i << " to " << k << " value: " << node_sub[i][k].get(GRB_DoubleAttr_X) 
                    << " objective coefficient: " << node_sub[i][k].get(GRB_DoubleAttr_Obj) << std::endl;
                }
            }
            if (deleted_i)
                std::cout << "deleted node " << i << std::endl;
        }
        for (int k = 0; k < n_h; ++k){
            if (inserted_k[k])
                std::cout << "inserted node " << k << std::endl;
        }
        for (int ij = 0; ij < m_g; ij++) {
            for (int kl = 0; kl < m_h; kl++) {
                if (edge_sub[ij][kl].get(GRB_DoubleAttr_X) > 0){
                    std::cout << "edge " << g.get_edge(ij).first << g.get_edge(ij).second << " to " 
                    << h.get_edge(kl).first << h.get_edge(kl).second 
                    << " value: " << edge_sub[ij][kl].get(GRB_DoubleAttr_X) 
                    << " objective coefficient: " << edge_sub[ij][kl].get(GRB_DoubleAttr_Obj) << std::endl;
                }
            }
        }

        std::vector<double>ij_sum = std::vector<double>(m_g, 0);
        std::vector<double>kl_sum = std::vector<double>(m_h, 0);
        for (int ij = 0; ij < m_g; ij++) {
            for (int kl = 0; kl < m_h; kl++) {
                ij_sum[ij] += edge_sub[ij][kl].get(GRB_DoubleAttr_X);
                kl_sum[kl] += edge_sub[ij][kl].get(GRB_DoubleAttr_X);
            }
        }
        for (int ij = 0; ij < m_g; ij++) {
            if(ij_sum[ij] == 0)
                std::cout << "deleted edge " << g.get_edge(ij).first << g.get_edge(ij).second << std::endl;
        }
        for (int kl = 0; kl < m_h; kl++) {
            if(kl_sum[kl] == 0)
                std::cout << "inserted edge " << h.get_edge(kl).first << h.get_edge(kl).second << std::endl;
        }

        idx counter_g = 0;
        for (node j = 0; j < n_g; j++){
            auto neighbors_j = g.get_neighbors(j);
            for (idx idx1 = 0; idx1 + 1 < neighbors_j.size(); idx1++){
                for (idx idx2 = idx1 + 1; idx2 < neighbors_j.size(); idx2++){
                    node i = neighbors_j[idx1];
                    node n = neighbors_j[idx2];

                    idx counter_h = 0;
                    for (node l = 0; l < n_h; l++){
                        auto neighbors_l = h.get_neighbors(l);
                        for (idx idx3 = 0; idx3 < neighbors_l.size(); idx3++){
                            for (idx idx4 = idx3 + 1; idx4 < neighbors_l.size(); idx4++){
                                node k = neighbors_l[idx3];
                                node m = neighbors_l[idx4];
                                if (path_sub[counter_g][counter_h].get(GRB_DoubleAttr_X) > 0){
                                    std::cout << "Path " << std::to_string(i) << std::to_string(j) << std::to_string(n) 
                                              << " to " << std::to_string(k) << std::to_string(l) << std::to_string(m)
                                              << " value: " << path_sub[counter_g][counter_h].get(GRB_DoubleAttr_X) 
                                              << " objective coefficient: " << path_sub[counter_g][counter_h].get(GRB_DoubleAttr_Obj) 
                                              << std::endl;
                                }
                                counter_h++;
                            }
                        }
                    }
                    counter_g++;
                }
            }
        }

        //const auto output = IO::create_output(model, options());
        // delete &model;
        // delete &env;

        return objective_value;

    }     
    catch (GRBException &e) {
        std::cerr << "Gurobi Error: " << e.getMessage() << std::endl;
        return -1;
    }
}
#endif // MODGED_HPP