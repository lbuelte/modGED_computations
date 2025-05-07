#include "../include/distance.h"


GRBLinExpr insdel(const std::vector<GRBVar> &theta, const Graph_boost &G, const Graph_boost &H) {
    GRBLinExpr insdelsum = 0;
    int index = 0;

    for (const auto &e: boost::make_iterator_range(edges(G))) {
        insdelsum += G[e].weight * theta[index];
        index++;
    }


    for (const auto &e: boost::make_iterator_range(edges(H))) {
        insdelsum += H[e].weight * theta[index];
        index++;
    }
    //std::cout << "testpoint 2.1.2" << std::endl;
    return insdelsum;
}

double angleBetweenLines(const Point_boost &a, const Point_boost &v, const Point_boost &b) {

    /*std::cout << "doofe Punkte: " << a.a << "  " << a.b << " , " << v.a << "  " << v.b << " , " << b.a << "  " << b.b
              << std::endl;*/
    // Vectors
    double d1x = v.a - a.a;
    double d1y = v.b - a.b;

    double d2x = v.a - b.a;
    double d2y = v.b - b.b;

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
}


GRBLinExpr colSum(const std::vector<std::vector<GRBVar>> &vars, int index) {
    GRBLinExpr sum = 0;
    for (const auto &vec: vars) sum += vec[index];
    return sum;
}

GRBLinExpr vecSum(const std::vector<GRBVar> &vec) {
    GRBLinExpr sum = 0;
    for (const auto &val: vec) sum += val;
    return sum;
}

GRBLinExpr anglediff_old(const Graph_boost &G, const Graph_boost &H,
                         std::pair<Vertex_boost, Vertex_boost> edgeG,
                         boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long> edgeH,
                         const std::vector<std::vector<GRBVar>> &y,
                         const std::vector<std::vector<GRBVar>> &x,
                         const std::vector<boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long>> yDecodeH,
                         const std::vector<std::vector<int>> yDecodeG) {

    GRBLinExpr angeldiffSum = 0;
    Vertex_boost v, u, z, vH, uH, zH;
    Point_boost vp, up, zp, vHp, uHp, zHp;

    v = edgeG.first;
    u = edgeG.second;
    vH = edgeH.m_source;
    uH = edgeH.m_target;
    vp = G[v].p;
    up = G[u].p;
    vHp = H[vH].p;
    uHp = H[uH].p;
    int i;

    for (const auto &neighbor: boost::make_iterator_range(out_edges(u, G))) {
        if ((neighbor.m_source == v && neighbor.m_target != u) || (neighbor.m_source != u && neighbor.m_target == v)) {
            i = yDecodeG[neighbor.m_source][neighbor.m_target];
            z = neighbor.m_source == u ? neighbor.m_target : neighbor.m_source;
            zp = G[z].p;
            double angle1 = angleBetweenLines(vp, up, zp);
            for (int j = 0; j < num_edges(H); j++) {
                auto neighborEq = yDecodeH[j];
                zH = neighborEq.m_source;
                if (edgeH.m_source == neighborEq.m_source) {
                    zH = neighborEq.m_target;
                }
                double angle2 = angleBetweenLines(vHp, uHp, H[zH].p);

                angeldiffSum += abs(angle1 - angle2) * y[i][j];
            }
        }
    }

    return angeldiffSum;
}

GRBLinExpr anglediff(const Graph_boost &G, const Graph_boost &H,
                     std::pair<Vertex_boost, Vertex_boost> edgeG,
                     boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long> edgeH,
                     const std::vector<std::vector<GRBVar>> &y,
                     const std::vector<std::vector<GRBVar>> &x,
                     const std::vector<boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long>> yDecodeH,
                     const std::vector<std::vector<int>> yDecodeG) {

    GRBLinExpr angeldiffSum = 0;
    Vertex_boost v, u, z, vH, uH, zH;
    Point_boost vp, up, zp, vHp, uHp, zHp;

    v = edgeG.first;
    u = edgeG.second;
    vH = edgeH.m_source;
    uH = edgeH.m_target;
    vp = G[v].p;
    up = G[u].p;
    vHp = H[vH].p;
    uHp = H[uH].p;
    int i;

    for (const auto &neighbor: boost::make_iterator_range(out_edges(u, G))) {
        if ((neighbor.m_source == u && neighbor.m_target != v) || (neighbor.m_source != v && neighbor.m_target == u)) {
            i = yDecodeG[neighbor.m_source][neighbor.m_target];
            z = neighbor.m_source == u ? neighbor.m_target : neighbor.m_source;
            zp = G[z].p;
            double angle1 = angleBetweenLines(vp, up, zp);
            for (int j = 0; j < num_edges(H); j++) {
                auto neighborEq = yDecodeH[j];
                zH = neighborEq.m_source;
                if (edgeH == neighborEq || !(neighborEq.m_source == uH || neighborEq.m_target == uH)) continue;
                if (uH == neighborEq.m_source) {
                    zH = neighborEq.m_target;
                }
                double angle2 = angleBetweenLines(vHp, uHp, H[zH].p);

                angeldiffSum += abs(angle1 - angle2) * y[i][j];
            }
        }
    }

    return angeldiffSum;
}

GRBLinExpr anglediff(const Graph_boost &G, const Graph_boost &H,
                     std::pair<Vertex_boost, Vertex_boost> edgeG,
                     boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long> edgeH,
                     const std::vector<std::vector<int>> &y,
                     const std::vector<std::vector<GRBVar>> &x,
                     const std::vector<boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long>> yDecodeH,
                     const std::vector<std::vector<int>> yDecodeG) {

    GRBLinExpr angeldiffSum = 0;
    Vertex_boost v, u, z, vH, uH, zH;
    Point_boost vp, up, zp, vHp, uHp, zHp;

    v = edgeG.first;
    u = edgeG.second;
    vH = edgeH.m_source;
    uH = edgeH.m_target;
    vp = G[v].p;
    up = G[u].p;
    vHp = H[vH].p;
    uHp = H[uH].p;
    int i;
    /*std::cout << "Das kommt an1: (" << vp.a << ", " << vp.b << "), (" << up.a << ", " << up.b << ") mapped to " << "("
              << vHp.a << ", " << vHp.b << "), (" << uHp.a << ", " << uHp.b << ")"
              << std::endl;*/
    for (auto neighbor: boost::make_iterator_range(out_edges(u, G))) {
        //std::cout << "So ein Müll " << std::endl;
        if ((neighbor.m_source == u && neighbor.m_target != v) || (neighbor.m_source != v && neighbor.m_target == u)) {
            //std::cout << "So ein Müll 2" << std::endl;
            i = yDecodeG[neighbor.m_source][neighbor.m_target];
            z = neighbor.m_source;
            if (neighbor.m_source == u) {
                z = neighbor.m_target;
            }
            zp = G[z].p;
            for (int j = 0; j < num_edges(H); j++) {
                auto neighborEq = yDecodeH[j];
                zH = neighborEq.m_source;
                if (edgeH == neighborEq || !(neighborEq.m_source == uH || neighborEq.m_target == uH)) continue;
                if (uH == neighborEq.m_source) {
                    zH = neighborEq.m_target;
                }
                /*std::cout << "Das kommt an2: (" << vp.a << ", " << vp.b << "), (" << up.a << ", " << up.b << "), ("
                          << zp.a << ", " << zp.b << ") mapped to " << "(" << vHp.a << ", " << vHp.b << "), (" << uHp.a
                          << ", " << uHp.b
                          << "), ("
                          << H[zH].p.a << ", " << H[zH].p.b << ")"
                          << std::endl;
                std::cout << "aka2: " << v << ", " << u << ", " << z << " to " << vH << ", "
                          << uH << ", " << zH << std::endl;
                std::cout << "G sagt: (" << edgeG.first << ", " << edgeG.second << "), (" << neighbor.m_source << ", "
                          << neighbor.m_target
                          << ") " << std::endl;
                std::cout << "H sagt: (" << edgeH.m_source << ", " << edgeH.m_target << "), (" << neighborEq.m_source
                          << ", " << neighborEq.m_target
                          << ") " << std::endl;*/
                double angle1 = angleBetweenLines(vp, up, zp);
                double angle2 = angleBetweenLines(vHp, uHp, H[zH].p);

                //std::cout << "winkel: " << angle1 << "  " << angle2 << std::endl;

                angeldiffSum += abs(angle1 - angle2) * y[i][j];
            }
        }
    }

    return angeldiffSum;
}

GRBQuadExpr
mapcost(const std::vector<std::vector<GRBVar>> &y, const std::vector<std::vector<GRBVar>> &x, const Graph_boost &G,
        const Graph_boost &H, double cLocal,
        double cGlobal, std::vector<std::vector<int>> yDecodeG,
        std::vector<boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long>> yDecodeH) {
    GRBQuadExpr costsum = 0;
    int i = 0, j;
    //std::cout << "testpoint 2.2\n";
    for (const auto &eG: boost::make_iterator_range(edges(G))) {
        auto edgeG = G[eG];
        j = 0;
        for (const auto &eH: boost::make_iterator_range(edges(H))) {
            auto edgeH = H[eH];
            costsum += abs(edgeG.weight - edgeH.weight) * cLocal * y[i][j];
            costsum += abs(edgeG.weight - edgeH.weight) * cLocal * y[i + 1][j];
            //std::cout << edgeG.weight << std::endl;
            //std::cout << edgeH.weight << std::endl;
            //std::cout << eG.m_source << " " << eG.m_target << " mapped to " << eH.m_source << " " << eH.m_target<< std::endl;

            auto angeldiff1 = anglediff(G, H, {eG.m_source, eG.m_target}, eH, y, x, yDecodeH, yDecodeG),
                    angeldiff2 = anglediff(G, H, {eG.m_target, eG.m_source}, eH, y, x, yDecodeH, yDecodeG);

            //std::cout << "dieser kack winkel " << angeldiff1 << " " << angeldiff2 << std::endl;

            costsum += (angeldiff1 + angeldiff2) * cGlobal * y[i][j];
            costsum += (angeldiff1 + angeldiff2) * cGlobal * y[i + 1][j];
            j++;
        }
        i += 2;
    }

    return costsum;
}

void printResult(GRBModel &model) {
    // Prüfen, ob eine optimale Lösung gefunden wurde
    int status = model.get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL) {
        double objVal = model.get(GRB_DoubleAttr_ObjVal);  // Optimalen Zielfunktionswert holen
        std::cout << "Optimal objective value: " << objVal << std::endl;


    } else {
        std::cout << "No optimal Solution found. Status: " << status << std::endl;
        if (status == GRB_INFEASIBLE) {
            std::cout << "Model is infeasible. Compute IIS..." << std::endl;
            model.computeIIS();
            model.write("infeasible.ilp");  // Speichere das IIS zur Analyse
        } else if (status == GRB_UNBOUNDED) {
            std::cout << "Model is not bounded." << std::endl;
        }
    }
}

void Testtest(const Graph_boost &G, const Graph_boost &H, const std::vector<std::vector<GRBVar>> &y,
              const std::vector<std::vector<GRBVar>> &x, double cLocal, double cGlobal,
              const std::vector<boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long>> yDecodeH,
              const std::vector<std::vector<int>> yDecodeG) {
    GRBQuadExpr costsum = 0;
    /*std::vector<std::vector<int>> yTest =
            {{1, 0, 0, 0, 0},
             {0, 0, 0, 0, 0},
             {0, 1, 0, 0, 0},
             {0, 0, 0, 0, 0},
             {0, 0, 1, 0, 0},
             {0, 0, 0, 0, 0},
             {0, 0, 0, 1, 0},
             {0, 0, 0, 0, 0},
             {0, 0, 0, 0, 1},
             {0, 0, 0, 0, 0},
             {0, 0, 0, 0, 0},
             {0, 0, 0, 0, 0}
            };*/
    std::vector<std::vector<int>> yTest =
            {{0, 0, 0, 0, 0},
             {0, 0, 0, 0, 0},
             {1, 0, 0, 0, 0},
             {0, 0, 0, 0, 0},
             {0, 0, 1, 0, 0},
             {0, 0, 0, 0, 0},
             {0, 0, 0, 1, 0},
             {0, 0, 0, 0, 0},
             {0, 0, 0, 0, 1},
             {0, 0, 0, 0, 0},
             {0, 0, 0, 0, 0},
             {0, 1, 0, 0, 0}
            };
    int i = 0, j;
    for (auto eG: boost::make_iterator_range(edges(G))) {
        auto edgeG = G[eG];
        j = 0;
        for (auto eH: boost::make_iterator_range(edges(H))) {
            auto edgeH = H[eH];
            costsum += abs(edgeG.weight - edgeH.weight) * cLocal * yTest[i][j];
            costsum += abs(edgeG.weight - edgeH.weight) * cLocal * yTest[i + 1][j];
            std::cout << edgeG.weight << std::endl;
            std::cout << edgeH.weight << std::endl;
            //if(yTest[i][j]==1||yTest[i + 1][j] ==1) {
                std::cout << eG.m_source << " " << eG.m_target << " mapped to " << eH.m_source << " " << eH.m_target
                          << " "
                          << yTest[i][j]
                          << " " << yTest[i + 1][j] << std::endl;
                //std::cout<<"i, j: "<< i<< ", "<<j<<std::endl;

                /*std::cout << "(" << G[eG.m_source].p.a << ", " << G[eG.m_source].p.b << "), (" << G[eG.m_target].p.a
                          << ", "
                          << G[eG.m_target].p.b << ") mapped to " << "(" << H[eH.m_source].p.a << ", "
                          << H[eH.m_source].p.b
                          << "), (" << H[eH.m_target].p.a << ", " << H[eH.m_target].p.b << ")"
                          << std::endl;*/
            //}

            auto angeldiff1 = anglediff(G, H, {eG.m_source, eG.m_target}, eH, yTest, x, yDecodeH, yDecodeG),
                    angeldiff2 = anglediff(G, H, {eG.m_target, eG.m_source}, eH, yTest, x, yDecodeH, yDecodeG);

            costsum += (angeldiff1 + angeldiff2) * cGlobal * yTest[i][j];
            //std::cout << costsum.getValue() << std::endl;
            costsum += (angeldiff1 + angeldiff2) * cGlobal * yTest[i + 1][j];
            std::cout << costsum.getValue() << std::endl;

            std::cout << " -------------------------" << std::endl;
            j++;
        }
        i += 2;
    }
}

void distance(Graph_boost &G, Graph_boost &H, double cGlobal, double cLocal, double cInsDel) {
    try {
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "gurobi.log");
        env.start();

        GRBModel model = GRBModel(env);

        int nh = num_vertices(H);
        int ng = num_vertices(G);
        int eh = num_edges(H);
        int eg = num_edges(G);

        std::vector<std::vector<GRBVar>> x(ng, std::vector<GRBVar>(nh));
        for (int i = 0; i < ng; i++) {
            for (int j = 0; j < nh; j++) {
                //std::string varName = "x_" + to_string(i) + "_" + to_string(j);
                x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        std::vector<std::vector<GRBVar>> y(eg * 2, std::vector<GRBVar>(eh));
        std::vector<std::vector<int>> yDecodeG(ng, std::vector<int>(ng));
        std::vector<boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long>> yDecodeH(eh);
        int i = 0, j = 0;
        for (auto edgeG: boost::make_iterator_range(edges(G))) {
            //std::cout << "testpoint 1" << std::endl;
            yDecodeG[edgeG.m_source][edgeG.m_target] = i;
            yDecodeG[edgeG.m_target][edgeG.m_source] = i + 1;
            j = 0;
            for (auto edgeH: boost::make_iterator_range(edges(H))) {
                //std::cout << "testpoint 1.2" << std::endl;
                if (i == 0) {
                    yDecodeH[j] = edgeH;
                }
                //std::string varName = "x_" + to_string(i) + "_" + to_string(j);
                y[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                y[i + 1][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                //std::cout << "testpoint 1.3" << std::endl;
                j++;
            }
            i += 2;
        }


        //std::cout << "testpoint 2" << std::endl;
        std::vector<GRBVar> theta(eg + eh);
        for (i = 0; i < eg + eh; i++) {
            theta[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }

        /*
         *
         * Constraints:
         *
         */

        // 1. Insert OR Map each edge...
        // ... in G
        int index = 0;
        i = 0, j = 0;
        for (const auto &e: boost::make_iterator_range(edges(G))) {
            std::cout << e.m_source << "    " << e.m_target << std::endl;
            model.addConstr(vecSum(y[index]) + vecSum(y[index + 1]) + theta[index / 2] == 1);
            index += 2;
        }
        // ... in H
        index = eg;
        for (const auto &e: boost::make_iterator_range(edges(H))) {
            model.addConstr(colSum(y, index - eg) + theta[index] == 1);
            index++;
        }


        // 2. Map Vertices at most once
        // ... in G
        index = 0;
        for (const auto &v: boost::make_iterator_range(vertices(G))) {
            model.addConstr(vecSum(x[index]) <= 1);
            index++;
        }

        // ... in H

        for (const auto &v: boost::make_iterator_range(vertices(H))) {
            model.addConstr(colSum(x, index - ng) <= 1);
            index++;
        }
        //std::cout << "testpoint 4" << std::endl;
        // 3. Ensure, that node and edge matching are compatible
        unsigned long sourceG, targetG, sourceH, targetH;
        for (const auto &eG: boost::make_iterator_range(edges(G))) {
            sourceG = eG.m_source;
            targetG = eG.m_target;
            j = 0;
            for (const auto &eH: boost::make_iterator_range(edges(H))) {
                sourceH = eH.m_source;
                targetH = eH.m_target;
                model.addQConstr(y[i][j] <= x[sourceG][sourceH]);
                model.addConstr(y[i + 1][j] <= x[sourceG][targetH]);
                model.addQConstr(y[i][j] <= x[targetG][targetH]);
                model.addQConstr(y[i + 1][j] <= x[targetG][sourceH]);
                j++;
            }
            i += 2;
        }

        std::cout << "testpoint 5" << std::endl;
        // Set cost function as objective
        model.setObjective(insdel(theta, G, H) * cInsDel + mapcost(y, x, G, H, cGlobal, cLocal, yDecodeG, yDecodeH),
                           GRB_MINIMIZE);
        std::cout << "testpoint 6" << std::endl;
        // Optimieren
        model.optimize();
        printResult(model);
        i = 0, index = 0;
        for (auto eG: boost::make_iterator_range(edges(G))) {
            auto &edgeG = G[eG];
            j = 0;
            for (auto eH: boost::make_iterator_range(edges(H))) {
                auto &edgeH = H[eH];
                if (y[i][j].get(GRB_DoubleAttr_X) == 1) {
                    std::cout << eG.m_source << " " << eG.m_target << " mapped to " << eH.m_source << " " << eH.m_target
                              << ": " << y[i][j].get(GRB_DoubleAttr_X) << std::endl;
                    std::cout << "(" << G[eG.m_source].p.a << ", " << G[eG.m_source].p.b << "), (" << G[eG.m_target].p.a
                              << ", "
                              << G[eG.m_target].p.b << ") mapped to " << "(" << H[eH.m_source].p.a << ", "
                              << H[eH.m_source].p.b
                              << "), (" << H[eH.m_target].p.a << ", " << H[eH.m_target].p.b << ")"
                              << std::endl;
                }
                if (y[i + 1][j].get(GRB_DoubleAttr_X) == 1) {
                    std::cout << eG.m_target << " " << eG.m_source << " mapped to " << eH.m_source << " " << eH.m_target
                              << ": " << y[i + 1][j].get(GRB_DoubleAttr_X) << std::endl;
                    std::cout << "(" << G[eG.m_source].p.a << ", " << G[eG.m_source].p.b << "), (" << G[eG.m_target].p.a
                              << ", "
                              << G[eG.m_target].p.b << ") mapped to " << "(" << H[eH.m_source].p.a << ", "
                              << H[eH.m_source].p.b
                              << "), (" << H[eH.m_target].p.a << ", " << H[eH.m_target].p.b << ")"
                              << std::endl;
                }

                if (y[i][j].get(GRB_DoubleAttr_X) == 1) {
                    edgeG.label = index;
                    edgeH.label = index;
                    index++;
                }
                if (y[i + 1][j].get(GRB_DoubleAttr_X) == 1) {
                    edgeG.label = index;
                    edgeH.label = index;
                    index++;
                }
                //std::cout << edgeG.label << "  " << edgeH.label << std::endl;
                //costsum += anglediff(G, H, eG, eH, y, x, yDecodeH, yDecodeG) * cGlobal;
                j++;
            }
            i += 2;
        }

        std::cout << "-------testest-------" << std::endl;
        std::cout<<"insdel: "<<insdel(theta, G, H).getValue()<<std::endl;
        std::cout<<"mapcost: "<<mapcost(y, x, G, H, cGlobal, cLocal, yDecodeG, yDecodeH).getValue()<<std::endl;

        //Testtest(G, H, y, x, cLocal, cGlobal, yDecodeH, yDecodeG);

    } catch (GRBException &e) {
        std::cerr << "Gurobi Error: " << e.getMessage() << std::endl;
    }
}

double distanceVal(const Graph_boost &G, const Graph_boost &H, double cGlobal, double cLocal, double cInsDel) {
    try {
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "gurobi.log");
        env.set("OutputFlag", "0");
        env.start();

        GRBModel model = GRBModel(env);

        int nh = num_vertices(H);
        int ng = num_vertices(G);
        int eh = num_edges(H);
        int eg = num_edges(G);

        std::vector<std::vector<GRBVar>> x(ng, std::vector<GRBVar>(nh));
        for (int i = 0; i < ng; i++) {
            for (int j = 0; j < nh; j++) {
                //std::string varName = "x_" + to_string(i) + "_" + to_string(j);
                x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        std::vector<std::vector<GRBVar>> y(eg * 2, std::vector<GRBVar>(eh));
        std::vector<std::vector<int>> yDecodeG(ng, std::vector<int>(ng));
        std::vector<boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long>> yDecodeH(eh);
        int i = 0, j = 0;
        for (auto edgeG: boost::make_iterator_range(edges(G))) {
            //std::cout << "testpoint 1" << std::endl;
            yDecodeG[edgeG.m_source][edgeG.m_target] = i;
            yDecodeG[edgeG.m_target][edgeG.m_source] = i + 1;
            j = 0;
            for (auto edgeH: boost::make_iterator_range(edges(H))) {
                //std::cout << "testpoint 1.2" << std::endl;
                if (i == 0) {
                    yDecodeH[j] = edgeH;
                }
                //std::string varName = "x_" + to_string(i) + "_" + to_string(j);
                y[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                y[i + 1][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                //std::cout << "testpoint 1.3" << std::endl;
                j++;
            }
            i += 2;
        }


        //std::cout << "testpoint 2" << std::endl;
        std::vector<GRBVar> theta(eg + eh);
        for (i = 0; i < eg + eh; i++) {
            theta[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        }

        /*
         *
         * Constraints:
         *
         */

        // 1. Insert OR Map each edge...
        // ... in G
        int index = 0;
        i = 0, j = 0;
        for (const auto &e: boost::make_iterator_range(edges(G))) {
            //std::cout << e.m_source << "    " << e.m_target << std::endl;
            model.addConstr(vecSum(y[index]) + vecSum(y[index + 1]) + theta[index / 2] == 1);
            index += 2;
        }
        // ... in H
        index = eg;
        for (const auto &e: boost::make_iterator_range(edges(H))) {
            model.addConstr(colSum(y, index - eg) + theta[index] == 1);
            index++;
        }


        // 2. Map Vertices at most once
        // ... in G
        index = 0;
        for (const auto &v: boost::make_iterator_range(vertices(G))) {
            model.addConstr(vecSum(x[index]) <= 1);
            index++;
        }

        // ... in H

        for (const auto &v: boost::make_iterator_range(vertices(H))) {
            model.addConstr(colSum(x, index - ng) <= 1);
            index++;
        }
        //std::cout << "testpoint 4" << std::endl;
        // 3. Ensure, that node and edge matching are compatible
        unsigned long sourceG, targetG, sourceH, targetH;
        for (const auto &eG: boost::make_iterator_range(edges(G))) {
            sourceG = eG.m_source;
            targetG = eG.m_target;
            j = 0;
            for (const auto &eH: boost::make_iterator_range(edges(H))) {
                sourceH = eH.m_source;
                targetH = eH.m_target;
                model.addQConstr(y[i][j] <= x[sourceG][sourceH]);
                model.addConstr(y[i + 1][j] <= x[sourceG][targetH]);
                model.addQConstr(y[i][j] <= x[targetG][targetH]);
                model.addQConstr(y[i + 1][j] <= x[targetG][sourceH]);
                j++;
            }
            i += 2;
        }

        //std::cout << "testpoint 5" << std::endl;
        // Set cost function as objective
        model.setObjective(insdel(theta, G, H) * cInsDel + mapcost(y, x, G, H, cGlobal, cLocal, yDecodeG, yDecodeH),
                           GRB_MINIMIZE);
        //std::cout << "testpoint 6" << std::endl;
        // Optimieren
        model.optimize();
        //printResult(model);
        /*i = 0, index = 0;
        for (auto eG: boost::make_iterator_range(edges(G))) {
            auto &edgeG = G[eG];
            j = 0;
            for (auto eH: boost::make_iterator_range(edges(H))) {
                auto &edgeH = H[eH];
                if (y[i][j].get(GRB_DoubleAttr_X) == 1) {
                    std::cout << eG.m_source << " " << eG.m_target << " mapped to " << eH.m_source << " " << eH.m_target
                              << ": " << y[i][j].get(GRB_DoubleAttr_X) << std::endl;
                    std::cout << "(" << G[eG.m_source].p.a << ", " << G[eG.m_source].p.b << "), (" << G[eG.m_target].p.a
                              << ", "
                              << G[eG.m_target].p.b << ") mapped to " << "(" << H[eH.m_source].p.a << ", "
                              << H[eH.m_source].p.b
                              << "), (" << H[eH.m_target].p.a << ", " << H[eH.m_target].p.b << ")"
                              << std::endl;
                }
                if (y[i + 1][j].get(GRB_DoubleAttr_X) == 1) {
                    std::cout << eG.m_target << " " << eG.m_source << " mapped to " << eH.m_source << " " << eH.m_target
                              << ": " << y[i + 1][j].get(GRB_DoubleAttr_X) << std::endl;
                    std::cout << "(" << G[eG.m_source].p.a << ", " << G[eG.m_source].p.b << "), (" << G[eG.m_target].p.a
                              << ", "
                              << G[eG.m_target].p.b << ") mapped to " << "(" << H[eH.m_source].p.a << ", "
                              << H[eH.m_source].p.b
                              << "), (" << H[eH.m_target].p.a << ", " << H[eH.m_target].p.b << ")"
                              << std::endl;
                }

                if (y[i][j].get(GRB_DoubleAttr_X) == 1) {
                    edgeG.label = index;
                    edgeH.label = index;
                    index++;
                }
                if (y[i + 1][j].get(GRB_DoubleAttr_X) == 1) {
                    edgeG.label = index;
                    edgeH.label = index;
                    index++;
                }
                //std::cout << edgeG.label << "  " << edgeH.label << std::endl;
                //costsum += anglediff(G, H, eG, eH, y, x, yDecodeH, yDecodeG) * cGlobal;
                j++;
            }
            i += 2;
        }

        //std::cout << "-------testest-------" << std::endl;
        std::cout<<"insdel: "<<insdel(theta, G, H).getValue()<<std::endl;
        std::cout<<"mapcost: "<<mapcost(y, x, G, H, cGlobal, cLocal, yDecodeG, yDecodeH).getValue()<<std::endl;*/

        //Testtest(G, H, y, x, cLocal, cGlobal, yDecodeH, yDecodeG);
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) return model.get(GRB_DoubleAttr_ObjVal);
        else return -1;

    } catch (GRBException &e) {
        std::cerr << "Gurobi Error: " << e.getMessage() << std::endl;
        return -1;
    }
}