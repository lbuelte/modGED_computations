#include "../include/hierachicalClustering3.h"
#include <chrono>
#include <fstream>


struct NearestSecondPair {
    int nearest;
    int second;
    double nearestDist;
    double secondDist;
};

double distTimePam, pamTime, distTimeCLARA;

std::vector<std::vector<double>>
computeDistanceMatrix(const std::vector<Graph_boost> &graphs, const std::vector<int> indexSet,
                      const std::vector<double> &edgeSums,
                      const std::vector<double> &longest) {
    int n = indexSet.size();
    std::vector<std::vector<double>> D(n, std::vector<double>(n, 0.0));
    //std::cout << "distance: ";
    for (int i = 0; i < n; i++) {
        int indexi = indexSet[i];
        for (int j = 0; j < i; j++) {
            int indexj = indexSet[j];
            if (i == j) continue;
            double dist = distVal_est(graphs[indexi], graphs[indexj], edgeSums[indexi], edgeSums[indexj],
                                      longest[indexi], longest[indexj]);
            D[i][j] = dist;
            D[j][i] = dist;
            //std::cout << ", "<<dist;
        }
    }
    //std::cout <<  std::endl;
    return D;
}


void saveClusterAssignments(
        const std::vector<int> &labels,
        const std::string &filename = "cluster_assignments.txt"
) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        throw std::runtime_error("Konnte Datei nicht öffnen: " + filename);
    }

    for (int i = 0; i < labels.size(); ++i) {
        out << i << " " << labels[i] << "\n";
    }

    out.close();
}

std::vector<NearestSecondPair>
nearest_second_pairs(std::vector<std::vector<double>> distMatrix, std::vector<int> medoids) {
    std::vector<NearestSecondPair> nearests;

    for (int i = 0; i < distMatrix.size(); i++) {
        auto nearestDist = DBL_MAX;
        auto secondDist = DBL_MAX;
        int nearest = -1, second = -1;

        for (int m = 0; m < medoids.size(); m++) {
            double d = distMatrix[i][medoids[m]];
            if (d < nearestDist) {
                second = nearest;
                secondDist = nearestDist;
                nearestDist = d;
                nearest = m;
            } else if (d < secondDist) {
                secondDist = d;
                second = m;
            }
        }

        nearests.push_back({nearest, second, nearestDist, secondDist});
    }

    return nearests;
}


std::vector<int> sampleIndices(const std::vector<int> &indices, int s) {
    std::vector<int> result;
    std::sample(indices.begin(), indices.end(), std::back_inserter(result),
                s, std::mt19937{std::random_device{}()});
    return result;
}


std::vector<int> sampleIndices(int n, int s) {
    std::vector<int> result;
    auto range = std::vector<int>(n);
    std::iota(range.begin(), range.end(), 0);
    std::sample(range.begin(), range.end(), std::back_inserter(result),
                s, std::mt19937{std::random_device{}()});
    return result;
}

std::vector<int> PAMBuild(const std::vector<std::vector<double>> &D, int k) {
    int n = D.size();
    auto TD = DBL_MAX;
    int m1 = -1;

    for (int xc = 0; xc < n; xc++) {
        double TDj = 0.0;
        for (int xo = 0; xo < n; xo++) {
            if (xo != xc) TDj += D[xo][xc];
        }
        if (TDj < TD) {
            TD = TDj;
            m1 = xc;
        }
    }

    std::vector<int> medoids = {m1};
    std::vector<double> dnearest(n, DBL_MAX);

    for (int xo = 0; xo < n; ++xo) {
        if (xo != m1) {
            dnearest[xo] = D[m1][xo];
        }
    }

    for (int i = 1; i < k; ++i) {
        double deltaTDBest = DBL_MAX;
        int xBest = 0;

        for (int xc = 0; xc < n; xc++) {
            if (find(medoids.begin(), medoids.end(), xc) != medoids.end()) continue;

            double deltaTD = 0.0;

            for (int xo = 0; xo < n; xo++) {
                if (find(medoids.begin(), medoids.end(), xo) != medoids.end() || xo == xc) continue;

                double delta = D[xo][xc] - dnearest[xo];
                if (delta < 0) deltaTD += delta;
            }

            if (deltaTD < deltaTDBest) {
                deltaTDBest = deltaTD;
                xBest = xc;
            }
        }

        medoids.push_back(xBest);
        TD += deltaTDBest;


        for (int xo = 0; xo < n; ++xo) {
            if (find(medoids.begin(), medoids.end(), xo) == medoids.end()) {
                dnearest[xo] = std::min(dnearest[xo], D[medoids.back()][xo]);
            }
        }
    }

    return medoids;
}

void print_vector(std::vector<int> v) {
    for (auto i: v) {
        std::cout << i << ", ";
    }
    std::cout << std::endl;
}


std::vector<int> fast_pam1(const std::vector<Graph_boost> &graphs,
                           const std::vector<int> &dataIndices,
                           int k, std::vector<double> EdgeSums, std::vector<double> longest, bool randomized) {


    std::vector<int> medoids;

    //compute distance matrix
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double>> distMatrix = computeDistanceMatrix(graphs, dataIndices, EdgeSums, longest);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    distTimePam += duration.count();
    // compute intial medoids
    if (randomized) {
        medoids = sampleIndices(dataIndices.size(), k);
    } else {
        medoids = PAMBuild(distMatrix, k);
    }

    std::cout << "medoids: " << medoids.size() << " " << medoids[0] << std::endl;


    auto nearests = nearest_second_pairs(distMatrix, medoids);
    std::cout << "nearests: " << nearests.size() << std::endl;
    int kill = 0;
    while (true) {
        kill++;
        //if(kill> 50) break;
        std::vector<double> deltaTD_minus(k, 0.0);
        double TDbest = 0;
        int mbest = -1, xbest = -1;

        //compute TD^-m_i
        for (int xo = 0; xo < dataIndices.size(); xo++) {
            auto TD_xo = nearests[xo];
            deltaTD_minus[TD_xo.nearest] += TD_xo.secondDist - TD_xo.nearestDist;
        }

        //std::cout << "TD^-m_i: " << deltaTD_minus[0] << " " << deltaTD_minus[1] << std::endl;

        for (int xc = 0; xc < dataIndices.size(); xc++) {
            //only iterate over xc that are not already medoids
            if (find(medoids.begin(), medoids.end(), xc) != medoids.end())
                continue;
            double TDplusxc = 0;

            auto deltaTD = deltaTD_minus;
            for (int xo = 0; xo < dataIndices.size(); xo++) {
                //if (xo == xc) continue;
                double d_oj = distMatrix[xo][xc];

                const auto &localTD = nearests[xo];

                if (d_oj < localTD.nearestDist) {
                    TDplusxc += d_oj - localTD.nearestDist;
                    deltaTD[localTD.nearest] += (localTD.nearestDist - localTD.secondDist);
                    //std::cout << "nearest"<<localTD.nearestDist-localTD.secondDist<<std::endl;
                } else if (d_oj < localTD.secondDist) {
                    deltaTD[localTD.nearest] += (d_oj - localTD.secondDist);
                    //std::cout << "second"<<d_oj-localTD.secondDist<<std::endl;;
                }

            }
            auto minTD = std::min_element(deltaTD.begin(), deltaTD.end());
            int i = std::distance(deltaTD.begin(), minTD);

            //std::cout << "TD^-m_i: " << deltaTD[i] << "   " << deltaTD_minus[i] << "    " << TDplusxc << std::endl;
            deltaTD[i] += TDplusxc;
            if (deltaTD[i] < TDbest) {
                TDbest = deltaTD[i];
                mbest = i;
                xbest = xc;
            }
            //std::cout << "TDbest: " << TDbest << " " << mbest << " " << xbest << std::endl;

        }
        std::cout << "kill: " << kill << "   " << TDbest << std::endl;
        if (TDbest >= -0.00000000001) break;
        medoids[mbest] = xbest;
        nearests = nearest_second_pairs(distMatrix, medoids);
    }

    std::vector<int> finalMedoids;
    for (int m: medoids) {
        finalMedoids.push_back(dataIndices[m]);
    }
    std::cout << "finalMedoids: ";
    print_vector(finalMedoids);
    return finalMedoids;

}

std::pair<std::vector<int>, std::vector<int>> claraClustering(const std::vector<Graph_boost> &graphs,
                                                              const std::vector<int> &options,
                                                              std::vector<double> EdgeSums,
                                                              std::vector<double> longest,
                                                              int k, int sampleSize, int numSamples) {
    std::pair<std::vector<int>, std::vector<int>> bestResult;
    double bestCost = std::numeric_limits<double>::max();


    for (int s = 0; s < numSamples; ++s) {
        // random sample
        std::vector<int> sample = sampleIndices(options, sampleSize);
        std::cout << "sample: " << sample.size() << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<int> medoids = fast_pam1(graphs, sample, k, EdgeSums, longest, false);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        pamTime += duration.count();
        //std::vector<int> medoids = pamParkJunSampled(graphs, k, EdgeSums, longest, sample, 100).first;
        std::cout << medoids[0] << " medoid " << medoids[1] << std::endl;

        //calculate total costs
        // 3. Jedem Graph den nächsten Medoid zuweisen (auf kompletter Menge)
        std::cout << "sample: " << sample.size() << std::endl;
        std::vector<int> labels(options.size(), 0);
        double total_cost = 0.0;
        bool notWorse = true;
        for (int i = 0; i < options.size(); ++i) {
            double min_dist = std::numeric_limits<double>::max();
            int best_medoid = 0;
            for (int j = 0; j < k; ++j) {
                if (abs(EdgeSums[options[i]] - EdgeSums[medoids[j]]) > min_dist) continue;
                start = std::chrono::high_resolution_clock::now();
                double dist = distVal_est(graphs[options[i]], graphs[medoids[j]], EdgeSums[options[i]],
                                          EdgeSums[medoids[j]], longest[i],
                                          longest[medoids[j]]);
                end = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
                distTimeCLARA += duration.count();
                if (dist < min_dist) {
                    min_dist = dist;
                    best_medoid = j;
                }
            }
            labels[i] = best_medoid;
            total_cost += min_dist;
            if (total_cost >= bestCost) {
                notWorse = false;
                break;
            }
        }



        // 4. Save best
        if (notWorse && total_cost < bestCost) {
            bestCost = total_cost;
            bestResult = {medoids, labels};
        }
    }
    std::cout << "bestResult: ";
    print_vector(bestResult.first);
    return bestResult;
}


void refineCluster(std::vector<int> &clusterassignment,
                   const std::vector<int> &refinement, const std::vector<int> &options, int hierarchyFactor) {
    for (auto i = 0; i < options.size(); i++) {
        int graphIndex = options[i];
        clusterassignment[graphIndex] += (refinement[i] + 1) * hierarchyFactor;
        std::cout << "clusterassignment[i] = " << clusterassignment[graphIndex] << std::endl;
        std::cout << "ref[i] = " << refinement[i] << std::endl;
    }

}

void refineCluster2(std::vector<int> &clusterassignment,
                    const std::vector<int> &refinement, const std::vector<int> &options, int hierarchyFactor) {
    for (auto i: options) {
        clusterassignment[i] += (refinement[i] + 1) * hierarchyFactor;
        //std::cout << "clusterassignment[i] = " << clusterassignment[i] << std::endl;
        //std::cout << "ref[i] = " << refinement[i] << std::endl;
    }

}

void
hierarchicalClustering(const std::vector<Graph_boost> &graphs, std::vector<int> &clusterassignment, int h, int k,
                       double epsilon) {
    std::vector<int> options;
    std::vector<double> EdgeSums, longest;
    double edgesum = 0, longestedge = 0;
    std::vector<Graph_boost> graphsCurr;
    int clusterID = 0;
    int hierarchyFactor = std::pow(10, h);
    for (int hCurr = 0; hCurr < h; hCurr++) {
        //lower hierarchyfactor

        hierarchyFactor /= 10;

        // Filter and smooth graphs
        std::cout << "hCurr: " << hCurr << std::endl;
        for (auto g: graphs) {
            graphsCurr.push_back(
                    simplifyAndFilterGraphPathsWithEdgeSum(g, epsilon, edgesum, longestedge, 0.87 - (hCurr * 0.175)));
            EdgeSums.push_back(edgesum);
            longest.push_back(longestedge);
            edgesum = 0;
            longestedge = 0;
        }
        std::cout << "Vereinfacht, alter" << std::endl;
        // Refine each cluster if size at least k

        for (int i = 1; i < std::pow(k, hCurr) + 1; i++) {
            std::cout << "ID: " << clusterID << std::endl;
            for (int index = 0; index < clusterassignment.size(); index++) {
                if (clusterassignment[index] == clusterID) {
                    options.push_back(index);
                }
            }
            if (options.size() < k) {
                continue;
            }
            std::cout << graphs.size() << "   " << options.size() << std::endl;
            std::vector<int> result;
            if (options.size() >= 3 * k)
                result = claraClustering(graphsCurr, options, EdgeSums, longest, k, 3 * k, 2).second;
                //if(options.size() >= 40+2*k) result = claraClustering(graphsCurr, options, EdgeSums, longest, k, 40+2*k, 5).second;
                //else result = claraClustering(graphsCurr, options, EdgeSums, longest, k, options.size(), 5).second;
            else result = claraClustering(graphsCurr, options, EdgeSums, longest, k, options.size(), 2).second;
            std::cout << "geclustert, alter" << std::endl;
            refineCluster(clusterassignment, result, options, hierarchyFactor);
            clusterID += i * (static_cast<int>(std::pow(10, h - hCurr - 1)));
            options.clear();
        }
        graphsCurr.clear();
        EdgeSums.clear();
        longest.clear();
    }
}

void generateHierarchical(int h, int k, std::vector<int> &result) {
    // Tiefe bis h-1, damit letzte Stelle 0 bleibt
    for (int depth = 1; depth < h; ++depth) {
        print_vector(result);
        int total = std::pow(k, depth);
        for (int i = 0; i < total; ++i) {
            int id = 0;
            int tmp = i;
            for (int d = 0; d < depth; ++d) {
                int digit = tmp % k + 1;
                id += digit * std::pow(10, h - d - 1);
                tmp /= k;
            }
            // Letzte Ziffer ist immer 0 → nichts hinzufügen
            result.push_back(id);
        }
        print_vector(result);
        result.clear();

    }
}

std::vector<int> generateHierarchicalIDs(int h, int k) {
    std::vector<int> ids;
    ids.push_back(0);
    generateHierarchical(h, k, ids);
    return ids;
}


void
hierarchicalClustering2(const std::vector<Graph_boost> &graphs, std::vector<int> &clusterassignment, int h, int k,
                        double epsilon, double initCutOff) {

    std::vector<int> ids;
    ids.push_back(0);
    std::vector<Graph_boost> graphsCurr;
    std::vector<int> options;
    std::vector<double> EdgeSums, longest;
    double edgesum = 0, longestedge = 0;
    for (int depth = 0; depth < h; ++depth) {
        for (auto g: graphs) {
            graphsCurr.push_back(
                    simplifyAndFilterGraphPathsWithEdgeSum(g, epsilon, edgesum, longestedge,
                                                           initCutOff - ((depth) * 0.175)));
            EdgeSums.push_back(edgesum);
            longest.push_back(longestedge);
            edgesum = 0;
            longestedge = 0;
        }
        int total = std::pow(k, depth);
        for (int i = 0; i < total; ++i) {
            int id = 0;
            int tmp = i;
            for (int d = 0; d < depth; ++d) {
                int digit = tmp % k + 1;
                id += digit * std::pow(10, h - d - 1);
                tmp /= k;
            }

            ids.push_back(id);
        }

        for (auto id: ids) {

            std::cout << "ID: " << id << std::endl;
            for (int index = 0; index < clusterassignment.size(); index++) {
                if (clusterassignment[index] == id) {
                    options.push_back(index);
                }
            }
            if (options.size() < k) {
                continue;
            }
            std::cout << graphs.size() << "   " << options.size() << std::endl;
            std::vector<int> result;
            //if(options.size() >= 3*k) result = claraClustering(graphsCurr, options, EdgeSums, longest, k, 3*k, 2).second;
            if (options.size() >= 40 + 2 * k)
                result = claraClustering(graphsCurr, options, EdgeSums, longest, k, 40 + 2 * k, 5).second;
            else result = claraClustering(graphsCurr, options, EdgeSums, longest, k, options.size(), 5).second;
            //else result = claraClustering(graphsCurr, options, EdgeSums, longest, k, options.size(), 2).second;
            std::cout << "geclustert, alter" << std::endl;
            refineCluster(clusterassignment, result, options, std::pow(10, h - depth - 1));
            options.clear();
        }


        graphsCurr.clear();
        EdgeSums.clear();
        longest.clear();
        ids.clear();

    }
    //print_vector(clusterassignment);
}


std::vector<int>
get_hierarchicalClustering(const std::vector<Graph_boost> &input_graphs, int h, int k, double epsilon,
                           double initCutOff, std::string filename) {
    distTimePam = 0;
    distTimeCLARA=0;
    pamTime = 0;
    std::vector<int> result = std::vector<int>(input_graphs.size(), 0);
    hierarchicalClustering2(input_graphs, result, h, k, epsilon, initCutOff);
    std::ofstream datei("../output/timelogh2k2zentrum.txt", std::ios::app);
    if (!datei) {
        std::cerr << "Fehler beim Öffnen der Datei!\n";
    }
    datei << initCutOff << "," << distTimeCLARA<<","<< distTimePam<<","<< pamTime<<std::endl;
    datei.close();
    //print_vector(result);
    saveClusterAssignments(result, filename);
    return result;
}