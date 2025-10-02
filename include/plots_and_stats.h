#ifndef PLOTS_AND_STATS_H
#define PLOTS_AND_STATS_H

#include "graph_boost.h"

namespace{

void write_svg(Graph_boost const &graph, std::string const &filename){

    std::cout << "Plotting graph with vertices " << std::endl;
    for (auto const &v : graph.m_vertices){
        std::cout << v.m_property.p.a << " " << v.m_property.p.b << std::endl;
    }
    std::cout << " and edges " << std::endl;
    for (auto const &e : graph.m_edges){
        std::cout << e.m_source << e.m_target << " - " << e.m_property.weight << std::endl;
    }

    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Error: cannot open " << filename << "\n";
        return;
    }

    ofs << "<svg xmlns='http://www.w3.org/2000/svg' "
        << "width='1000' height='1000' viewBox='0 0 200 200'>\n";

    ofs << "<rect width=\"100%\" height=\"100%\" fill=\"lightgrey\" />\n";

    long double min_x = std::numeric_limits<long double>::max();
    long double max_x = std::numeric_limits<long double>::lowest();
    long double min_y = std::numeric_limits<long double>::max();
    long double max_y = std::numeric_limits<long double>::lowest();

    for (auto& v : graph.m_vertices) {
        min_x = std::min(min_x, v.m_property.p.a);
        max_x = std::max(max_x, v.m_property.p.a);
        min_y = std::min(min_y, v.m_property.p.b);
        max_y = std::max(max_y, v.m_property.p.b);
    }

    double width  = max_x - min_x;
    double height = max_y - min_y;
    double scale  = 180.0 / std::max(width, height);

    for (auto const &edge : graph.m_edges){
        auto const v = graph.m_vertices[edge.m_source];
        auto const w = graph.m_vertices[edge.m_target];

        double x1 = (v.m_property.p.a - min_x) * scale + 10;
        double y1 = (v.m_property.p.b - min_y) * scale + 10;
        double x2 = (w.m_property.p.a - min_x) * scale + 10;
        double y2 = (w.m_property.p.b - min_y) * scale + 10;

        ofs << "<line x1='" << x1 << "' y1='" << y1
           << "' x2='" << x2 << "' y2='" << y2
           << "' stroke='black' stroke-width='1'/>\n";

        double mx = ((v.m_property.p.a + w.m_property.p.a) / 2.0 - min_x) * scale + 20;
        double my = ((v.m_property.p.b + w.m_property.p.b) / 2.0 - min_y) * scale + 10;

        ofs << "<text x='" << mx
            << "' y='" << my
            << "' font-size='4' text-anchor='middle' "
               "dominant-baseline='middle' fill='blue'>"
            << edge.m_property.weight
            << "</text>\n";
    }

    int i = 0;
    for (auto const &v : graph.m_vertices) {
        double x = (v.m_property.p.a - min_x) * scale + 10;
        double y = (v.m_property.p.b - min_y) * scale + 10;
        ofs << "<circle cx='" << x << "' cy='" << y
            << "' r='3' fill='red'/>\n";
        ofs << "<text x='" << x << "' y='" << y 
            << "' font-size='4' text-anchor='middle' dominant-baseline='middle' "
               "fill='black' stroke='white' stroke-width='0.8' paint-order='stroke'>"
            << i << "</text>\n";
        i++;
    }
    
    ofs << "<text x='" << 195
        << "' y='" << 195
        << "' font-size='4' text-anchor='end' "
        "dominant-baseline='auto' fill='black'>"
        << filename
        << "</text>\n";

    ofs << "</svg>\n";
}

}

void plot_graphs_as_svg(std::vector<Graph_boost> const &graphs, std::string const &foldername, size_t counter_start){
    std::cout << "Plotting " << graphs.size() << " graphs ..." << std::endl;
    int i = counter_start;
    for (auto const &graph : graphs){
        std::string filename = foldername + "/graph_" + std::to_string(i) + ".svg";
        write_svg(graph, filename);
        i++;
    }
}


struct Stats{
    size_t num_runs = 0;
    size_t num_timeouts_f2plus = 0;
    size_t num_timeouts_fori = 0;
    size_t num_timeouts_both = 0;

    size_t num_fori_smaller = 0;
    size_t num_f2plus_smaller = 0;

    size_t duration_count_f2plus = 0;
    size_t duration_count_fori = 0;

    std::vector<std::pair<int,int>> TO_instances_f2plus;
    std::vector<std::pair<int,int>> TO_instances_fori;

    void merge(const Stats &other){
        num_runs += other.num_runs;
        num_timeouts_f2plus += other.num_timeouts_f2plus;
        num_timeouts_fori += other.num_timeouts_fori;
        num_timeouts_both += other.num_timeouts_both;
        num_fori_smaller += other.num_fori_smaller;
        num_f2plus_smaller += other.num_f2plus_smaller;
        duration_count_f2plus += other.duration_count_f2plus;
        duration_count_fori += other.duration_count_fori;
        TO_instances_f2plus.insert(TO_instances_f2plus.end(), other.TO_instances_f2plus.begin(), other.TO_instances_f2plus.end());
        TO_instances_fori.insert(TO_instances_fori.end(), other.TO_instances_fori.begin(), other.TO_instances_fori.end());        
    }

    void print_stats();
    void print_fori_stats();
    void print_TO_instances();
    void print_fori_TO_instances();
    void print_f2plus_TO_instances();
};

void Stats::print_stats(){
    std::cout << std::endl << "Printing statistics:" << std::endl
              << "In total " << num_runs << " runs, there were " << num_timeouts_f2plus
              << " timeout of F2Plus, " << num_timeouts_fori << " timeout of FORI and "
              << num_timeouts_both << " timeouts of both." << std::endl
              << "The total duration for F2Plus was " << duration_count_f2plus 
              << " milliseconds, the total duration for FORI was " << duration_count_fori 
              << " milliseconds." << std::endl
              << "The solution of FORI was smaller " << num_fori_smaller 
              << " times, the solution of F2Plus was smaller " << num_f2plus_smaller << " times."
              << std::endl;
}

void Stats::print_fori_stats(){
    std::cout << std::endl << "Printing FORI statistics:" << std::endl
              << "In total " << num_runs << " runs, there were " << num_timeouts_fori 
              << " timeout of FORI." << std::endl
              << "The total duration for FORI was " << duration_count_fori 
              << " milliseconds, which amount to " 
              << static_cast<double>(duration_count_fori)/static_cast<double>(num_runs) 
              << " milliseconds per run." << std::endl;
}

void Stats::print_TO_instances(){
    print_fori_TO_instances();
    print_f2plus_TO_instances();
}

void Stats::print_fori_TO_instances(){
    std::cout << "FORI timed out on the instances: ";
    for (auto const & instance_pair : TO_instances_fori){
        std::cout << "(" << instance_pair.first << "," << instance_pair.second << ") ";
    }
    std::cout << std::endl;
}

void Stats::print_f2plus_TO_instances(){
    std::cout << "F2Plus timed out on the instances: ";
    for (auto const & instance_pair : TO_instances_f2plus){
        std::cout << "(" << instance_pair.first << "," << instance_pair.second << ") ";
    }    
}

#endif //PLOT_INSTANCES_H