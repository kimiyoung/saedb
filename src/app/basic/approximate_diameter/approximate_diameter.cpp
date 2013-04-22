#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include "sae_include.hpp"

const double FM_COEFFICIENT = 0.77351;
const double EPSILON = 1e-4;
const int MAX_ITERATION_TIMES = 100;
const int DUPLICATION_OF_BITMASKS = 10;

int hash_func(){
    int ret = 0;
    while ((double)rand() / RAND_MAX < 0.5) ret ++;
    return ret;
}

typedef std::vector<std::vector<bool> > bitmask_type;

struct vertex_data_type{
    bitmask_type bitmask;

    vertex_data_type(): bitmask() {}

    void create_hashed_bitmask(){
        for (int i=0; i<DUPLICATION_OF_BITMASKS; i++){
            int hash_val = hash_func();
            std::vector<bool> mask(hash_val + 2, 0);
            mask[hash_val] = 1;
            bitmask.push_back(mask);
        }
    }

    // TODO save and load
};

typedef saedb::sae_graph<vertex_data_type, double> graph_type;

void bitwise_or(bitmask_type& v1, const bitmask_type& v2){
    for (int i=0; i<v1.size(); i++){
        while (v1[i].size() <v2[i].size()) 
            v1[i].push_back(0);

        for (int j=0; j<v2[i].size(); j++)
            v1[i][j] = v1[i][j] || v2[i][j];
    }
}

struct gather_type {
    bitmask_type bitmask;

    gather_type(): bitmask() {}

    explicit gather_type(const bitmask_type& v): bitmask(){
        for (int i=0; i<v.size(); i++)
            bitmask.push_back(v[i]);
    }

    gather_type& operator +=(const gather_type& other){
        bitwise_or(bitmask, other.bitmask);
        return *this;
    }

    // TODO save and load
};

class InitGraph: public saedb::IAlgorithm<graph_type, saedb::empty> {
public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const{
        return saedb::NO_EDGES;
    }

    gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const{   
        return gather_type();
    }

    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total){
        vertex.data().create_hashed_bitmask();
    }

    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
        return saedb::NO_EDGES;
    }

    void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const {
    }
};

class ApproximateDiameter: public saedb::IAlgorithm<graph_type, gather_type> {
public:
    edge_dir_type gather_edges(icontext_type& context, const vertex_type& vertex) const{
        return saedb::OUT_EDGES;
    }

    gather_type gather(icontext_type& context, const vertex_type& vertex, edge_type& edge) const{
        return gather_type(edge.target().data().bitmask);
    }

    void apply(icontext_type& context, vertex_type& vertex, const gather_type& total){
        if (total.bitmask.size() > 0)
            bitwise_or(vertex.data().bitmask, total.bitmask);
    }

    edge_dir_type scatter_edges(icontext_type& context, const vertex_type& vertex) const {
        return saedb::NO_EDGES;
    }

    void scatter(icontext_type& context, const vertex_type& vertex, edge_type& edge) const{
    }
};

double approximate_pair_number(ApproximateDiameter::icontext_type& context,
    const graph_type::vertex_type& vertex){
    double sum = 0.0;
    const bitmask_type& bitmask_temp = vertex.data().bitmask;

    for (int i=0; i<bitmask_temp.size(); i++) {
        for (int j=0; j<=bitmask_temp[i].size(); j++) {
            if (j == bitmask_temp[i].size() || !bitmask_temp[i][j]) {
                sum += (double) i;
                break;
            }
        }
    }
    return pow(2.0, sum / (double) bitmask_temp.size()) / FM_COEFFICIENT;
}

/*\reference
HADI: Fast Diameter Estimation and Mining in Massive Graphs with Hadoop, U Kang, 2008
Probabilistic Counting Algorithms for Data Base Applications, P.Flajolet, G.N.Martin, 1985
*/

int main() {
    graph_type graph;
    graph.load_format("test_graph");

    saedb::IEngine<InitGraph> *initEngine = new saedb::EngineDelegate<InitGraph>(graph);
    initEngine->signalAll();
    initEngine->start();
    delete initEngine;

    saedb::IEngine<ApproximateDiameter> *engine = new saedb::EngineDelegate<ApproximateDiameter>(graph);

    double previous_cnt = 0.0;
    int diameter;
    for (int i=0; i<MAX_ITERATION_TIMES; i++){
        engine->signalAll();
        engine->start();
        double current_cnt = engine->map_reduce_vertices<double>(approximate_pair_number);

        if (i > 0 && current_cnt < previous_cnt * (1.0 + EPSILON)){
            diameter = i;
            break;
        }
    }

    std::cout << "Approximate Diameter is " << diameter << std::endl;

    delete engine;
    return EXIT_SUCCESS;
}
