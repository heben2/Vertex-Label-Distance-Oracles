#include <unordered_set>
#include <algorithm>
#include <array>

#include "oracle.hpp"
#include "lib/structs.hpp"


#ifndef CHECHIK_ORA
#define CHECHIK_ORA

//Note that B holds for label ids
class ChechikOracle: public virtual VertexLabelDistanceOracle {
protected:
    // // std::unique_ptr<BMap> B;
    // BMap B;
    // //i -> v.id -> p_i(v), dist(v,p_i(v))
    // std::unique_ptr< WitnessMap >* witnessDistArray;
    // size_t witnessDistArraySize;

    // BMap AkLabelDistMap;
    BMap C; //v.id -> label -> dist(v,label); inverse of B
    bool checkAllIndices; //when |L| is large
public:
    explicit ChechikOracle();
    ChechikOracle& operator=(const ChechikOracle&) = delete;
    ~ChechikOracle();
    
    //n=|G.V|, k=const
    int preproc(Graph &G, unsigned const k);

    //return (4k-5)-approximate distance between Vertex a and b
    double query(Vertex u, unsigned label);

    size_t getSize();

};

#endif /* CHECHIK_ORA */