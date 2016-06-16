#include <unordered_set>
#include <algorithm>
#include <array>

#include "oracle.hpp"
#include "lib/structs.hpp"


#ifndef THORUP_ZWICK_ORA
#define THORUP_ZWICK_ORA

//vertex_id w -> {v.id -> dist(w,v) | v\in V and dist(w,v) < dist(A_i+1,v)}
// typedef std::unordered_map < unsigned, std::unique_ptr< 
//                                 std::unordered_map< unsigned, double> > > BMap;
//v.id -> u.id -> dist(u,v)
// typedef std::unordered_map < unsigned, std::unique_ptr< VertexDistMap > > BMap;

class ThorupZwickOracle: public DistanceOracle {
private:
    //reminder of what is here:
    // BMap B;
    //i -> v.id -> p_i(v), dist(v,p_i(v))
    // std::unique_ptr< WitnessMap >* witnessDistArray;
    // size_t witnessDistArraySize;

public:
    //init stuff
    explicit ThorupZwickOracle();
    ThorupZwickOracle& operator=(const ThorupZwickOracle&);
    ~ThorupZwickOracle();
    
    //n=|G.V|, k=const
    int preproc(Graph &G, unsigned const k);

    //return (2k-1)-approximate distance between Vertex a and b
    double query(Vertex u, Vertex v);

    //return spanner constructed from the preprocessing of graph G.
    void constructSpanner(Graph &GIn, Graph &spannerOut, unsigned const k);
};

#endif /* THORUP_ZWICK_ORA */