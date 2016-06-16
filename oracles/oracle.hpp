#include <map>
#include <vector>
#include <memory>
#include "lib/structs.hpp"
/* Abstract classes for distance oracles

*/
#ifndef ORACLE
#define ORACLE

class Oracle {
protected:
    //unsigned (e.g. v.id) -> u.id -> dist(u,v)
    BMap B;     // std::unique_ptr<BMap> B;
    //i -> v.id -> p_i(v), dist(v,p_i(v))
    std::vector< std::unique_ptr< WitnessMap > > witnessDistArray;
    size_t witnessDistArraySize = 0;

    /* For all v\in V, compute distances and witnesses
       Add source vertex s with zero-edges to all of A_i.
       Run Dijkstra's single-source shortest paths alg defined in DijkstraSSSP.
    */
    void setupWitnessesBMap(Graph &G, AdjacentEdgeMap& adjEMap, unsigned const k, 
                            std::vector< std::map<unsigned, Vertex> >& A);
    /* Set witnesses for witnessDistArray[i].
    */
    void setWitnesses(Graph &G, std::map<unsigned, Vertex>& a, unsigned i);
public:
    //return 0 for unsuccessful run, any other for success.
    virtual int preproc(Graph &G, unsigned const k) = 0;
    // virtual void preproc(Graph &G, unsigned const k);
    ~Oracle();
    void clear();
};


class DistanceOracle: public Oracle {
public:
    virtual double query(Vertex u, Vertex v) = 0;
};


class VertexLabelDistanceOracle: public Oracle {
public:
    virtual double query(Vertex u, unsigned labelId) = 0;
};

#endif /* ORACLE */
