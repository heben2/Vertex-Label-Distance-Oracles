#include <unordered_map>

#include "oracle.hpp"
#include "lib/structs.hpp"
#include "chechik_oracle.hpp"
#include "lib/std_methods.hpp"

#ifndef WITH_DIJKSTRA_ORA
#define WITH_DIJKSTRA_ORA

/* If k <= 1, then k'=k. Else k' = k-1.
*/
class WithDijkstraOracle: public VertexLabelDistanceOracle {
private:
    ChechikOracle *OGS;
    //TODO this should be fixed size of |S|*ell, not just a hash map!! (hash maps are possibly larger). But how to??
    std::unordered_map<unsigned, std::vector< double> > *lookupTableHandle;
    std::unordered_map<unsigned, std::vector< double> > lookupTable;

public:
    //Assumes all labels are used and has integer ids in range [0,|L|).
    //n=|G.V|, k=const
    WithDijkstraOracle();
    ~WithDijkstraOracle();

    int preproc(Graph &G, unsigned const k);

    //return (4k-5)-approximate distance between Vertex a and b
    double query(Vertex u, unsigned label);

    size_t getSize();
};

#endif /* WITH_DIJKSTRA_ORA */