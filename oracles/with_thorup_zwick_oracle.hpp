#include <unordered_map>

#include "oracle.hpp"
#include "lib/structs.hpp"
#include "chechik_oracle.hpp"
#include "lib/std_methods.hpp"


#ifndef WITH_THORUP_ZWICK_ORA
#define WITH_THORUP_ZWICK_ORA

/* If k <= 1, then k'=k. Else k' = k-1.
*/
class WithThorupZwickOracle: public VertexLabelDistanceOracle {
private:
    ChechikOracle* OGS;
    std::unordered_map<unsigned, std::vector< double> > *lookupTableHandle;
    std::unordered_map<unsigned, std::vector< double> > lookupTable;

public:
    //Assumes all labels are used and has integer ids in range [0,|L|).
    //n=|G.V|, k=const
    WithThorupZwickOracle();
    ~WithThorupZwickOracle();

    int preproc(Graph &G, unsigned const k);

    //return (4k-5)-approximate distance between Vertex a and b
    double query(Vertex u, unsigned label);

    size_t getSize();
};

#endif /* WITH_THORUP_ZWICK_ORA */