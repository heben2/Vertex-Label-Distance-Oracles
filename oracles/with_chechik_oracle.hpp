#include <unordered_map>

#include "oracle.hpp"
#include "lib/structs.hpp"
#include "chechik_oracle.hpp"
#include "restricted_chechik_oracle.hpp"
#include "lib/std_methods.hpp"


#ifndef WITH_CHECHIK_ORA
#define WITH_CHECHIK_ORA

class WithChechikOracle: public VertexLabelDistanceOracle {
private:
    ChechikOracle *OGS;
    RestrictedChechikOracle *RCO;
    bool RCODone = true;
public:
    WithChechikOracle();
    ~WithChechikOracle();
    //Assumes all labels are used and has integer ids in range [0,|L|).
    //n=|G.V|, k=const
    //returns 0 for unsuccessful run, any other for success
    int preproc(Graph &G, unsigned const k);

    //return (4k-5)-approximate distance between Vertex a and b
    double query(Vertex u, unsigned label);

    size_t getSize();
};

#endif /* WITH_CHECHIK_ORA */