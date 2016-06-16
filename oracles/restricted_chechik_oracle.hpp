#include <unordered_set>
#include <algorithm>
#include <array>

#include "oracle.hpp"
#include "chechik_oracle.hpp"
#include "lib/structs.hpp"


#ifndef RESTRICTED_CHECHIK_ORA
#define RESTRICTED_CHECHIK_ORA

/*
    A restricted Chechik's oracle used as a sub-oracle. Thus expects data 
    structures from its super oracle.
*/
class RestrictedChechikOracle: public ChechikOracle {
public:
    /*
        k = k'', superk = k (the original k for the original ora)
        superWitnessDistArray is the witness dist map for S on the original ora.
    */
    int preproc(Graph& G, Graph& H, std::map<unsigned, Vertex> S, 
                 unsigned const k, unsigned const superk, 
                 WitnessMap& superWitnessDistArray);

};

#endif /* RESTRICTED_CHECHIK_ORA */