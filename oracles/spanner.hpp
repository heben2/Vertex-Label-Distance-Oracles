
#include "lib/structs.hpp"


#ifndef SPANNER
#define SPANNER

namespace spanner{
    /* (2k-1)-spanner based on Baswana and Sen's spanner.
       expected O(km) construction time and O(kn^{1+1/k}) edges in spanner.
    */
    int constructSpanner(Graph &G, Graph& spanner, unsigned const k);
}

#endif /* SPANNER */