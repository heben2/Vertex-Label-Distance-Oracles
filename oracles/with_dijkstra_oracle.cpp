#include <iostream>
#include <cmath>
#include "with_dijkstra_oracle.hpp"
// #include "thorup_zwick_oracle.hpp"
#include "spanner.hpp"
#include "lib/structs.hpp"
#include "lib/std_methods.hpp"
#include "lib/sssp.hpp"

#include <chrono>
WithDijkstraOracle::WithDijkstraOracle() {
    OGS = new ChechikOracle();
    lookupTableHandle = new std::unordered_map<unsigned, std::vector< double> >();
    lookupTable = *lookupTableHandle;
}

WithDijkstraOracle::~WithDijkstraOracle() {
    delete OGS;
    delete lookupTableHandle;
}

/*
Make S-set, sampling with probability p from V.
Compute E_S, make graph G_S = (V,E_S).
Compute witnesses and their distances for set S.
Make chechik-ora on G_S.

Make k'-spanner H on G.
Make |L|*|V|-lookup table.
Run Dijkstra on H for each v\in S. Store in table:
Each vertex u_λ\in V added to the SP tree using Dijkstra’s algorithm, add 
distance (v, u_λ) to the lookup table if the entry of λ is empty. 
Terminate when all labels for entry v have been filled out or when the 
complete SP tree is constructed.
*/
int WithDijkstraOracle::preproc(Graph &G, unsigned const k){
    unsigned const c = k > 2? k-1 : k; //k'
    unsigned n = G.getVSize();
    unsigned l = G.getLSize();
    if(k < 3 || n < 1 || l < 1)
        return 0;
    std::map<unsigned, Vertex>& V = G.getV();
    std::unordered_map<std::string, Edge>& E = G.getE();
    std::unordered_set<unsigned>& L = G.getL();
    /* TODO if equation (3) holds, we can choose better/optimal p
       THIS CAN BE TESTED HERE if we know k'
       Might be better to assign p from caller -- to test for other p's, then only test if within size bound
    */
    double p = k*std::pow(l, 1./k-1.);


    if (l <= std::pow(n, double(k) / (2.*k - 1.))) {
        bool const b = pow(l, 1. - 2. / (3.*k)) <= pow(k, 2. / 3.)*pow(c, 1. / 3.)*pow(n, 1. / 3. + 1. / (3.*c));
        std::cout << "small l, optimal p? " << b << std::endl;
        if (b) {
            p = pow(n, -1. / 2. - 1. / (2.*c))*sqrt(l)*pow(c,-1./2.);
        }
    }
    else {
        bool const b = pow(n, double(k)/ (4.* k - 2.) - 1./2. - 1./(2. * c)) / sqrt(c) <= k*pow(l, 1. /double(k) - 1.);
        std::cout << "large l, optimal p? " << b << std::endl;
        if (b) {
            p = pow(n, double(k) / (4.*k - 2.) - 1. / 2. - 1. / (2.*c)) * pow(c, -1. / 2.);// sqrt(l) / sqrt(c);
        }
    }

    std::cout << currentDateTime() << " " << "WithDijkstraOracle::preproc begin, k = " << k << ", p = " << p << std::endl;
    std::map<unsigned, Vertex> S;
    makeRandSubset(p, V, S);


    //get/set witnesses
    witnessDistArraySize = 1;
	witnessDistArray.resize(witnessDistArraySize);
    setWitnesses(G, S, 0);

    //make E_S = UNION_v\in V E_S(v)
    //E_S(v) = the set of edges incident to v with weight < dist_G(v, p_S(v)).
    std::unordered_map<std::string, Edge> ES;
    getEdgeSMap(E, *witnessDistArray[0], ES);


    //Make Chechik's oracle on G_S
    Graph* GS = new Graph(V, ES, L);
    OGS->preproc(*GS, k);
    delete GS;


    //Make k'-spanner of Benswana and Sen
    // Graph* H = new Graph();
    Graph* H = new Graph(); //TODO place in mem?
    auto t1 = std::chrono::high_resolution_clock::now();
    // ThorupZwickOracle TZO;
    // TZO.constructSpanner(G, *H, c);
    spanner::constructSpanner(G, *H, c);
    auto t2 = std::chrono::high_resolution_clock::now();

    AdjacentEdgeMap* adjHEMap = new AdjacentEdgeMap();
    getIncidentEdgesMap(H->getE(), *adjHEMap);
    lookupTable.reserve(S.size());
    for(auto& v : S){
        AdjacentEdgeMap* adjHEMapTmp = new AdjacentEdgeMap(*adjHEMap);
        sssp::dijkstraSourceLabelsDists(*H, *adjHEMapTmp, v.second,
                                        lookupTable[v.first]);
        delete adjHEMapTmp;
    }
    delete adjHEMap;
    delete H;
    return 1;
}

//d_G(u, p_S(u)) and the distance d_H(p_S(u), λ)
double WithDijkstraOracle::query(Vertex u, unsigned label){
    double dist1 = OGS->query(u,label);
    Vertex witness;
    double distWitness;
    std::tie(witness,distWitness) = (*witnessDistArray[0])[u.id];
    double dist2 = distWitness + lookupTable[witness.id][label];

    return std::min(dist1, dist2);
}




size_t WithDijkstraOracle::getSize(){
    size_t s = 0;
    for(auto e : lookupTable)
        s += e.second.size();

    return std::max(s, OGS->getSize());
}