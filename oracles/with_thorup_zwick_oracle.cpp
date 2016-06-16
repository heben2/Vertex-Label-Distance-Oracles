#include <iostream>
#include <cmath>
#include <float.h>
#include <algorithm>
#include "with_thorup_zwick_oracle.hpp"
#include "thorup_zwick_oracle.hpp"
#include "spanner.hpp"
#include "lib/structs.hpp"
#include "lib/std_methods.hpp"
#include "lib/sssp.hpp"

#include <chrono>



WithThorupZwickOracle::WithThorupZwickOracle() {
    OGS = new ChechikOracle();
    lookupTableHandle = new std::unordered_map<unsigned, std::vector< double> >();
    lookupTable = *lookupTableHandle;
}
WithThorupZwickOracle::~WithThorupZwickOracle() {
    delete OGS;
    delete lookupTableHandle;
}

int WithThorupZwickOracle::preproc(Graph &G, unsigned const k){
    unsigned c = floor(double(sqrt(4*k-7) + 2)/4.); //k' for spanner
    unsigned kk = c; //k'' for TZ-oracle
    unsigned n = G.getVSize();
    unsigned l = G.getLSize();
    auto boundHoldsLambda = [&](unsigned kk_in, unsigned c_in){
            return c_in >= 1 && kk_in >= 1 
                   && 2 + 4*(2*kk_in - 1)*(2*c_in-1) <= 4*k-5;};
    if(k < 4 || n < 1 || l < 1 || !boundHoldsLambda(kk, c))
        return 0;

    //fill out what is left of the stretch bound
    while(boundHoldsLambda(kk, c+1)) //decrease edge size of spanner
        ++c;

    double const p = 1./(std::cbrt(kk)*std::cbrt(n));
    std::map<unsigned, Vertex>& V = G.getV();
    std::unordered_map<std::string, Edge>& E = G.getE();
    std::unordered_set<unsigned>& L = G.getL();

    // std::map<unsigned, Vertex> S = *new std::map<unsigned, Vertex>();
    std::map<unsigned, Vertex> S;
    makeRandSubset(p, V, S);


    witnessDistArraySize = 1;
	witnessDistArray.resize(witnessDistArraySize);

    std::cout << S.size() << std::endl;

    //init lookup table
    lookupTable.reserve(S.size());
    for(auto& v : S){
        lookupTable[v.first].resize(l, DBL_MAX);
        lookupTable[v.first][v.second.label] = 0;
    }
    /* 
    For each v\in V and lable lambda, when computing witnesses p_S(v) and 
    dist_G(v,p_S(v)), update 
    lookupTable[p_S(v)][lambda] = min(lookupTable[p_S(v)][lambda], dist(v,p_S(v)));
    */
    setWitnesses(G, S, 0);

    for(auto& vp : *witnessDistArray[0]){
        unsigned l = V[vp.first].label;

        lookupTable[vp.second.first.id][l] = 
                    std::min(lookupTable[vp.second.first.id][l], vp.second.second);
    }

    //make E_S = UNION_v\in V E_S(v)
    //E_S(v) = the set of edges incident to v with weight < dist_G(v, p_S(v)).
    std::unordered_map<std::string, Edge> ES;
    getEdgeSMap(E, *witnessDistArray[0], ES);


    //Make Chechik's oracle on G_S
    Graph GS(V, ES, L);
    OGS->preproc(GS, k);

    //Graph* H = new Graph();
    Graph* H = new Graph();
    auto t1 = std::chrono::high_resolution_clock::now();
    spanner::constructSpanner(G, *H, c);
    auto t2 = std::chrono::high_resolution_clock::now();


    ThorupZwickOracle* TZO = new ThorupZwickOracle();
    TZO->preproc(*H, kk);
    delete H;

    for(auto& vp : S)
        for(auto& up : S){
            double d = TZO->query(vp.first, up.first);
            for(auto& l : L){
                lookupTable[vp.first][l] = std::min(lookupTable[vp.first][l], 
                                                    d+lookupTable[up.first][l]);
            }
        }
    delete TZO;

    return 1;
}


double WithThorupZwickOracle::query(Vertex u, unsigned label){
    double dist1 = OGS->query(u,label);
    Vertex witness;
    double distWitness;
    std::tie(witness,distWitness) = (*witnessDistArray[0])[u.id];
    double dist2 = distWitness + lookupTable[witness.id][label]; //overflow protection (no wrap, simply stops at max)
    return std::min(dist1, dist2);
}


size_t WithThorupZwickOracle::getSize(){
    size_t s = 0;
    for(auto e : lookupTable)
        s += e.second.size();

    return std::max(s, OGS->getSize()); //s + OGS.getSize();
}