#include <iostream>
#include <random>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <exception>
#include <algorithm>
#include <float.h>
#include <iterator>


#include <array>
#include <string>
#include <set>
#include <map>
#include <climits>

#include "thorup_zwick_oracle.hpp"
#include "lib/structs.hpp"
#include "lib/sssp.hpp"
#include "lib/heap.hpp"
#include "lib/std_methods.hpp"

/*
    We use unordered sets to get better average access time compared to sets:
    O(1) avg vs O(log n).

    Assumes vertex ids are below UINT_MAX (unsigned integer max value)!
*/


ThorupZwickOracle::ThorupZwickOracle(){}
ThorupZwickOracle::~ThorupZwickOracle(){}


int ThorupZwickOracle::preproc(Graph &G, unsigned const k){
    unsigned n = G.getVSize();
    if(k < 1 || n < 1)
        return 0;//TODO throw exception()
    const double p = pow(n,-1./k);
    std::map<unsigned, Vertex>& V = G.getV();
    witnessDistArraySize = k+1;
	witnessDistArray.resize(witnessDistArraySize);

    //vertex_id w -> {{v.id, dist(w,v)} | v\in V and dist(w,v) < dist(A_i+1,v)}
    // BMap C;
    witnessDistArray[k] = std::unique_ptr< WitnessMap > (new WitnessMap());
    for(auto vp : V){ //set all dists to infty (i.e. double max)
        std::pair< Vertex, double > tmpP = 
                                std::pair< Vertex, double >(vp.second, DBL_MAX);
        (*witnessDistArray[k].get())[vp.first] = tmpP;//{vp.second, DBL_MAX};
    }

    /*TODO this is ineffecient as each vertex added to A is copied. 
      Consider pointers to v\in V. (Though Vertex is a small object).*/
    //need to be sorted set, if std::set_difference is to be used!!
    // std::map<unsigned, Vertex> A[k+1]; //Malloc here (needed on heap?)
    std::vector< std::map<unsigned, Vertex> > A;
    A.resize(k+1);
    A[0] = std::map<unsigned, Vertex>(V);
    for(auto & vp : A[0]){
        vp.second.setWitnessSelf();
    }
    

    //1 <= i <= k-1
    constructASets(p, k-1, A);
    /*
        TODO probably should test/check if BALANCE is good enough, or make it DETERMINISTIC!! 
        (see p 13 in Thorup-Zwick paper) 
    */
    
    //init vector pointers
    for(auto vp : V){
        unsigned vid = vp.second.id;
        B[vid] = std::unique_ptr< VertexDistMap >(new VertexDistMap());
    }
    AdjacentEdgeMap* adjEMap = new AdjacentEdgeMap();
    getIncidentEdgesMap(G.getE(), *adjEMap);
    setupWitnessesBMap(G, *adjEMap, k-1, A);
    delete adjEMap;

    return 1;
}



double ThorupZwickOracle::query(Vertex u, Vertex v){
    unsigned i = 0;
    Vertex w = u;

    //while(w not in B(v) )
    while(B[v.id]->find(w.id) == B[v.id]->end() && ++i < witnessDistArraySize){
        //swap u,v
        Vertex tmp = u;
        u = v;
        v = tmp;

        w = (*witnessDistArray[i])[u.id].first;
    }
    if(i >= witnessDistArraySize)
        return DBL_MAX; // This should never happen in connected graph, but might in unconnected.

    //dist(w,u) = dist(p_i(u), u), dist (w,v) from B
    return (*witnessDistArray[i])[u.id].second + (*B[v.id])[w.id];
}


void ThorupZwickOracle::constructSpanner(Graph &GIn, Graph& spannerOut, 
                                         unsigned const k){
    unsigned n = GIn.getVSize();
    if(k < 1 || n < 1)
        return;//TODO throw exception()
    const double p = pow(n,-1./k);
    std::map<unsigned, Vertex>& V = GIn.getV();
    witnessDistArraySize = k+1;
	witnessDistArray.resize(witnessDistArraySize);
    witnessDistArray[k] =
        std::unique_ptr< WitnessMap >(new WitnessMap());
    for (auto vp : V) { //set all dists to infty (i.e. double max)
        std::pair< Vertex, double > tmpP =
            std::pair< Vertex, double >(vp.second, DBL_MAX);
        (*witnessDistArray[k].get())[vp.first] = tmpP;//{vp.second, DBL_MAX};
    }
    std::vector< std::map<unsigned, Vertex> > A;
    A.resize(k+1);
    A[0] = std::map<unsigned, Vertex>(V);

    AdjacentEdgeMap* adjEMap = new AdjacentEdgeMap();
    getIncidentEdgesMap(GIn.getE(), *adjEMap);
    
    
    for(auto & vp : A[0])
        vp.second.setWitnessSelf();
    constructASets(p, k-1, A);
    for(int i=k-1; i >= 0; --i){
        setWitnesses(GIn, A[i], i);
        //(*) if dist(A[i],v) == dist(A[i+1],v), set p_i(v) = p_i+1(v)
        for(auto vp : *(witnessDistArray[i]) ){
            if(i+1 <= k 
               && (*witnessDistArray[i+1])[vp.first].second == vp.second.second)
                (*witnessDistArray[i])[vp.first] 
                                        = (*witnessDistArray[i+1])[vp.first];
        }
        std::map<unsigned, Vertex> ATmp; // v.id -> v (because sets are stupid)
        set_difference(A[i].begin(), A[i].end(),
                       A[i+1].begin(), A[i+1].end(),
                       std::inserter(ATmp, ATmp.begin()));

        //compute modified trees/SSSP per w
        //for every w \in A_i - A_i+1
        for(auto w : ATmp){
            Graph* T = new Graph();
            AdjacentEdgeMap* adjEMapTmp = new AdjacentEdgeMap(*adjEMap);
            sssp::dijkstraModTree(GIn, *adjEMapTmp, *T, w.second,
                                  *witnessDistArray[i+1] );
            delete adjEMapTmp;
            //take union over all trees T(w); resulting edges are the spanner
            for(auto se : T->getE())
                spannerOut.getE()[se.first] = se.second;//
            delete T;
        }
        spannerOut.setV(V);
        spannerOut.setL(GIn.getL());
    }
    delete adjEMap;
}
