#include <iostream>
#include <map>
#include <unordered_map>
#include <float.h>
#include <iterator>

#include "oracle.hpp"
#include "lib/structs.hpp"
#include "lib/sssp.hpp"
#include "lib/std_methods.hpp"

#include <thread>
#include <chrono>

Oracle::~Oracle(){
    clear();
}

void Oracle::clear(){
    for(auto & a : B)
        a.second.reset();

    for(unsigned i=0; i < witnessDistArraySize; ++i){
        witnessDistArray[i].reset(nullptr);
    }
}


void Oracle::setWitnesses(Graph &G, std::map<unsigned, Vertex>& a, unsigned i){
    if (a.size() < G.getVSize()) { 
        std::vector<Edge> zeroEdges;
        Vertex s = getZeroVertex(a, zeroEdges);

        Graph* GTmp = new Graph(G);
        (*GTmp->V)[s.id] = s;
        for (Edge e : zeroEdges) {
            (*GTmp->E)[e.getKey()] = e;
        }
        for (auto & ap : a) {
            GTmp->getV()[ap.first].setWitnessSelf();
        }
        AdjacentEdgeMap* adjEMap = new AdjacentEdgeMap();;
        getIncidentEdgesMap(GTmp->getE(), *adjEMap);

        auto wmap = sssp::dijkstra(*GTmp, *adjEMap, s);
        delete GTmp;
        delete adjEMap;
        wmap->erase(s.id);

        //// set witnesses for those who witness them selves (are disconnected) and are not in a.
        for(auto& v : *wmap){
            unsigned vid = v.second.first.id;
            if (a.find(vid) == a.end() && a.size() > 0) {
                v.second.first = a.begin()->second;
            }
        }
        
        witnessDistArray[i].reset(wmap);
    } else { //(to speed up Djijkstra's for all V; they are their own witness for all vertices)
        std::map<unsigned, Vertex>& V = G.getV();
        WitnessMap* wMap = new WitnessMap();
        for (auto& vp : V) {
            (*wMap)[vp.second.id] = std::pair<Vertex, double>(vp.second, 0);
        }

        witnessDistArray[i].reset(wMap);
    }
}


void Oracle::setupWitnessesBMap(Graph &G, AdjacentEdgeMap& adjEMap, 
                                unsigned const k, 
                                std::vector< std::map<unsigned, Vertex> >& A){
    for(int i=k; i >= 0; --i){
        setWitnesses(G, A[i], i);
        if((*witnessDistArray[i]).size() > 0){
            for(auto vp : *(witnessDistArray[i]) ){
                if( (*witnessDistArray[i+1]).find(vp.first) != (*witnessDistArray[i+1]).end() && //i+1 <= k &&
                   (*witnessDistArray[i+1])[vp.first].second == vp.second.second)
                    (*witnessDistArray[i])[vp.first] 
                                            = (*witnessDistArray[i+1])[vp.first];
            }
        }

        /* For every w\in A[i] - A[i+1], compute 
           C(w)={v\in V | dist(w,v) < dist(A[i+1],v)}.
           Done by computing a modified sssp tree.
        */
        std::map<unsigned, Vertex> ATmp; // v.id -> v (because sets are stupid)
        set_difference(A[i].begin(), A[i].end(),
                       A[i+1].begin(), A[i+1].end(),
                       std::inserter(ATmp, ATmp.begin()));

        //compute modified trees/SSSP per w
        //for every w \in A_i - A_i+1
        for(auto& w : ATmp){
            AdjacentEdgeMap* adjEMapTmp = new AdjacentEdgeMap(adjEMap);
            WitnessMap* l = sssp::dijkstraMod(G, *adjEMapTmp, w.second,
                                              *witnessDistArray[i+1] );
            delete adjEMapTmp;
            for(auto& p : *l){
                if(p.second.second < DBL_MAX){ //distances of DBL_MAX is regarded as infty
                    (*B[p.first])[w.first] = p.second.second;
                }
            }
            delete l;
        }
    }

}