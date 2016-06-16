#include <algorithm>
#include <unordered_set>
// #include <climits>
#include <float.h>
#include <map>
#include <functional>
#include <stdexcept>
#include <memory>

// #include <tuple>

#include "sssp.hpp"
#include "structs.hpp"
#include "std_methods.hpp"
#include "heap.hpp"


/*

  If we get into problems with stack size, move heap to dynamic mem (heap), 
  (i.e. new heap).

*/


void initHeap(std::map<unsigned, Vertex>& V, Heap<Vertex >& heap, Vertex& s){
    std::vector<std::pair<Vertex*, double> > inputVec;
    inputVec.reserve(V.size());
    for(auto & p : V){
        Vertex* v = &p.second;
        if(*v == s){
            inputVec.push_back({v, 0});
        } else
            inputVec.push_back({v, DBL_MAX});
    }
    heap.buildHeap(inputVec);
}
//expects s to point at s in V
void initHeap(Heap<Vertex >& heap, Vertex* s){
    std::vector<std::pair<Vertex*, double> > inputVec;
    inputVec.reserve(1);
    inputVec.push_back({s, 0});
    heap.buildHeap(inputVec);
}

////////////////////////// constructing distance map //////////////////////////

/*
  Returns a map with vertices and their their shortest path distance to s.
  Assumes maximum distances between vertices is of size int.
  Also updates vertex references with witness pointer, such that SSSP tree is 
  formed.
  
  Could have returned just a vertex-dist map, but instead returns a 
  vertex-{witness, dist} map for convenience.
*/
WitnessMap* SSSP(Graph &G, AdjacentEdgeMap& adjEMap, Vertex s,
                 std::function<bool (double, Vertex) > op){
    std::map<unsigned, Vertex>& V = G.getV();
    WitnessMap* retMap = new WitnessMap();
    Heap<Vertex > uvisitedVHeap = Heap<Vertex >(std::less<double>());
    initHeap(uvisitedVHeap, &V[s.id]);
    // Run SSSP
    while(!uvisitedVHeap.empty()){
        std::pair<Vertex,double> h = uvisitedVHeap.pop();

        Vertex& v = h.first;
        double weight = h.second;
        if(v.witness != nullptr)
            (*retMap)[v.id] = std::pair<Vertex, double>(*(v.witness), weight);
        else
            (*retMap)[v.id] = std::pair<Vertex, double>(v, weight);
        std::unordered_set<Edge, Hash> adjEdges = retrieveIncidentEdges(adjEMap,
                                                                        v.id);
        for(auto& e : adjEdges){
            Vertex& u = e.u == v? V[e.v.id] : V[e.u.id];
            double newWeight = weight+e.weight;
            if(op(newWeight, u)){
                if(!u.ptr){//new vertex
                    uvisitedVHeap.insert(&u, newWeight);
                    if(u.witness == NULL || u.id != u.witness->id ){
                        u.witness = v.witness;
                    }
                }else{ //if(newWeight < nu->key){
					//relax edge
                    Node<Vertex>* nu = reinterpret_cast<Node<Vertex>* >(u.ptr);
                    if(newWeight < nu->key){
    					uvisitedVHeap.alterKey(nu->index, newWeight); // makes "less" check for us
    					if(u.witness == NULL || u.id != u.witness->id ){
    						u.witness = v.witness;
    					}
                    }
				}
			}
        }
    }
    //The following is needed as we do not init with all vertices V with inft dist in heap (as dijkstraMod cannot run in this time)
    //fill out for all not connected to s with own witness and dist of 0
    if(retMap->size() < G.getVSize()){
        for(auto& vp : G.getV()){
            if(retMap->find(vp.first) == retMap->end()){
                (*retMap)[vp.first] = std::pair<Vertex, double>(vp.first, DBL_MAX);
            }
        }
    }
    return retMap;
}

/*
    Single Source Shortest-Paths problem. Implemented using Dijkstra and binary 
    heap.
*/
WitnessMap* sssp::dijkstra(Graph& G, AdjacentEdgeMap& adjEMap, Vertex s){
    auto l = [](double a, Vertex v) {return true;}; //dummy
    return SSSP(G, adjEMap, s, l);
}

/* 
   Modified version of Dijkstra's SSSP (above) based on Thorup-Zwick's article
   Approximate Distance Oracles [2005] (page 16).
   Relax edge (u,v) only if u,dist + (u,v).weight < delta(A_i+1,v).
   Assumes:
   ADistMap(v) = delta(A_i+1,v)
*/
WitnessMap* sssp::dijkstraMod(Graph& G, AdjacentEdgeMap& adjEMap, Vertex s,  WitnessMap& ADistMap){
    // auto l = [&AiDistMap](int a, Vertex v) {return a < AiDistMap[v];};
    std::function<bool (double, Vertex b) > l = 
        [&ADistMap](double a, Vertex v) {
            if(ADistMap.find(v.id) != ADistMap.end()){
                return a < ADistMap[v.id].second;
            } else{
                return true;
            }};
    return SSSP(G, adjEMap, s, l);
}

////////////////////// constructing vertex-label dist map //////////////////////
/*
Each vertex u_λ\in V added to the SP tree using Dijkstra’s algorithm, add 
distance (v, u_λ) to the lookup table if the entry of λ is empty. 
Terminate when all labels for entry v have been filled out or when the 
complete SP tree is constructed.
*/
void sssp::dijkstraSourceLabelsDists(Graph& G, AdjacentEdgeMap& adjEMap, 
                                     Vertex s, 
                                     std::vector<double>& distVectOut){
    std::map<unsigned, Vertex>& V = G.getV();
    unsigned l = G.getLSize();
    distVectOut.resize(l, DBL_MAX); //-1.);
    Heap<Vertex > uvisitedVHeap = Heap<Vertex >(std::less<double>());
    initHeap(uvisitedVHeap, &V[s.id]);
    unsigned filledDists = 0;
    // Run SSSP
    while(!uvisitedVHeap.empty()){
        std::pair<Vertex,double> h = uvisitedVHeap.pop();
        Vertex& v = h.first;
        double weight = h.second;
        if(weight < distVectOut[v.label]){
            distVectOut[v.label] = weight;
            ++filledDists;
            if(filledDists == l)
                break;
        }
        std::unordered_set<Edge, Hash> adjEdges = retrieveIncidentEdges(adjEMap,
                                                                        v.id);
        for(auto& e : adjEdges){
            Vertex& u = e.u == v? V[e.v.id] : V[e.u.id];
            double newWeight = weight+e.weight;
            //relax edge      
            if(!u.ptr){//new vertex
                uvisitedVHeap.insert(&u, newWeight);
            } else {
                Node<Vertex>* nu = reinterpret_cast<Node<Vertex>* >(u.ptr);
                // makes "less" check for us
                uvisitedVHeap.alterKey(nu->index, newWeight);
            }
        }
    }
}




////////////////////////////// constructing trees /////////////////////////////

void SSSPTree(Graph &GIn, AdjacentEdgeMap& adjEMap, Graph& TOut, Vertex s, 
              std::function<bool (double, Vertex) > op){
    std::map<unsigned, Vertex>& V = GIn.getV();
    TOut.setV(V);// = std::unique_ptr< std::map< unsigned, Vertex>(V) >;
    std::unordered_map<unsigned, Edge> vertexEdgeMap;

    Heap<Vertex > uvisitedVHeap = Heap<Vertex >(std::less<double>());
    //initialize heap
    initHeap(uvisitedVHeap, &V[s.id]);
    // Run SSSP
    while(!uvisitedVHeap.empty()){
        std::pair<Vertex,double> h = uvisitedVHeap.pop();
        Vertex& v = h.first;
        double weight = h.second;
        std::unordered_set<Edge, Hash> adjEdges = retrieveIncidentEdges(adjEMap,
                                                                        v.id);
        for(auto& e : adjEdges){
            Vertex& u = e.u == v? V[e.v.id] : V[e.u.id];
            double newWeight = weight+e.weight;
            if(op(newWeight, u)){
                vertexEdgeMap[u.id] = e;
                
                if(!u.ptr){//new vertex
                    uvisitedVHeap.insert(&u, newWeight);
                }else {
                    //relax edge
                    Node<Vertex>* nu = reinterpret_cast<Node<Vertex>* >(u.ptr);
                    uvisitedVHeap.alterKey(nu->index, newWeight); // makes "less" check for us
                }
            }
        }
    }
    //input the edges of the tree T(s)
    std::unordered_map<std::string, Edge>& E = TOut.getE();
    for(auto ve : vertexEdgeMap){
        E[ve.second.getKey()] = ve.second;
    }
}


void sssp::dijkstraModTree(Graph& GIn, AdjacentEdgeMap& adjEMap, Graph& TOut,
                           Vertex s, WitnessMap& ADistMap){
    std::function<bool (double, Vertex) > l = 
        [&ADistMap](double a, Vertex v) {
            if(ADistMap.find(v.id) != ADistMap.end()){
                return a < ADistMap[v.id].second;
            } else{
                return true;
            } };
    return SSSPTree(GIn, adjEMap, TOut, s, l);
}

