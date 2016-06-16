#include <unordered_set>
#include <string>
#include <vector>
#include <map>
#include "structs.hpp"
#include <iostream>


#ifndef STD_METHODS_H
#define STD_METHODS_H

//construct vertex.id -> incident edges for all vertices and edges
void getIncidentEdgesMap(EdgeMap& E, AdjacentEdgeMap& adjEMap);

//retrieve and delete (from whole map) incident edges of vertex
std::unordered_set<Edge, Hash> retrieveIncidentEdges(AdjacentEdgeMap& EMap,
                                                     unsigned vid);


/* Simple string split functions.
   Copied from http://stackoverflow.com/questions/236129/split-a-string-in-c
*/
std::vector<std::string> &split(const std::string &s, char delim,
                                std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);


/* Fills up outMap with a random subset of inMap based on probability parameter
   p.
*/
void makeRandSubset(const double p, std::map<unsigned, Vertex>& inMap, 
                    std::map<unsigned, Vertex>& outMap);

/* Constructs A such that A[i] contains each element of A[i-1] independently 
   with probability p for 1<= i <= k. (using constructRandSubset)
   Expects A[0] to be non-empty.
*/
void constructASets(const double p, const unsigned k, 
                    std::vector< std::map<unsigned, Vertex> >& A);


/* Constructs a dummy-vertex with id=UINT_MAX, and creates zero-edges to all 
   vertices of A. These are placed in zeroEdges.

*/
Vertex getZeroVertex(std::map<unsigned, Vertex>& A, 
                     std::vector<Edge>& zeroEdges);


/* Get E_S(v): the set of edges incident to v with weight < dist_G(v, p_S(v)).
*/
void getEdgeSMap(std::unordered_map<std::string, Edge>& E,
                 WitnessMap& witnessDistMap,
                 std::unordered_map<std::string, Edge>& EOut);


void printG(Graph& G);


template <typename K, typename V>
void updateMin(std::unordered_map<K, V>& m, K key, V& val){
    if(m.find(key) != m.end())
        m[key] = std::min(m[key], val);
    else
        m[key] = val;
}
void updateMinEdge(std::unordered_map<unsigned, Edge>& m, unsigned key, 
                   Edge& e);


const std::string currentDateTime();

template <typename K, typename V>
bool inMap(std::unordered_map<K, V>& m, K key){
    return m.find(key) != m.end();
}
template <typename K, typename V>
bool inMap(std::map<K, V>& m, K key){
    return m.find(key) != m.end();
}

#endif /* STD_METHODS_H */