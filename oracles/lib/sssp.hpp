#include <map>
#include <vector>
#include "structs.hpp"
/*
  Return map vertex id -> {witness, dist}, where dist is the shortest path 
  distances to s.

  If s only has zero-edges, dist is also distance to the witness, i.e.:
  Return map from id (unsigned) of vertex v to {p(v),distance(v,p(v))}, 
  p(v) is the witness of v.

  Assumes maximum distances between vertices is of size int.
  Updates parent attribute of vertex references with pointer to actual parent, 
  such that SSSP tree is formed. Root (s) has NULL pointer.
*/

#ifndef SSSP_H
#define SSSP_H


namespace sssp {
    WitnessMap* dijkstra(Graph& G, AdjacentEdgeMap& adjEMap, Vertex s);

    WitnessMap* dijkstraMod(Graph& G, AdjacentEdgeMap& adjEMap, Vertex s, 
                            WitnessMap& ADistMap);

    /* 
    Constructs a shortest path tree from s in TOut.
    The tree TOut is spanning cluster 
    C(s) = {v\in V | dist(s,v) < dist(A[i+1],v)}
    where A[i+1] = AiDistMap.
    */
    void dijkstraModTree(Graph& GIn, AdjacentEdgeMap& adjEMap, Graph& TOut, 
                         Vertex s, WitnessMap& ADistMap);


    /*
    Each vertex u_λ\in V added to the SP tree using Dijkstra’s algorithm, add 
    distance (v, u_λ) to the lookup table if the entry of λ is empty. 
    Terminate when all labels for entry v have been filled out or when the 
    complete SP tree is constructed.
    */
    /*
    For s in G.V, construct SP tree using Dijkstra's and add dist(v, u_λ) to the
    output vector if the entry of λ is empty (i.e. the shortest distance).
    Terminate when all labels for entry s have been filled out or when the 
    complete SP tree is constructed.
    Assumes labels are unsigned from 0 to |L|-1, i.e. can work as indices for 
    distVectOut.
    */
    void dijkstraSourceLabelsDists(Graph& G, AdjacentEdgeMap& adjEMap, Vertex s, 
                                   std::vector<double>& distVectOut);

}

#endif /* SSSP_H */