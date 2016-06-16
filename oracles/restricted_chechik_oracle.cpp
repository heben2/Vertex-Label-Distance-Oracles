#include <iostream>
#include <float.h>
#include <cmath>

#include "restricted_chechik_oracle.hpp"
#include "lib/structs.hpp"
#include "lib/sssp.hpp"
#include "lib/std_methods.hpp"




/*
    compute lookup table T = |S|\ell, init with infty distances.
    For each p_S(v), fill out 
    T[p_S(v)][label(v)] = min(T[p_S(v)][label(v)], dist(p_s(v), v))

    Let G_D = H (*copy)
    Assign dummy-labels (unsigned max) to all vertices of V_D-S. (not added to L)
    
    For each vertex u\in S and each label in L, 
    add vertex v with label and edge (u, v) with weight (T[u][label]) to G_D.
    Let D be set of all such added vertices
    Let S' = S union D
    
    Construct k A-sets from S' (A[k-1] to A[0] = S')
    construct bunches and all as in Chechik's oracle from this.


    Store only:
    - B[label] for every label in L (excluding dummy-label)
    - Compute A[k-1] distances to all labels
    - pivots for vertices in S (note this is not vertices in S')
*/
//k = k''
//superk = k
int RestrictedChechikOracle::preproc(Graph& G, Graph& H, 
        std::map<unsigned, Vertex> S, unsigned const k, unsigned const superk,
        WitnessMap& superWitnessDistArray){
    unsigned n = G.getVSize();
    unsigned labelSize = G.getLSize();
    if(k < 2 || n < 1 || labelSize < 1 || superWitnessDistArray.size() < 1
       || S.size() < 1){
        std::cout << "RestrictedChechikOracle::preproc bad input: superWitnessDistArray.size() = " << superWitnessDistArray.size() << ", S.size() = " << S.size() << ", labelSize = " << labelSize << ", n = " << n << ", k = " << k << std::endl;
        return 0;
    }
    std::map<unsigned, Vertex>& V = G.getV();
    std::unordered_map<std::string, Edge>& E = G.getE();
    std::unordered_set<unsigned>& L = G.getL();

    double const p = std::pow(labelSize, -1./k);

    checkAllIndices = labelSize > std::pow(n,double(k)/(2.*k-1.));


    //init lookup table with 'infty' dists
    std::unordered_map< unsigned, std::vector< double > > T;
    for(auto& vp : S)
        T[vp.first].resize(labelSize, DBL_MAX);


    //update T with witness dists
    //Could be optimised by doing while finding witness dists
    for(auto& wdp : superWitnessDistArray)
        T[wdp.second.first.id][V[wdp.first].label] 
                        = std::min(T[wdp.second.first.id][V[wdp.first].label], 
                                   wdp.second.second);

    //init graph D
    unsigned dummyLabel = L.size();
    Graph* D = new Graph(H); //copy H
    VertexMap& VD = D->getV();
    EdgeMap& ED = D->getE();
    for(auto& vp : VD){
        unsigned vid = vp.first;
        if(S.find(vid) == S.end())
            vp.second.label = dummyLabel;
    }

    //Make new nodes. This is errorprone, as we just choose new ids from VD.size().
    //Also make S' (here called W)
    std::map<unsigned, Vertex> W;
    unsigned maxVId = VD.size();
    for(auto& vp : T){
        std::vector< double >& labelList = vp.second;
        for(unsigned l = 0; l < labelSize; ++l){
            double dist = labelList[l];
            if(dist < DBL_MAX){
                Vertex s = Vertex(maxVId++, l);
                Edge e = Edge(VD[vp.first], s, dist);
                VD[s.id] = s;
                W[s.id] = s;
                ED[e.getKey()] = e;
            }
        }
        W[vp.first] = S[vp.first];
    }

    //Init stuff as in Chechik's
    witnessDistArraySize = k;
	witnessDistArray.resize(witnessDistArraySize);
    BMap BLabel;
    for(auto l : L){
        BLabel[l] = std::unique_ptr< VertexDistMap >(new VertexDistMap());
    }
    for(auto vp : VD){
        B[vp.second.id] = std::unique_ptr< VertexDistMap >(new VertexDistMap());
    }

    // Construct k A-sets from S' (A[k-1] to A[0] = S')
    std::vector< std::map<unsigned, Vertex> > A;
    A.resize(k);
    A[0] = std::map<unsigned, Vertex>(S);
    for(auto & vp : A[0]){
        vp.second.setWitnessSelf();
    }

    constructASets(p, k-1, A);
    setWitnesses(*D, A[k-1], k-1);

    // construct bunches and all as in Chechik's oracle from this.
    AdjacentEdgeMap* adjEDMap = new AdjacentEdgeMap();
    getIncidentEdgesMap(D->getE(), *adjEDMap);
    setupWitnessesBMap(*D, *adjEDMap, k-2, A);


    if(checkAllIndices){ //if l>n^{k/(2k-1)}, run over A_k-1 instead over L
        // for every x\in B[label] store dist(x,label) already in B[v]
        // BLabel[label] = union(B[v_label]) for v_label\in VD
        for(auto vp : VD){
            unsigned l = vp.second.label;
            if(l != dummyLabel)
                for(auto bp : *B[vp.first]){
                    updateMin<unsigned, double>(*BLabel[l], bp.first, 
                                                bp.second);
                }
        }

        for(auto & vp : A[k-1]){
            Vertex v = vp.second;
            Graph* GTmp = new Graph(D);
            std::vector<double> distVectOut;
            AdjacentEdgeMap* tmpEDMap = new AdjacentEdgeMap(*adjEDMap);
            sssp::dijkstraSourceLabelsDists(*GTmp, *tmpEDMap, v, distVectOut);
            delete GTmp;
            delete tmpEDMap;
            //update BLabel with distances from v to all labels
            for(unsigned l=0; l < labelSize; ++l)
                (*BLabel[l])[v.id] = distVectOut[l];
        }
        delete adjEDMap;
    } else{
        delete adjEDMap;
        for(auto& vp : VD){
            unsigned l = vp.second.label;
            if(l != dummyLabel)
                for(auto& bp : *B[vp.first])
                    (*BLabel[l])[bp.first] = -1;

        }

        //Need the vertice set for each label
        //label-vertex map. Takes more space (temporarily), but O(n) running time instead of O(nl)
        std::unordered_map<unsigned, std::map<unsigned, Vertex> > LV;
        for(auto& vp : VD){
            unsigned l = vp.second.label;
            if(l != dummyLabel)
                LV[l][vp.first] = vp.second;
        }

        for(auto l : L){ //store dist(v,l) for every label l\in L in A[k-1]
            //for every label l, add zero-edges to v_l \in VD then run dijkstra
            std::vector<Edge> zeroEdges;
            
            Vertex s = getZeroVertex(LV[l], zeroEdges); //VD[l]
            Graph* GTmp = new Graph(D); 
            (*GTmp->V)[s.id] = s;
            for(Edge e : zeroEdges){
                (*GTmp->E)[e.getKey()] = e;
            }
            AdjacentEdgeMap* adjEMap = new AdjacentEdgeMap();
            getIncidentEdgesMap(GTmp->getE(), *adjEMap);
            WitnessMap* tmpMap = sssp::dijkstra(*GTmp, *adjEMap, s);
            delete GTmp;
            delete adjEMap;
            for(auto blp : *BLabel[l] ){
                unsigned v = blp.first;
                (*BLabel[l])[v] = (*tmpMap)[v].second;
            }

            for(auto & vp : A[k-1]){
                unsigned v = vp.first;
                (*BLabel[l])[v] = (*tmpMap)[v].second;
            }
            delete tmpMap;
        }
    }

    //Construct clusters, inverse of B.
    for(auto& v : VD){
        unsigned vid = v.first;
        C[vid] = std::unique_ptr< VertexDistMap >(new VertexDistMap()); //Note this is a labelDistMap
        for(auto &u : *B[vid]){
            unsigned l = VD[u.first].label;
            if(l != dummyLabel)
                updateMin<unsigned, double>(*C[vid], l, u.second);
        }
    }
    for(auto & a : B) //clean up B -- need to be reassigned to BLabel
        a.second.reset();
    B = std::move(BLabel);
    delete D;

    return 1;
}



// double RestrictedChechikOracle::query(Vertex u, unsigned label){
// //TODO
//     return -1.;
// }