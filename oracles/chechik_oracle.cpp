#include <iostream>
#include <float.h>
#include <cmath>

#include "chechik_oracle.hpp"
#include "lib/structs.hpp"
#include "lib/sssp.hpp"
#include "lib/std_methods.hpp"



ChechikOracle::ChechikOracle(){}
ChechikOracle::~ChechikOracle(){}

int ChechikOracle::preproc(Graph &G, unsigned const k){
std::cout << currentDateTime() << " " << "ChechikOracle::preproc begin n = " << G.getVSize() << std::endl;
    unsigned const n = G.getVSize(), labelSize = G.getLSize();
    if(k < 2 || n < 1 || labelSize < 1)
        return 0;//TODO throw exception()

    checkAllIndices = labelSize > std::pow(n,double(k)/(2.*k-1.)); //for query
    const double p = std::pow(labelSize,-1./k);
    std::map<unsigned, Vertex>& V = G.getV();
    std::unordered_set<unsigned>& L = G.getL();

    std::cout << "ChechikOracle::preproc p =" << p<< std::endl;


    witnessDistArraySize = k;
	witnessDistArray.resize(witnessDistArraySize);

    BMap BLabel;
    for(auto l : L){
        BLabel[l] = std::unique_ptr< VertexDistMap >(new VertexDistMap());
    }
    for(auto vp : V){
        B[vp.second.id] = std::unique_ptr< VertexDistMap >(new VertexDistMap());
    }

    /*TODO this is ineffecient as each vertex added to A is copied. 
      Consider pointers to v\in V. (Though Vertex is a small object).*/
    //need to be sorted set, if std::set_difference is to be used!!
    // std::map<unsigned, Vertex> A[k]; //Malloc here (needed on heap?)
    std::vector< std::map<unsigned, Vertex> > A;
    A.resize(k);
    A[0] = std::map<unsigned, Vertex>(V);
    for(auto & vp : A[0]){
        vp.second.setWitnessSelf();
    }
    constructASets(p, k-1, A);
    
    setWitnesses(G, A[k-1], k-1); //set witnesses for last k (must not affect bunches B, so we do it seperately)

    AdjacentEdgeMap* adjEMap = new AdjacentEdgeMap();;
    getIncidentEdgesMap(G.getE(), *adjEMap);
    setupWitnessesBMap(G, *adjEMap, k-2, A);

    //Compute B(lambda) aka BLabel along with distances from A[k-1] to all labels.
    //Take special care of A[k-1], that is, compute distances to all labels. Just add this to BLabel.
    if(checkAllIndices){ //if l>n^{k/(2k-1)}, run over A_k-1 instead over L
        // for every x\in B[label] store dist(x,label) already in B[v]
        // BLabel[label] = union(B[v_label]) for v_label\in V
        for(auto vp : V){
            unsigned l = vp.second.label;
            for(auto bp : *B[vp.first]){
                updateMin<unsigned, double>(*BLabel[l], bp.first, bp.second);
            }
        }
        for(auto & vp : A[k-1]){
            Vertex v = vp.second;
            std::vector<double> distVectOut;
            AdjacentEdgeMap* tmpEM = new  AdjacentEdgeMap(*adjEMap);
            sssp::dijkstraSourceLabelsDists(G, *tmpEM, v, distVectOut);
            delete tmpEM;
            //update BLabel with distances from v to all labels
            for(unsigned l=0; l < labelSize; ++l)
                (*BLabel[l])[v.id] = distVectOut[l];
        }
        delete adjEMap;
    } else{
        delete adjEMap;
        //init BLabel[label] = union(B[v]) for v\in V_label
        for(auto vp : V){
            unsigned l = vp.second.label;
            for(auto bp : *B[vp.first]){
                (*BLabel[l])[bp.first] = DBL_MAX;//-1;
            }
        }

        //Need the vertice set for each label
        //label-vertex map. Takes more space (temporarily), but O(n) running time instead of O(nl)
        std::unordered_map<unsigned, VertexMap > LV;
        for(auto& vp : V){
            LV[vp.second.label][vp.first] = vp.second;
        }

        for(auto l : L){ //store dist(v,l) for every label l\in L in A[k-1]
            //for every label l, add zero-edges to v_l \in V then run dijkstra
            std::vector<Edge> zeroEdges;
            
            Vertex s = getZeroVertex(LV[l], zeroEdges); //V[l]
            Graph* GTmp = new Graph(G); 
            (*GTmp->V)[s.id] = s;
            for(Edge e : zeroEdges){
                (*GTmp->E)[e.getKey()] = e;
            }
            adjEMap = new AdjacentEdgeMap();
            getIncidentEdgesMap(GTmp->getE(), *adjEMap);
            WitnessMap* tmpMap = sssp::dijkstra(*GTmp, *adjEMap, s);
            delete adjEMap;
            delete GTmp;

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


//////////////////////////////////////////////
    //Construct clusters, inverse of B.
    for(auto& v : V){
        unsigned vid = v.first;
        C[vid] = std::unique_ptr< VertexDistMap >(new VertexDistMap()); //Note this is a labelDistMap
        for(auto &u : *B[vid]){
            unsigned l = V[u.first].label;
            updateMin<unsigned, double>(*C[vid], l, u.second);
        }
    }


    for(auto & a : B) //clean up B -- need to be reassigned to BLabel
        a.second.reset();
    B = std::move(BLabel);

    return 1;
}



double ChechikOracle::query(Vertex u, unsigned label){
    int i = -1;
    double dist1 = DBL_MAX;
    double dist2 = DBL_MAX;

    Vertex witness = u;
    double distWitness = 0;

    if(C[u.id]->find(label) != C[u.id]->end())
        dist1 = (*C[u.id])[label];

    //TODO should always expect witnessDistArray to contain all vertices, but check is nice I suppose (especially for bad input)
    //instead of splitting this up, just include the last iteration with i=k-1
    if(checkAllIndices){ // When l > n^(k/(2k-1))
        while(++i < witnessDistArraySize && (*witnessDistArray[i]).find(u.id) 
                                         != (*witnessDistArray[i]).end()) {
            std::tie(witness,distWitness) = (*witnessDistArray[i])[u.id];
            if(B[label]->find(witness.id) != B[label]->end())
                dist2 = std::min(dist2, distWitness + (*B[label])[witness.id]);
        }
    } else{
        auto it = B[label]->find(witness.id);
        while(it == B[label]->end() && ++i < witnessDistArraySize
              && (*witnessDistArray[i]).find(u.id) 
              != (*witnessDistArray[i]).end()){
            std::tie(witness,distWitness) = (*witnessDistArray[i])[u.id];
            it = B[label]->find(witness.id);
        }
        if(distWitness < DBL_MAX && i != witnessDistArraySize ) //&& it != B[label]->end() should never happen
            dist2 = distWitness + it->second;
    }

    return std::min(dist1, dist2);
}

size_t ChechikOracle::getSize(){
    size_t s = 0;
    for(auto& b : B){
        s += b.second->size();
    }
    return s;
}