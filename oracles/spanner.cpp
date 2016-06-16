#include <vector>
#include <float.h>
#include <unordered_map>
#include <iterator>
#include <random>
#include <iostream>
#include <climits>
#include <unordered_set>

#include "spanner.hpp"
#include "lib/structs.hpp"
#include "lib/std_methods.hpp"


struct CVertex{
    Vertex v;

    //AdjacentEdges can just be unsigned (target_id), double (weight) -- all other info is already present
    std::unordered_map<unsigned, Edge> AdjacentEdges; //u.id -> edge(u.id,v)

    CVertex(Vertex& v):v(v){
        AdjacentEdges = std::unordered_map<unsigned, Edge>();
    }

    CVertex(){
        AdjacentEdges = std::unordered_map<unsigned, Edge>();
    }

};

struct Cluster{
    unsigned vid; //center id
    std::unordered_set<unsigned> c; //the cluster (ids)

    Cluster(Vertex& v):vid(v.id){
        c = std::unordered_set<unsigned>();
        c.insert(v.id);
    }
    Cluster(unsigned vid):vid(vid){
        c = std::unordered_set<unsigned>();
        c.insert(vid);
    }
    Cluster(){
        c = std::unordered_set<unsigned>();
    }
};

typedef std::unordered_map<unsigned, Cluster> ClusterList; //c.id -> c

typedef std::unordered_map<unsigned, CVertex> ClusterVertexMap; //v.id -> CVertex (with adjacency list)

typedef std::unordered_map<unsigned, unsigned> ClusterMap; //v.id -> cluster.center.id

ClusterVertexMap V; //V', like G.V but with each vertex having adjacency list

ClusterMap CMapPrevi; //v.id -> cluster.center.id for C[i-1] NOT R_i-1!!
ClusterMap CMap;
ClusterMap ClusterMapR; // for i, works as R_i

// EdgeMap ClusterEdgesPrevi; //epsilon[i-1]
EdgeMap ES, //E_S, output spanner edges
        E, //E'
        ClusterEdges; //epsilon[i]


///////////////////////////////// DEBUGGING
void printDebug(){
std::cout << "printDebug\n";
std::cout << "V.size() = " << V.size() << std::endl;
std::cout << "CMap.size() = " << CMap.size() << std::endl;

std::cout << "CMapPrevi.size() = " << CMapPrevi.size() << std::endl;
std::cout << "ClusterMapR.size() = " << ClusterMapR.size() << std::endl;
std::cout << "ES.size() = " << ES.size() << std::endl;
std::cout << "E.size() = " << E.size() << std::endl;
std::cout << "ClusterEdges.size() = " << ClusterEdges.size() << std::endl;
}
void printE(){
    std::cout << "E = [" << std::endl;
    for(auto& e : E){
        std::cout << "e.u.id = " << e.second.u.id << ", e.v.id = " << e.second.v.id << ", weight = " << e.second.weight << std::endl;
    }
    std::cout << "]" << std::endl;
}

/////////////////////////////////

//Used if adjacent edges must represent E'
void removeAdjEdge(unsigned uid, unsigned vid){
    if(V.find(uid) != V.end())
        V[uid].AdjacentEdges.erase(vid);
    if(V.find(vid) != V.end())
        V[vid].AdjacentEdges.erase(uid);
}


/* Discard from E' edges E'(vid,c), that is, edges between vid and cluster with 
   center cid.
*/
void discardEdges(unsigned vid, unsigned cid, ClusterMap cm){
    for(auto& e : V[vid].AdjacentEdges){
        if(inMap<unsigned,unsigned>(cm, e.first) && cm[e.first] == cid){
            auto it = E.find(e.second.getKey());
            if(it != E.end()){
                E.erase(it);
            }
        }
    }
}


/* Get lightest edges from E' between vid and each cluster of cm.
   For each such edge, add it to ES and discard from E' all edges of E' between 
   vid and relavant cluster.
*/
void findAdjClusters(unsigned vid, ClusterMap cm, double w){
// std::cout << "spanner::findAdjClusters begin" << std::endl;
    std::unordered_map<unsigned, Edge> adjCLightestEdges;
    for(auto& ep : V[vid].AdjacentEdges){
        unsigned uid = ep.first;
        Edge& e = ep.second;

        if(inMap<unsigned,unsigned>(cm, uid) 
           && inMap<std::string,Edge>(E, e.getKey()) 
           && e.weight < w){
            updateMinEdge(adjCLightestEdges, cm[uid], e);
        }
    }
    for(auto& ep : adjCLightestEdges){
        ES[ep.second.getKey()] = ep.second;
        //Discard E'(v,c') from E'
        discardEdges(vid, ep.first, cm);
    }
}
 

// unsigned iterationDebugging = 0;
/* 1. Forming a sample of clusters:
*/
void sampleClusters(double const p){
    // For randomization
    std::random_device rd;
    std::ranlux48_base e2(rd());
    std::uniform_real_distribution<> dist(0, 1);

    ClusterMapR.clear();
    CMapPrevi = CMap;
    CMap.clear();
    EdgeMap ClusterEdgesPrevi = ClusterEdges;
    ClusterEdges.clear();
    ClusterList tmpClusterList;

    //get cluster centers
    for(auto kv : CMapPrevi){
        if(tmpClusterList.find(kv.second) != tmpClusterList.end())
            tmpClusterList[kv.second].c.insert(kv.first);
        else{
            tmpClusterList[kv.second] = Cluster(kv.second);
            tmpClusterList[kv.second].c.insert(kv.first);
        }
    }

    while(p > 0 && ClusterMapR.size() == 0) //TODO does it work without this?
        for(auto& cp : tmpClusterList){ //O(n)
            Cluster const& c = cp.second;
            if(dist(e2) <= p){
                for(auto& vid : c.c){
                    ClusterMapR[vid] = c.vid;
                }
            }
        }
    CMap = ClusterMap(ClusterMapR);

    //ClusterEdges[i] is defined of edges of ClusterEdges[i-1] defining C[i].
    for(auto& ep : ClusterEdgesPrevi){
        Edge const& e = ep.second;
        unsigned vid = e.v.id,
                 uid = e.u.id;
        if(ClusterMapR.find(vid) != ClusterMapR.end() 
           && ClusterMapR.find(uid) != ClusterMapR.end()){
            ClusterEdges[ep.first] = e;
        }
    }
}

/* 2+3. Finding nearest neighboring sampled cluster for each vertex and adding
   edges to spanner.
*/
void findNearestNeighboringCluster(){
// std::cout << "spanner::findNearestNeighboringCluster begin" << std::endl;

    /*For each v\in V' not belonging to any sampled cluster, compute its nearest
      neighboring cluster (if any) from R[i].
      I.e. the cluster incident on v with lightest edge
    */
    //v.id -> cluster c, edge(v,c)   //edge with min weight
    for(auto& vp : V){//V'
        unsigned vid = vp.first;
        if(!inMap<unsigned,unsigned>(ClusterMapR, vid)){ //not belonging to any sampled cluster
            //get nearest neighboring sampled cluster
            //scan adjecency list
            Edge eMinWeight = Edge(Vertex(UINT_MAX-1),Vertex(UINT_MAX),DBL_MAX);
            unsigned cid;
            for(auto& e : vp.second.AdjacentEdges){
                //get cluster of
                if(inMap<unsigned, unsigned>(ClusterMapR, e.first)
                   && eMinWeight.weight > e.second.weight){ //edge to cluster
                    cid = ClusterMapR[e.first];
                    eMinWeight = e.second;
                }
            }
           
            /*(3) Adding edges to the spanner */
            if(eMinWeight.weight < DBL_MAX){ //(b) if not dummy then v has neighboring cluster
                //add eMinWeight to ES and Epsilon[i]
                std::string ek = eMinWeight.getKey();
                ES[ek] = eMinWeight;
                ClusterEdges[ek] = eMinWeight;
                CMap[vid] = cid;

                //discard E'(v,c) from E'
                discardEdges(vid, cid, ClusterMapR);

                /* For each cluster c'\in Cprevi adjecent to v with 
                   weight < eMinWeight, add the least weight edge from E'(v,c')
                   to ES, then remove E'(v,c') from E'. 
                */
                findAdjClusters(vid, CMapPrevi, eMinWeight.weight);
            } else{ //(a), not adjecent to any sampled cluster
                findAdjClusters(vid, CMapPrevi, DBL_MAX);
            }
        }
    }
}


/* 4. Removing intra-cluster edges
   Remove all edges with both endpoints belonging to same cluster of C from E'
*/
void removeIntraClusterEdges(){
    for(auto it = E.begin(); it != E.end(); ){
        unsigned uid = it->second.u.id,
                 vid = it->second.v.id;
        if(CMap.find(uid) != CMap.end() 
           && CMap.find(vid) != CMap.end()
           && CMap[uid] == CMap[vid]){
            it = E.erase(it);
        } else{
            ++it;
        }
    }
}


/*  Phase 1 Forming the clusters
*/
void phase1(Graph& G, unsigned const k){
    // EdgeMap& E = G.getE();
    // VertexMap& V = G.getV();
    unsigned n = G.getVSize();
    double const p = pow(n,-1./k);
    
    //init
    ES = EdgeMap(); //E_S, output spanner edges
    E = EdgeMap(G.getE()); //E'
    ClusterEdges = EdgeMap();

    for(auto& vp : G.getV()){
        unsigned vid = vp.first;
        Vertex& v = vp.second;
        // C[vid] = Cluster(v);
        CMap[vid] = vid;
        V[vid] = CVertex(v);
        ClusterMapR[vid] = vid;
    }
    for(auto& ep : E){ //O(m) (E == G.getE())
        Edge& e = ep.second;
        V[e.v.id].AdjacentEdges[e.u.id] = e;
        V[e.u.id].AdjacentEdges[e.v.id] = e;
    }


    for(unsigned i=1; i<k; ++i){ //k-1 iterations
        sampleClusters(p);
        findNearestNeighboringCluster();
        removeIntraClusterEdges();
    }
}

/*   Phase 2: Vertex-Cluster joining
  For each v\in V' and each cluster c\in C_k-1
  add the least weight edge from set E'(v,c) to the spanner E_S,
  and discard E'(v,c) from E'.
*/
void phase2(){
    for(auto& vp : V){
        std::unordered_map<unsigned, Edge> adjCLightestEdges;
        unsigned vid = vp.first;
        for(auto& ep : vp.second.AdjacentEdges){ //TODO make sure adjecency represents E!!!!!
            Edge& e = ep.second;
            unsigned uid = ep.first;
            if(inMap<std::string,Edge>(E, e.getKey()) 
               && inMap<unsigned,unsigned>(CMap, uid)){
                unsigned cid = CMap[uid];
                updateMinEdge(adjCLightestEdges, CMap[uid], e);
            }
        }
        for(auto& ep : adjCLightestEdges){
            ES[ep.second.getKey()] = ep.second;
            /* Discard E'(v,c') from E'
               Maybe needed to hold bounds on running time??? 
               Nope, as we do not loop on E (we loop over adjecencylist) 
               and map is constant access time.
            */
            // discardEdges(vid, ep.first, CMap); 
        }
    }
}


/* (2k-1)-spanner based on Baswana and Sen's spanner.
   expected O(km) construction time and O(kn^{1+1/k}) edges in spanner.
*/

int spanner::constructSpanner(Graph &G, Graph& spanner, unsigned const k){
    unsigned n = G.getVSize();
    if(k < 1 || n < 1)
        return 0;//TODO throw exception()

    phase1(G, k);
    phase2();
    spanner.reset(G.getV(), ES, G.getL());

    return 1;
}