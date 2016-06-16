#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <random>
#include <algorithm>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <climits>
#include <time.h>
#include <iterator>
#include "std_methods.hpp"

void getIncidentEdgesMap(EdgeMap& E, AdjacentEdgeMap& adjEMap){
    AdjacentEdgeMap out_E;
    for (auto& ep : E){
        // std::string& key = ep.first;
        Edge& e = ep.second;
        adjEMap[e.v.id].insert(e);
        adjEMap[e.u.id].insert(e);
    }
}

std::unordered_set<Edge, Hash> retrieveIncidentEdges(AdjacentEdgeMap& EMap,
                                                     unsigned vid){
    //no check for vid needed; EMap[vid] initializes if value not present (which is returned)
    std::unordered_set<Edge, Hash> Es = EMap[vid];
    EMap.erase(vid);

    for(auto& e : Es){
        unsigned uid = e.v.id == vid ? e.u.id : e.v.id;
         EMap[uid].erase(e);       
    }
    return Es;
}


//copied from http://stackoverflow.com/questions/236129/split-a-string-in-c
std::vector<std::string> &split(const std::string &s, char delim,
                                std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}



void makeRandSubset(const double p, std::map<unsigned, Vertex>& inMap, 
                    std::map<unsigned, Vertex>& outMap){
    // For randomization
    std::random_device rd;
    std::ranlux48_base e2(rd());
    std::uniform_real_distribution<> dist(0, 1);

    for(auto& ap : inMap){
        if(dist(e2) <= p){
            unsigned i = ap.first;
            outMap.emplace(std::piecewise_construct,
              std::forward_as_tuple(i),
              std::forward_as_tuple(ap.second, true)); //selfwitness=true
        }
    }
}



//Assumes A[0] is already set
void constructASets(const double p, const unsigned k, 
                    std::vector< std::map<unsigned, Vertex> >& A){
    // For randomization
    std::random_device rd;
    std::ranlux48_base e2(rd());
    std::uniform_real_distribution<> dist(0, 1);
    for(unsigned i=1; i<=k; ++i){
        makeRandSubset(p, A[i-1], A[i]);
    }
}



/*
    MIGHT BE NEEDED TO MAKE A-SETS as references/pointers on heap 
    (i.e. not stack, if stack is overflown).
    See http://stackoverflow.com/questions/10063037/unordered-set-storing-elements-as-pointers
*/
Vertex getZeroVertex(std::map<unsigned, Vertex>& A, 
                     std::vector<Edge>& zeroEdges){
    Vertex s = Vertex(UINT_MAX); //dummy id
    for(auto a : A){
        zeroEdges.push_back(Edge(s,a.second,0));
    }
    return s;
}




void getEdgeSMap(std::unordered_map<std::string, Edge>& E,
                 WitnessMap& witnessDistMap, 
                 std::unordered_map<std::string, Edge>& EOut){
    //Loop edges, update output set with edge if weight < dist_G(v, p_S(v)).
    for(auto& se : E){
        Edge e = se.second;
        
        // TODO make sure check if witnessDistMap.find() etc.
        // If nothing, then add edge (empty S set)
        if(witnessDistMap.find(e.u.id) == witnessDistMap.end()
           || witnessDistMap.find(e.v.id) == witnessDistMap.end()
           || e.weight < witnessDistMap[e.u.id].second 
           || e.weight < witnessDistMap[e.v.id].second){
            EOut[se.first] = e;

        }
    }
}


void printG(Graph& G){
    std::cout << "Printing G:" << std::endl;
    for(auto& vp : G.getV()){
        Vertex& v = vp.second;
        std::cout << "v.id = " << v.id << ", label = " << v.label << std::endl;
    }
    for(auto& ep : G.getE()){
        Edge& e = ep.second;
        std::cout << "e.u.id = " << e.u.id << ", e.v.id = " << e.v.id << ", weight = " << e.weight << std::endl;
    }
}




// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[10];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%X", &tstruct);

	return buf;
}



void updateMinEdge(std::unordered_map<unsigned, Edge>& m, unsigned key, 
                   Edge& e){
    if(m.find(key) != m.end()){
        if(m[key].weight > e.weight){
            m[key] = e;
        }
    }else{
        m[key] = e;
    }
}