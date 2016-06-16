#include <iostream>
#include <cassert>

#include "tst_methods.hpp"
#include "structs.hpp"
#include "sssp.hpp"
#include "std_methods.hpp"

// #ifndef TST_METHODS
// #define TST_METHODS

TestParams::TestParams(unsigned minN, unsigned maxN, unsigned stepSizeN,
                       unsigned minK, unsigned maxK, unsigned stepSizeK,
                       unsigned const l, unsigned repeatIter):
minN(minN), maxN(maxN), stepSizeN(stepSizeN), minK(minK), maxK(maxK), 
stepSizeK(stepSizeK), l(l), repeatIter(repeatIter) {}

// TestParams::TestParams(unsigned minL, unsigned maxL, unsigned minK, 
//                        unsigned maxK, unsigned stepSizeK, unsigned stepSizeL,
//                        unsigned const n, unsigned repeatIter):
// minK(minK), maxK(maxK), stepSizeK(stepSizeK), minL(minL), maxL(maxL), 
// stepSizeL(stepSizeL), n(n), repeatIter(repeatIter) {}

TestParams::TestParams(unsigned minN, unsigned maxN, unsigned stepSizeN,
                       unsigned minK, unsigned maxK, unsigned stepSizeK,
                       unsigned minL, unsigned maxL, unsigned stepSizeL,
                       unsigned const n, unsigned repeatIter):
minN(minN), maxN(maxN), stepSizeN(stepSizeN), minK(minK), maxK(maxK), 
stepSizeK(stepSizeK), minL(minL), maxL(maxL), stepSizeL(stepSizeL), n(n),
repeatIter(repeatIter) {}


// TestParams::TestParams(unsigned minN, unsigned maxN, unsigned stepSizeN,
//                        unsigned minK, unsigned maxK, unsigned stepSizeK,
//                        unsigned const l, unsigned repeatIter):
// minN(minN), maxN(maxN), stepSizeN(stepSizeN), minK(minK), maxK(maxK), 
// stepSizeK(stepSizeK), l(l), repeatIter(repeatIter) {}



//compute the actual shortest distances.
void getShortestDistances(Graph& G, ResultList& resultsTrue){
    std::map<unsigned, Vertex>& V = G.getV();

    AdjacentEdgeMap* adjEMap = new AdjacentEdgeMap();
    getIncidentEdgesMap(G.getE(),*adjEMap);
    for(auto& vp : V){
        // Graph GTmp = Graph(G);//TODO remove cpy?
        AdjacentEdgeMap* adjEMapTmp = new AdjacentEdgeMap(*adjEMap);
        WitnessMap* wmap = sssp::dijkstra(G, *adjEMapTmp, vp.second);
        delete adjEMapTmp;
        for(auto& wp : *wmap){
            resultsTrue[vp.first][wp.first] = wp.second.second;
        }
        delete wmap;
    }
    delete adjEMap;
}


//TODO run through results and compare
void verifyStretchBound(ResultList& resultsTrue, ResultList& resultsOther, 
                        const unsigned stretch, bool print){
std::cout << "verifyStretchBound begin, stretch = " << stretch << std::endl;
    unsigned const n = resultsTrue.size();
std::cout << "verifyStretchBound begin, n = " << n << std::endl;
    for(unsigned i=0; i < n; ++i){
        for(unsigned j=0; j < n; ++j){
            double trueResult = resultsTrue[i][j];
            double otherResult = resultsOther[i][j];
            if(print){
                std::cout << "i = " << i << ", j = " << j << ", trueResult = " << trueResult << ", otherResult = " << otherResult << std::endl;
            }
            assert(trueResult <= otherResult);
            assert(otherResult <= stretch*trueResult);
        }
    }
// std::cout << "verifyStretchBound end" << std::endl;
}


void initResultsVector(ResultList* results, unsigned const n){
// std::cout << "initResultsVector begin" << std::endl;
    results->resize(n);
    for(auto& l : *results){
        l.resize(n,-1);
    }
}

void outputParams(std::string path, TestParams tp){
    std::ofstream f;
    f.open(path, std::fstream::out | std::fstream::app);

    f << "minK : " << tp.minK << std::endl;
    f << "maxK : " << tp.maxK << std::endl;
    f << "stepSizeK : " << tp.stepSizeK << std::endl;
    if(tp.n > 0){//iterating for fixed n
        f << "minL : " << tp.minL << std::endl;
        f << "maxL : " << tp.maxL << std::endl;
        f << "stepSizeL : " << tp.stepSizeL << std::endl;
        f << "n : " << tp.n << std::endl;
    }
    if(tp.l > 0){ //for fixed l
        f << "minN : " << tp.minN << std::endl;
        f << "maxN : " << tp.maxN << std::endl;
        f << "stepSizeN : " << tp.stepSizeN << std::endl;
        f << "l : " << tp.l << std::endl;
    }
    f << "repeatIter : " << tp.repeatIter << std::endl;

    f.close();
}


void printResults(std::string path, 
                  std::vector< std::vector< double> > results,
                  TestParams tstParams){
    std::ofstream f;
    f.open(path);
    for(unsigned i = 0; i < results.size(); ++i){
        f << i << std::endl;
        for(auto& d : results[i])
            f << d << ",";
        f << std::endl;
    }
    f.close();
}



void outputResultPairs(std::string path,
                std::vector< std::vector< std::pair<double, double>> >& list){
    std::ofstream f;
    f.open(path, std::fstream::out | std::fstream::app);

    std::string s = "";
    f << "[";
    for(auto a : list){
        s += "[";
        for(auto b : a){
            s += "(" + std::to_string(b.first) + ", " 
                     + std::to_string(b.second) + "), ";
        }
        s = s.substr(0, s.size()-2);
        s += "], \n";
    }
    s = s.substr(0, s.size()-3);
    f << s + "]" << std::endl;
    f.close();
}
