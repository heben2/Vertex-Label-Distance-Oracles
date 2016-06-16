#include <iostream>
#include <cassert>
#include <vector>
#include <float.h>
// #include <time.h> 
#include <chrono>

#include "../oracles/lib/structs.hpp"
#include "../oracles/lib/sssp.hpp"
#include "../oracles/lib/std_methods.hpp"
#include "../oracles/lib/tst_methods.hpp"
#include "../oracles/lib/parser.hpp"
#include "../oracles/thorup_zwick_oracle.hpp"
#include "../oracles/chechik_oracle.hpp"
#include "../oracles/with_dijkstra_oracle.hpp"
#include "../oracles/with_thorup_zwick_oracle.hpp"
#include "../oracles/with_chechik_oracle.hpp"

#include "space_measurement.hpp"


namespace vertexLabelDistOras{
typedef std::vector<std::vector<double> > ResultList;

//MAYBE USE CHRONO INSTEAD?? http://en.cppreference.com/w/cpp/chrono/duration/duration_cast
// int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
// {
//     unsigned int resolution=1000000;
//     long int diff = (t2->tv_usec + resolution * t2->tv_sec) - (t1->tv_usec + resolution * t1->tv_sec);
//     result->tv_sec = diff / resolution;
//     result->tv_usec = diff % resolution;
//     return (diff<0);
// }


// wallclock
double runOracle(Oracle& O, Graph& G, unsigned const k, size_t& s){
    // {
    //     int tSize = 0, resident = 0, share = 0;
    //     std::ifstream buffer("/proc/self/statm");
    //     buffer >> tSize >> resident >> share;
    //     buffer.close();

    //     long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    //     double rss = resident * page_size_kb;
    //     std::cout << "RSS - " << rss << " kB\n";

    //     double shared_mem = share * page_size_kb;
    //     std::cout << "Shared Memory - " << shared_mem << " kB\n";

    //     std::cout << "Private Memory - " << rss - shared_mem << "kB\n";
    // }
    Graph GTmp(G); //make sure each preproc gets its own graph copy
    size_t s1 = getCurrentRSS();
    auto t1 = std::chrono::high_resolution_clock::now();
    int r = O.preproc(GTmp, k);
    auto t2 = std::chrono::high_resolution_clock::now();
    size_t s2 = getCurrentRSS();

    std::cout << "runOracle after preproc" << std::endl;
    // std::cout << "s2 = " << s2 << std::endl;
    // std::cout << "s2-s1 = " << s2-s1 << std::endl;
    s = s2 - s1;
    // s = 0;
    // {
    //     int tSize = 0, resident = 0, share = 0;
    //     std::ifstream buffer("/proc/self/statm");
    //     buffer >> tSize >> resident >> share;
    //     buffer.close();

    //     long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    //     double rss = resident * page_size_kb;
    //     std::cout << "RSS - " << rss << " kB\n";

    //     double shared_mem = share * page_size_kb;
    //     std::cout << "Shared Memory - " << shared_mem << " kB\n";

    //     std::cout << "Private Memory - " << rss - shared_mem << "kB\n";
    // }

    if(!r)
        return 0.;
    // auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
    
    return fp_ms.count();
}

//// cpu clock
// double runOracle(Oracle& O, Graph& G, unsigned const k){
//     std::clock_t c_start = std::clock();
//     O.preproc(G, k);
//     std::clock_t c_end = std::clock();

//     double fp_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
//     return fp_ms;
// }



double testWithDijkstra(Graph& G, unsigned const k, ResultList& results, 
                        size_t& s){
    WithDijkstraOracle wdo;// = ThorupZwickOracle<k>();
    // wdo.preproc(G, k);
    // size_t s = 0;
    double estimatedTime = runOracle(wdo, G, k, s);

    // thorup_zwick_oracle::thorup_zwick_oracle();
    //query all vertice pairs
    std::map<unsigned, Vertex>& V = G.getV();
    std::unordered_set<unsigned>& L = G.getL();
    for(auto& vp : V){
        for(auto l : L){
            results[vp.first][l] = wdo.query(vp.second, l);
            // std::cout << "main::query(" << vp.first << "," << l << ") = " << dist << std::endl;
        }
    }
    return estimatedTime;
}


double testWithThorupZwick(Graph& G, unsigned const k, ResultList& results,
                           size_t& s){
    WithThorupZwickOracle wtzo;// = ThorupZwickOracle<k>();
    // wtzo.preproc(G, k);
    double estimatedTime = runOracle(wtzo, G, k, s);
    // thorup_zwick_oracle::thorup_zwick_oracle();
    //query all vertice pairs
    std::map<unsigned, Vertex>& V = G.getV();
    std::unordered_set<unsigned>& L = G.getL();
    for(auto& vp : V){
        for(auto l : L){
            results[vp.first][l] = wtzo.query(vp.second, l);
            // std::cout << "main::query(" << vp.first << "," << l << ") = " << dist << std::endl;
        }
    }

    return estimatedTime;
}

double testWithChechik(Graph& G, unsigned const k, ResultList& results, 
                       size_t& s){
    WithChechikOracle wco;// = ThorupZwickOracle<k>();
    // wtzo.preproc(G, k);
    double estimatedTime = runOracle(wco, G, k, s);
    // thorup_zwick_oracle::thorup_zwick_oracle();
    //query all vertice pairs
    if(estimatedTime > 0){
        std::map<unsigned, Vertex>& V = G.getV();
        std::unordered_set<unsigned>& L = G.getL();
        for(auto& vp : V){
            for(auto l : L){
                results[vp.first][l] = wco.query(vp.second, l);
                // std::cout << "main::query(" << vp.first << "," << l << ") = " << dist << std::endl;
            }
        }
    }
    return estimatedTime;
}




double testChechik(Graph& G, unsigned const k, ResultList& results, size_t& s){
std::cout << "testChechik preproc begin" << std::endl;
    ChechikOracle co;// = ThorupZwickOracle<k>();
    // co.preproc(G, k);
    double estimatedTime = runOracle(co, G, k, s);

std::cout << "testChechik preproc done" << std::endl;
    // thorup_zwick_oracle::thorup_zwick_oracle();
    //query all vertice pairs
    std::map<unsigned, Vertex>& V = G.getV();
    std::unordered_set<unsigned>& L = G.getL();
    for(auto& vp : V){
        for(auto& l : L){
            results[vp.first][l] = co.query(vp.second, l);
            // std::cout << "Chechik query(" << vp.first << "," << l << ") = " << dist << std::endl;
        }
    }
std::cout << "testChechik done" << std::endl;
    return estimatedTime;
}




//compute the actual shortest distances.
void getShortestLabelDistances(Graph& G, ResultList& resultsTrue){
// std::cout << "getShortestLabelDistances begin" << std::endl;
    std::map<unsigned, Vertex>& V = G.getV();
    AdjacentEdgeMap* adjEMap = new AdjacentEdgeMap();
    getIncidentEdgesMap(G.getE(), *adjEMap);
    for(auto& vp : V){
        AdjacentEdgeMap* adjEMapTmp = new AdjacentEdgeMap(*adjEMap);
        sssp::dijkstraSourceLabelsDists(G, *adjEMapTmp, vp.second, 
                                        resultsTrue[vp.first]);
        delete adjEMapTmp;
        // for(auto d : resultsTrue[vp.first] )
            // std::cout << "getShortestLabelDistances v.id = " << vp.first << " with dist = " << d << std::endl;

    }
    delete adjEMap;
// std::cout << "getShortestLabelDistances end" << std::endl;
}


//TODO run through results and compare
void verifyStretchBound(ResultList& resultsTrue, ResultList& resultsOther, 
                        unsigned const stretch, unsigned const labelSize,
                        Graph& G){
// std::cout << "verifyStretchBound begin, stretch = " << stretch << std::endl;
    unsigned const n = resultsTrue.size();
// std::cout << "verifyStretchBound begin, n = " << n << std::endl;
    for(unsigned i=0; i < n; ++i){
        for(unsigned l=0; l < labelSize; ++l){
            double trueResult = resultsTrue[i][l];
            double otherResult = resultsOther[i][l];
            if(!(trueResult <= otherResult && otherResult <= stretch*trueResult)){
                std::cout << "i = " << i << ", l = " << l << ", trueResult = " << trueResult << ", otherResult = " << otherResult << std::endl;
                std::cout << "true label = " << G.getV()[i].label << std::endl;
            }
            // std::cout << "i = " << i << ", l = " << l << ", trueResult = " << trueResult << ", otherResult = " << otherResult << std::endl;
            assert(trueResult <= otherResult);
            assert(otherResult <= stretch*trueResult);
        }
    }
// std::cout << "verifyStretchBound end" << std::endl;
}


void initResultsVector(ResultList* results, unsigned const n, unsigned const l){
// std::cout << "initResultsVector begin" << std::endl;
    results->resize(n);
    for(auto& r : *results){
        r.resize(l,DBL_MAX);
    }
}


//Right now with fixed l-size
// void testVertexLabelDistOras(unsigned minN, unsigned maxN, unsigned stepSizeN,
//                              unsigned minK, unsigned maxK, unsigned stepSizeK,
//                              unsigned const l){
void testStretchCorrectness(TestParams tstParams, std::string outPath){
    std::cout << currentDateTime() << " " << "Testing STRETCH for Vertex-labeled distance oracles with l = " << tstParams.l << std::endl;
    std::string path = "tmp.txt";
    // std::string outPath = "timeResults.json";
    unsigned numOracles = 4;
    size_t s;

    outputParams(outPath, tstParams);

    std::vector< std::vector< double> > timeMeasurements;
    for(unsigned i = 0; i < numOracles; ++i)
        timeMeasurements.push_back(std::vector< double>());

    for(unsigned n = tstParams.minN; n <= tstParams.maxN; 
        n += tstParams.stepSizeN){
        for(unsigned k = tstParams.minK; k <= tstParams.maxK; 
            k += tstParams.stepSizeK){
            unsigned stretch = 4*k-5;
            // parser::createFullLabelGraph(path, n, tstParams.l);
            // Graph* G = parser::loadGraph(path);
            Graph* G = parser::returnFullLabelGraph(path, n, tstParams.l);
            // printG(*G);
    std::cout << currentDateTime() << " " << "Running with n = " << n << " and k = " << k << std::endl;

            for(unsigned i=0; i < tstParams.repeatIter; ++i){
                ResultList* resultsTrue = new ResultList(); 
                initResultsVector(resultsTrue, n, tstParams.l);
                getShortestLabelDistances(*G, *resultsTrue);
    std::cout << currentDateTime() << " " << "True results computed, now testing Chechik's oracle." << std::endl;


                ResultList* resultsChechik = new ResultList();
                initResultsVector(resultsChechik, n, tstParams.l);
                timeMeasurements[0].push_back(testChechik(*G, k, *resultsChechik, s));
                verifyStretchBound(*resultsTrue, *resultsChechik, stretch, 
                                   tstParams.l, *G);
    std::cout << currentDateTime() << " " << "Chechik's oracle verified, now testing With Dijkstra oracle." << std::endl;


                ResultList* resultsWithDijkstra = new ResultList();
                initResultsVector(resultsWithDijkstra, n, tstParams.l);
                timeMeasurements[1].push_back(testWithDijkstra(*G, k, *resultsWithDijkstra, s));
                verifyStretchBound(*resultsTrue, *resultsWithDijkstra, stretch, tstParams.l, *G);
    std::cout << currentDateTime() << " " << "With Dijkstra's oracle verified, now testing With Thorup-Zwick oracle." << std::endl;

                ResultList* resultsWithTZ = new ResultList();
                if(k >= 3){
                    initResultsVector(resultsWithTZ, n, tstParams.l);
                    timeMeasurements[2].push_back(testWithThorupZwick(*G, k, *resultsWithTZ, s));
                    verifyStretchBound(*resultsTrue, *resultsWithTZ, stretch, tstParams.l, *G);
    std::cout << currentDateTime() << " " << "With Thorup-Zwick's oracle verified" << std::endl;
                } else{
                    timeMeasurements[2].push_back(0);
                }

                ResultList* resultsWithC = new ResultList();
                if(k > 3){
                    initResultsVector(resultsWithC, n, tstParams.l);
                    double r = testWithChechik(*G, k, *resultsWithC, s);
                    timeMeasurements[3].push_back(r);
                    if(r > 0){
                        std::cout << currentDateTime() << " " << "Verifying With Chechik's oracle" << std::endl;
                        verifyStretchBound(*resultsTrue, *resultsWithC, stretch, tstParams.l, *G);
                        std::cout << currentDateTime() << " " << "With Chechik's oracle verified" << std::endl;
                    }
                } else{
                    timeMeasurements[3].push_back(0);
                }

    std::cout << currentDateTime() << " " << "Tests passed for n = " << n << ", k = " << k << ", repetition iteration " << i<< std::endl;
            
                delete resultsTrue;
                delete resultsChechik;
                delete resultsWithDijkstra;
                delete resultsWithTZ;
                delete resultsWithC;
                outputResults<double>(outPath, timeMeasurements);
            }
            delete G;
        }
    }
    // printResults(outPath, timeMeasurements, tstParams);
    outputResults<double>(outPath, timeMeasurements);
}


//Compute (mean,var) for stretch over whole run. Return resultsStretch
std::pair<double, double> computeStretchMean(ResultList& resultsTrue, 
                                             ResultList& resultsOther, 
                                             unsigned const stretch, 
                                             unsigned const labelSize,
                                             Graph& G){
    double mean = 0,
           var = 0;
    std::vector<double> stretches;
    unsigned const n = resultsTrue.size();
    for(unsigned i=0; i < n; ++i){
        for(unsigned l=0; l < labelSize; ++l){
            double trueResult = resultsTrue[i][l];
            double otherResult = resultsOther[i][l];
            if(trueResult != 0){
                double d = otherResult/trueResult;
                stretches.push_back(d);
                mean += d;
            }
        }
    }
    mean = mean/stretches.size();
    for(auto d : stretches)
        var += pow(d - mean, 2);
    var = var/stretches.size();
    return {mean, var};
// std::cout << "verifyStretchBound end" << std::endl;
}


//Testing stretch accuracy
void testStretch(TestParams tstParams, std::string outPath){
    std::cout << currentDateTime() << " " << "Testing STRETCH ACCURACY for Vertex-labeled distance oracles with l = " << tstParams.l << std::endl;
    std::string path = "tmp.txt";
    // std::string outPath = "timeResults.json";
    unsigned numOracles = 4;
    size_t s;

    outputParams(outPath, tstParams);

    std::vector< std::vector< std::pair<double, double> > > stretchResults;
    for(unsigned i = 0; i < numOracles; ++i)
        stretchResults.push_back(std::vector< std::pair<double, double> >());

    for(unsigned l = tstParams.minL; l <= tstParams.maxL; 
        l += tstParams.stepSizeL){
        for(unsigned k = tstParams.minK; k <= tstParams.maxK; 
            k += tstParams.stepSizeK){
            unsigned stretch = 4*k-5;
            // parser::createFullLabelGraph(path, tstParams.n, l);
            // Graph* G = parser::loadGraph(path);
            Graph* G = parser::returnFullLabelGraph(path, tstParams.n, l);
            // printG(*G);
    std::cout << currentDateTime() << " " << "Running with tstParams.n = " << tstParams.n << " and k = " << k << std::endl;

            for(unsigned i=0; i < tstParams.repeatIter; ++i){
                ResultList* resultsTrue = new ResultList(); 
                initResultsVector(resultsTrue, tstParams.n, l);
                getShortestLabelDistances(*G, *resultsTrue);
    std::cout << currentDateTime() << " " << "True results computed, now testing Chechik's oracle." << std::endl;


                ResultList* resultsChechik = new ResultList();
                initResultsVector(resultsChechik, tstParams.n, l);
                testChechik(*G, k, *resultsChechik, s);
                stretchResults[0].push_back(computeStretchMean(*resultsTrue, 
                                                               *resultsChechik, 
                                                               stretch, 
                                                               l, 
                                                               *G));
    std::cout << currentDateTime() << " " << "Chechik's oracle verified, now testing With Dijkstra oracle." << std::endl;


                ResultList* resultsWithDijkstra = new ResultList();
                initResultsVector(resultsWithDijkstra, tstParams.n, l);
                testWithDijkstra(*G, k, *resultsWithDijkstra, s);
                stretchResults[1].push_back(computeStretchMean(*resultsTrue, *resultsWithDijkstra, stretch, l, *G));
    std::cout << currentDateTime() << " " << "With Dijkstra's oracle verified, now testing With Thorup-Zwick oracle." << std::endl;

                ResultList* resultsWithTZ = new ResultList();
                if(k >= 3){
                    initResultsVector(resultsWithTZ, tstParams.n, l);
                    testWithThorupZwick(*G, k, *resultsWithTZ, s);
                    stretchResults[2].push_back(computeStretchMean(*resultsTrue, *resultsWithTZ, stretch, l, *G));
                    
    std::cout << currentDateTime() << " " << "With Thorup-Zwick's oracle verified" << std::endl;
                } else{
                    stretchResults[2].push_back({0,0});
                }

                ResultList* resultsWithC = new ResultList();
                if(k > 3){
                    initResultsVector(resultsWithC, tstParams.n, l);
                    double r = testWithChechik(*G, k, *resultsWithC, s);
                    if(r > 0){
                        stretchResults[3].push_back(computeStretchMean(*resultsTrue, *resultsWithC, stretch, l, *G));
                    } else{
                        stretchResults[3].push_back({0,0});
                    }
                } else{
                    stretchResults[3].push_back({0,0});
                }

    std::cout << currentDateTime() << " " << "Tests passed for tstParams.n = " << tstParams.n << ", k = " << k << ", repetition iteration " << i<< std::endl;
            
                delete resultsTrue;
                delete resultsChechik;
                delete resultsWithDijkstra;
                delete resultsWithTZ;
                delete resultsWithC;
                outputResultPairs(outPath, stretchResults);
            }
            delete G;
        }
    }
}



void testPreprocTime(TestParams tstParams, std::string outPath, 
                     std::string outPathSize){
    std::cout << "Testing Vertex-labeled distance oracles with l = " << tstParams.l << std::endl;
    std::string path = "tmp.txt";
    unsigned numOracles = 4;
    size_t s;
    outputParams(outPath, tstParams);
    outputParams(outPathSize, tstParams);
    
    std::vector< std::vector< double> > timeMeasurements;
    for(unsigned i = 0; i < numOracles; ++i)
        timeMeasurements.push_back(std::vector< double>());

    std::vector< std::vector< size_t> > sizeMeasurements;
    for(unsigned i = 0; i < numOracles; ++i)
        sizeMeasurements.push_back(std::vector< size_t>());

    for(unsigned n = tstParams.minN; n <= tstParams.maxN; 
        n += tstParams.stepSizeN){
        for(unsigned k = tstParams.minK; k <= tstParams.maxK; 
            k += tstParams.stepSizeK){
            unsigned stretch = 4*k-5;
			std::cout << currentDateTime() << " " << "Create full labeled graph" << std::endl;
            // parser::createFullLabelGraph(path, n, tstParams.l);
			// std::cout << "Load graph" << std::endl;
            // Graph* G = parser::loadGraph(path);
			Graph* G = parser::returnFullLabelGraph(path, n, tstParams.l);
            // printG(*G);
    std::cout << currentDateTime() << " " << "Running with n = " << n << " and k = " << k << std::endl;

            for(unsigned i=0; i < tstParams.repeatIter; ++i){
    std::cout << currentDateTime() << " " << "Testing Chechik's oracle." << std::endl;
                ChechikOracle co;
                timeMeasurements[0].push_back(runOracle(co, *G, k, s));
                sizeMeasurements[0].push_back(s);
    std::cout << currentDateTime() << " " << "Chechik's oracle done, now testing With Dijkstra oracle." << std::endl;

                WithDijkstraOracle wdo;
                timeMeasurements[1].push_back(runOracle(wdo, *G, k, s));
                sizeMeasurements[1].push_back(s);
    std::cout << currentDateTime() << " " << "With Dijkstra's oracle done, now testing With Thorup-Zwick oracle." << std::endl;

                if(k >= 3){
                    WithThorupZwickOracle wtzo;
                    timeMeasurements[2].push_back(runOracle(wtzo, *G, k, s));
                    sizeMeasurements[2].push_back(s);
                    std::cout << currentDateTime() << " " << "With Thorup-Zwick's oracle done, now testing With Chechik's." << std::endl;
                    WithChechikOracle wco;
                    timeMeasurements[3].push_back(runOracle(wco, *G, k, s));
                    sizeMeasurements[3].push_back(s);
    std::cout << currentDateTime() << " " << "With Chechik's oracle done." << std::endl;
                } else{
                    timeMeasurements[2].push_back(0);
                    timeMeasurements[3].push_back(0);

                    sizeMeasurements[2].push_back(0);
                    sizeMeasurements[3].push_back(0);
                }
    //std::cout << currentDateTime() << " " << "Tests passed for n = " << n << ", k = " << k << ", repetition iteration " << i<< std::endl;
                outputResults<double>(outPath, timeMeasurements);
                outputResults<size_t>(outPathSize, sizeMeasurements);
            }
            delete G;
        }
    }
    // printResults(outPath, timeMeasurements, tstParams);
    // outputResults<double>(outPath, timeMeasurements);
}



void testPreprocTimeRealGraph(std::string graphPath, unsigned k, unsigned numLabels
                     , std::string outPath){
    unsigned numOracles = 4;
    size_t s; //dummy for size
    double t = 0;
    
    std::vector< std::vector< double> > timeMeasurements;
    for(unsigned i = 0; i < numOracles; ++i)
        timeMeasurements.push_back(std::vector< double>());

    Graph* G = new Graph();
 std::cout << currentDateTime() << " " << "Loading Graph. Might take some time." << std::endl;
    parser::loadGraphFromFileGr(graphPath, *G, numLabels);
 std::cout << currentDateTime() << " " << "Graph loaded; testing Chechik's oracle." << std::endl;

    {
        WithDijkstraOracle wdo;
        t = runOracle(wdo, *G, k, s);
        timeMeasurements[1].push_back(t);
std::cout << currentDateTime() << " " << "With Dijkstra's oracle done in "<< t << " ms! Now testing With Thorup-Zwick oracle." << std::endl;
    }
    outputResults<double>(outPath, timeMeasurements);
    {
        WithThorupZwickOracle wtzo;
        t = runOracle(wtzo, *G, k, s);
        timeMeasurements[2].push_back(t);
std::cout << currentDateTime() << " " << "With Thorup-Zwick's oracle done in "<< t << " ms! Now testing With Chechik's." << std::endl;
    }
    outputResults<double>(outPath, timeMeasurements);
    {
        WithChechikOracle wco;
        t = runOracle(wco, *G, k, s);
        timeMeasurements[3].push_back(t);
std::cout << currentDateTime() << " " << "With Chechik's oracle done in "<< t << " ms! All done now." << std::endl;

    }
    outputResults<double>(outPath, timeMeasurements);
    {
        ChechikOracle co;
        t = runOracle(co, *G, k, s);
        timeMeasurements[0].push_back(t);
std::cout << currentDateTime() << " " << "Chechik's oracle done in "<< t << " ms! Now testing With Dijkstra oracle." << std::endl;
    }
    outputResults<double>(outPath, timeMeasurements);
}









void testPreprocSize(TestParams tstParams, std::string outPathSize){
    std::cout << "Testing SIZE of Vertex-labeled distance oracles with l = " << tstParams.l << std::endl;
    std::string path = "tmp.txt";
    unsigned numOracles = 4;
    size_t s; //dummy

    outputParams(outPathSize, tstParams);
    
    std::vector< std::vector< size_t> > sizeMeasurements;
    for(unsigned i = 0; i < numOracles; ++i)
        sizeMeasurements.push_back(std::vector< size_t>());

    for(unsigned n = tstParams.minN; n <= tstParams.maxN; 
        n += tstParams.stepSizeN){
        for(unsigned k = tstParams.minK; k <= tstParams.maxK; 
            k += tstParams.stepSizeK){
            std::cout << currentDateTime() << " " << "Create full labeled graph" << std::endl;
            Graph* G = parser::returnFullLabelGraph(path, n, tstParams.l);
    std::cout << currentDateTime() << " " << "Running with n = " << n << " and k = " << k << std::endl;

            for(unsigned i=0; i < tstParams.repeatIter; ++i){
    std::cout << currentDateTime() << " " << "Testing Chechik's oracle." << std::endl;
                ChechikOracle co;
                runOracle(co, *G, k, s);
                sizeMeasurements[0].push_back(co.getSize());
    std::cout << currentDateTime() << " " << "Chechik's oracle done, now testing With Dijkstra oracle." << std::endl;

                WithDijkstraOracle wdo;
                runOracle(wdo, *G, k, s);
                sizeMeasurements[1].push_back(wdo.getSize());
    std::cout << currentDateTime() << " " << "With Dijkstra's oracle done, now testing With Thorup-Zwick oracle." << std::endl;

                if(k >= 3){
                    WithThorupZwickOracle wtzo;
                    runOracle(wtzo, *G, k, s);
                    sizeMeasurements[2].push_back(wtzo.getSize());
                    std::cout << currentDateTime() << " " << "With Thorup-Zwick's oracle done, now testing With Chechik's." << std::endl;
                    WithChechikOracle wco;
                    runOracle(wco, *G, k, s);
                    sizeMeasurements[3].push_back(wco.getSize());
    std::cout << currentDateTime() << " " << "With Chechik's oracle done." << std::endl;
                } else{
                    sizeMeasurements[2].push_back(0);
                    sizeMeasurements[3].push_back(0);
                }
    //std::cout << currentDateTime() << " " << "Tests passed for n = " << n << ", k = " << k << ", repetition iteration " << i<< std::endl;
                outputResults<size_t>(outPathSize, sizeMeasurements);
            }
            delete G;
        }
    }
    // printResults(outPath, timeMeasurements, tstParams);
    // outputResults<double>(outPath, timeMeasurements);
}



}