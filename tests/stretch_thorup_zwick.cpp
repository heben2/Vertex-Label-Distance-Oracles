#include <iostream>
#include <cassert>
#include <vector>
#include <random>
#include <ctime>
#include <chrono>


#include "oracles/lib/structs.hpp"
#include "oracles/thorup_zwick_oracle.hpp"
#include "oracles/lib/sssp.hpp"
#include "oracles/lib/std_methods.hpp"
#include "oracles/lib/tst_methods.hpp"
#include "oracles/spanner.hpp"


/*
Test correctness of the Thorup-Zwick oracle and spanner.
Used for the other oracles, so important that it is correct.
*/


void testThorupZwickOra(Graph& G, unsigned const k, ResultList& results){
// std::cout << currentDateTime() << " " << "testThorupZwickOra begin" << std::endl;

    ThorupZwickOracle tzo;// = ThorupZwickOracle<k>();
    tzo.preproc(G, k);
    // thorup_zwick_oracle::thorup_zwick_oracle();
    //query all vertice pairs
    std::map<unsigned, Vertex>& V =  G.getV();
    // double dist = tzo.query(0, 1);
    for(auto& vp : V){
        for(auto& up : V){
            // if(vp.first != up.first){
                results[vp.first][up.first] = tzo.query(vp.second, up.second);
                // double dist = tzo.query(vp.second, up.second);
                // std::cout << currentDateTime() << " " << "main::query(" << vp.first << "," << up.first << ") = " << dist << std::endl;
        }
    }
// std::cout << currentDateTime() << " " << "testThorupZwickOra end" << std::endl;
}

//How to test? Same rules apply, just use correctness method on spanner S.
void testThorupZwickSpanner(Graph& G, unsigned const k, Graph& S){
// std::cout << currentDateTime() << " " << "testThorupZwickSpanner begin" << std::endl;

    ThorupZwickOracle tzo;// = ThorupZwickOracle<k>();
    tzo.constructSpanner(G, S, k);
    // printG(*S);
// std::cout << currentDateTime() << " " << "testThorupZwickSpanner end" << std::endl;
}




//Test both the Thorup-Zwick oracle and spanner for stretch bound.
void testThorupZwick(unsigned minN, unsigned maxN, unsigned stepSizeN,
                     unsigned minK, unsigned maxK, unsigned stepSizeK){
    std::cout << currentDateTime() << " " << "Testing Thorup-Zwick oracle stretch bound" << std::endl;
    std::string path = "tmp.txt";

    for(unsigned n = minN; n < maxN; n+=stepSizeN){
        for(unsigned k = minK; k < maxK; k+=stepSizeK){
std::cout << currentDateTime() << " " << "Running with n = " << n << " and k = " << k << std::endl;
            unsigned stretch = 2*k-1;
            parser::createFullGraph(path, n);
            Graph* G = parser::loadGraph(path);
            // printG(*G);
            ResultList* resultsTrue = new ResultList(); 
            initResultsVector(resultsTrue, n);
            getShortestDistances(*G, *resultsTrue);
            Graph* S = new Graph();
            testThorupZwickSpanner(*G, k, *S);
            ResultList* resultsSpanner = new ResultList();
            initResultsVector(resultsSpanner, n);
            getShortestDistances(*S, *resultsSpanner);
            verifyStretchBound(*resultsTrue, *resultsSpanner, stretch, false);

std::cout << currentDateTime() << " " << "Thorup-Zwick spanner verified, now testing Thorup-Zwick oracle." << std::endl;
            ResultList* resultsTZO = new ResultList();
            initResultsVector(resultsTZO, n);
            testThorupZwickOra(*G, k, *resultsTZO);
            verifyStretchBound(*resultsTrue, *resultsTZO, stretch, false);
std::cout << currentDateTime() << " " << "Tests passed for n = " << n << ", k = " << k << std::endl;
        
            delete G;
            delete S;
            delete resultsTrue;
            delete resultsSpanner;
            delete resultsTZO;
        }
    }
}





void testSpannerTime(TestParams tstParams){
    std::string path = "spanner.graph.txt",
                outPath = "spanner.timeResults.json";
    unsigned l = 1;
    //unsigned const k = 3; //> 1
    //unsigned stretch = 2*k-1;
    // parser::createFullGraph(path, n);
    // Graph* G = parser::loadGraph(path);
    // Graph* G =  parser::returnFullLabelGraph(path, n, l);
    // Graph* G = parser::loadGraph(path); //reuse graph for debugging
    Graph spanner_bs,
          spanner_tz;

    parser::clearFile(outPath);
    outputParams(outPath, tstParams);
    // unsigned //for short tests
    //     minN = 1000,
    //     maxN = 10000,
    //     stepSizeN = 1000,
    //     minK = 2,
    //     maxK = 12,
    //     stepSizeK = 1;

    std::chrono::duration<double, std::milli> fp_ms;
    
    unsigned numSpanners = 2;
    std::vector< std::vector< double> > timeMeasurements;
    for(unsigned i = 0; i < numSpanners; ++i)
        timeMeasurements.push_back(std::vector< double>());
    unsigned n = tstParams.minN;
    unsigned mfactor = 4;
    unsigned minM = floor((n*(n - 1) / 2) / mfactor),
        maxM = floor(n*(n - 1) / 2),
        stepSizeM = minM;
    /*for(unsigned n = tstParams.minN; n <= tstParams.maxN; 
        n += tstParams.stepSizeN){*/
    for(unsigned m = minM; m <= maxM; m += stepSizeM){
        for(unsigned k = tstParams.minK; k <= tstParams.maxK; 
            k += tstParams.stepSizeK){
            /*unsigned maxM = floor(n*(n-1)/2);
            unsigned minM = floor((n*(n-1)/2)/2);
            int m = minM + rand() % (maxM - minM +1);*/

            std::cout << currentDateTime() << " " << "Running with n = " << n << ", m = " << m << " and k = " << k << std::endl;
            unsigned stretch = 2*k-1;
            // Graph* G = parser::returnFullLabelGraph(path, n, l);
            
            // Graph* G = parser::loadGraph(path); //reuse graph for debugging
            // printG(*G);


            for(unsigned i=0; i < tstParams.repeatIter; ++i){
                Graph* G = parser::returnLabelGraph(path, n, m, l);
                std::cout << currentDateTime() << " " << "G.E.size() = " << G->getESize() << std::endl;

                std::cout << currentDateTime() << " " << "Constructing Baswanna Sen spanner\n";
                auto t1 = std::chrono::high_resolution_clock::now();
                spanner::constructSpanner(*G, spanner_bs, k);
                auto t2 = std::chrono::high_resolution_clock::now();
                std::cout << currentDateTime() << " " << "G.getE().size() = " << G->getE().size() << ", spanner_bs.getE().size() = " << spanner_bs.getE().size() << std::endl;
                fp_ms = t2 - t1;
                timeMeasurements[0].push_back(fp_ms.count());

                std::cout << currentDateTime() << " " << "Spanner_bs done in " << fp_ms.count() << "ms!\nNow Constructing Thorup-Zwick spanner\n";
                ThorupZwickOracle tzo;// = ThorupZwickOracle<k>();
                t1 = std::chrono::high_resolution_clock::now();
                tzo.constructSpanner(*G, spanner_tz, k);
                t2 = std::chrono::high_resolution_clock::now();
                std::cout << currentDateTime() << " " << "G.getE().size() = " << G->getE().size() << ", spanner_tz.getE().size() = " << spanner_tz.getE().size() << std::endl;
                fp_ms = t2 - t1;
                timeMeasurements[1].push_back(fp_ms.count());

                std::cout << currentDateTime() << " " << "Spanner_tz done in " << fp_ms.count() << "ms!\nNow verifying spanners." << std::endl;

                // printG(spanner);

                std::cout << currentDateTime() << " " << "Getting actual shortest distances on G." << std::endl;
                ResultList* resultsTrue = new ResultList(); 
                initResultsVector(resultsTrue, n);
                getShortestDistances(*G, *resultsTrue);

                ResultList* resultsSpanner_bs = new ResultList();
                ResultList* resultsSpanner_tz = new ResultList();

                std::cout << currentDateTime() << " " << "Verifying Thorup-Zwick spanner" << std::endl;
                initResultsVector(resultsSpanner_bs, n);
                getShortestDistances(spanner_bs, *resultsSpanner_bs);
                verifyStretchBound(*resultsTrue, *resultsSpanner_bs, stretch, false);
                
                std::cout << currentDateTime() << " " << "Verifying Thorup-Zwick spanner" << std::endl;
                initResultsVector(resultsSpanner_tz, n);
                getShortestDistances(spanner_bs, *resultsSpanner_tz);
                verifyStretchBound(*resultsTrue, *resultsSpanner_tz, stretch, false);

                outputResults<double>(outPath, timeMeasurements);
                delete resultsTrue;
                delete resultsSpanner_bs;
                delete resultsSpanner_tz;
                delete G;
            }
            
        }
    }
}



void testSpannerSize(unsigned n, unsigned minK, unsigned maxK, 
                     unsigned stepSizeK, unsigned repeatIter){
    std::cout << "Testing spanner sizes" << std::endl;
    std::string path = "spanner.graph.txt",
                outPath = "spanner.sizeResults.json";
    unsigned l = 1;

    parser::clearFile(outPath);
    // outputParams(outPath, tstParams);
    // unsigned //for short tests
    //     minN = 1000,
    //     maxN = 10000,
    //     stepSizeN = 1000,
    //     minK = 2,
    //     maxK = 12,
    //     stepSizeK = 1;

    unsigned numSpanners = 2;
    std::vector< std::vector< unsigned> > sizeMeasurements;
    for(unsigned i = 0; i < numSpanners; ++i)
        sizeMeasurements.push_back(std::vector< unsigned>());

    for(unsigned k = minK; k <= maxK; k += stepSizeK){
        outputParam<unsigned>(outPath, "k", k);
        outputParam<unsigned>(outPath, "n", n);

        std::cout << currentDateTime() << " " << "Running with n = " << n << " and k = " << k << std::endl;
        // Graph* G = parser::returnLabelGraph(path, n, m, l);

        for(unsigned i=0; i < repeatIter; ++i){
            Graph spanner_bs,
                  spanner_tz;
            Graph* G = parser::returnFullLabelGraph(path, n, l);
            std::cout << currentDateTime() << " " << "Constructing Baswanna Sen spanner" << std::endl;;
            spanner::constructSpanner(*G, spanner_bs, k);
            std::cout << currentDateTime() << " " << "BS done! G.getE().size() = " << G->getE().size() << ", spanner_bs.getE().size() = " << spanner_bs.getE().size() << std::endl;
            sizeMeasurements[0].push_back(spanner_bs.getESize());

            std::cout << currentDateTime() << " " << "Constructing Thorup-Zwick spanner" << std::endl;
            ThorupZwickOracle tzo;// = ThorupZwickOracle<k>();
            tzo.constructSpanner(*G, spanner_tz, k);
            std::cout << currentDateTime() << " " << "TZ done! G.getE().size() = " << G->getE().size() << ", spanner_tz.getE().size() = " << spanner_tz.getE().size() << std::endl;
            sizeMeasurements[1].push_back(spanner_tz.getESize());

            outputResults<unsigned>(outPath, sizeMeasurements);
            delete G;
        }
    }
}
