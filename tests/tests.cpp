#include <iostream>

#include <chrono>
#include <ctime>
#include <iomanip>
#include "../oracles/thorup_zwick_oracle.hpp"

#include "../oracles/lib/tst_methods.hpp"
#include "oracles/lib/heap.test.cpp"
#include "oracles/lib/parser.test.cpp"
#include "stretch_thorup_zwick.cpp"
#include "stretch_vertex_label_oras.cpp"
// #include "spanner.test.cpp"

void correctnessTests(){
    std::cout << "correctnessTests: init" << std::endl;
    // testHeap();
    std::cout << "runTestSuite: test heap done" << std::endl;
    // testParser();
    std::cout << "runTestSuite: test parser done" << std::endl;

    unsigned minN = 100, maxN = 1000, stepSizeN = 100,
             minK = 2,  maxK = 10,  stepSizeK = 2;
    testThorupZwick(minN,maxN,stepSizeN,minK,maxK,stepSizeK);
    std::cout << "runTestSuite: test Thorup-Zwick done" << std::endl;

    std::cout << "correctnessTests: done" << std::endl;
}



// Just for testing cpu time vs wall clock time
void tmpTestTime(){
    unsigned const k = 2;
    Graph* G = parser::loadGraph("tmp.txt");
    for(unsigned i=0; i < 10; ++i){
        auto t_start = std::chrono::high_resolution_clock::now();
        std::clock_t c_start = std::clock();
        ThorupZwickOracle tzo;// = ThorupZwickOracle<k>();
        tzo.preproc(*G, k);
        auto t_end = std::chrono::high_resolution_clock::now();
        std::clock_t c_end = std::clock();
    
        std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
                  << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n"
                  << "Wall clock time passed: "
                  << std::chrono::duration<double, std::milli>(t_end-t_start).count()
                  << " ms\n";
    }
    delete G;
}



void testStretch(){
        // // testing stretch accuracy. Here over n l k (not l n k as other tests)
        std::cout << "runTestSuite::testStretch" << std::endl;
        std::string outPathStretch = "stretchResults.json";
        unsigned minN = 800, // (2) increase to 1000
        maxN = 800, // (2) increase to 4000
        stepSizeN = 1,
        minK = 5,
        maxK = 12,
        stepSizeK = maxK-minK,
        repeatIter = 20,
        l = 0;

        TestParams tstParams(minN, maxN, stepSizeN, minK, maxK, stepSizeK, l, repeatIter);

        parser::clearFile(outPathStretch);
        for(unsigned n=minN; n <= maxN; n += stepSizeN){
            tstParams.n = n;
            // tstParams.minL = floor(n/2);
            // tstParams.maxL = floor(3*n/4);
            tstParams.minL = 8;//floor(double(n)/10.);
            tstParams.maxL = 8;//floor(double(3*n)/4.);//floor(n/3);
            tstParams.stepSizeL = 1;//tstParams.maxL - tstParams.minL;//only run for min and max L
            vertexLabelDistOras::testStretch(tstParams, outPathStretch);
        }

}

void testSpannerSize(){
    unsigned n = 1000,
             minK = 2, 
             maxK = 8, 
             stepSizeK = 6,
             repeatIter = 20;
    testSpannerSize(n,minK,maxK,stepSizeK,repeatIter);
}

void testOracleSize(){
    std::string outPathSize = "sizeResultsPreproc.json";
    unsigned //for longer tests (1) small maxN
         minN = 800,
         maxN = 800,
         stepSizeN = 1,
         minK = 5,
         maxK = 23,
         stepSizeK = 9,
         repeatIter = 20,
         l = 5;
    TestParams tstParams(minN, maxN, stepSizeN, minK, maxK, stepSizeK, l, repeatIter);

    // static const unsigned arr[] = {8,300,600};
    // std::vector<unsigned> L (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    std::vector<unsigned> L = {8,300,600};
    for(auto l : L){
        tstParams.l = l;
        tstParams.minN = tstParams.minN < l? l : tstParams.minN;
        vertexLabelDistOras::testPreprocSize(tstParams, outPathSize);
    }

}

void testRealGraph(){
    std::string graphPath = "USA-road-d.NY.gr",
                outPathTime = "timeResultsPreprocTime.json";
    parser::clearFile(outPathTime);
    unsigned numLabels = 100;
    unsigned k = 14;
    vertexLabelDistOras::testPreprocTimeRealGraph(graphPath, k, numLabels, outPathTime);
}


void runTestSuite(){
    std::cout << "runTestSuite init" << std::endl;

    std::string outPathStretch = "stretchResults.json",
                outPathTime = "timeResultsPreprocTime.json",
                outPathSize = "sizeResultsPreproc.json";
	

    testOracleSize();

	 // unsigned //for longer tests (1) small maxN
	 // 	minN = 200, // (2) increase to 1000
	 // 	maxN = 1000, // (2) increase to 4000
	 // 	stepSizeN = 100,
	 // 	minK = 5,
	 // 	maxK = 12,
	 // 	stepSizeK = 2,
	 // 	minL = 10,
	 // 	maxL = 200,
	 // 	stepSizeL = 190,
	 // 	repeatIter = 5,
	 // 	l = 5;
		
    // TestParams tstParams(minN, maxN, stepSizeN, minK, maxK, stepSizeK, l, repeatIter);
    

    //unsigned //for short tests
    //   minN = 50, 
    //   maxN = 1000, 
    //   stepSizeN = 50,
    //   minK = 7,
    //   maxK = 12,
    //   stepSizeK = 1,
    //   minL = 10,
    //   maxL = 40,
    //   stepSizeL = 10,
    //   repeatIter = 5,
    //   l = 5;
    //TODO increase repeatIter

    //// Iterate over l values
    // for(unsigned l = minL; l < lMax; ++l)
    //     TestParams tstParams(minN, maxN, stepSizeN, minK, maxK, stepSizeK, l, repeatIter);
    //     vertexLabelDistOras::testVertexLabelDistOras(tstParams);
    // }

    // correctnessTests(); //testing heap, parser and thorup-zwick (missing sssp tests)
    // parser::clearFile(outPathStretch);
    // for(unsigned l=minL; l < maxL; l += stepSizeL){
    //     tstParams.l = l;
    //     tstParams.minN = tstParams.minN < l? l : tstParams.minN;
    //     vertexLabelDistOras::testStretchCorrectness(tstParams, outPathStretch);
    // }


    // testStretch();



    // testRealGraph();
    

    // // testing only time
    // parser::clearFile(outPathTime);
    // parser::clearFile(outPathSize);
    // for(unsigned l=minL; l <= maxL; l += stepSizeL){
    //     tstParams.l = l;
    //     tstParams.minN = tstParams.minN < l? l : tstParams.minN;
    //     vertexLabelDistOras::testPreprocTime(tstParams, outPathTime, 
    //                                          outPathSize);
    // }


    //for short tests
    // minN = 400,
    // maxN = 2000,
    // stepSizeN = 100,
    // minK = 2,
    // maxK = 8,
    // stepSizeK = 1,
    // repeatIter = 5,
    // l = 5;
    // tstParams = TestParams(minN, maxN, stepSizeN, minK, maxK, stepSizeK, l, repeatIter);

    // testSpanner(tstParams);

    // tmpTestTime();

    std::cout << "runTestSuite: end" << std::endl;
}