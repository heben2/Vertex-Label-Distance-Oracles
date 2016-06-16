#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "oracles/lib/structs.hpp"

#ifndef TST_METHODS
#define TST_METHODS

typedef std::vector<std::vector<double> > ResultList;

struct TestParams{
    unsigned minN, maxN, stepSizeN, minK, maxK, stepSizeK, repeatIter,
             minL=0, maxL=0, stepSizeL=0,
             l=0, n=0;
    TestParams(unsigned minN, unsigned maxN, unsigned stepSizeN,
               unsigned minK, unsigned maxK, unsigned stepSizeK,
               unsigned const l, unsigned repeatIter);
    
    // TestParams(unsigned minK, unsigned maxK, unsigned stepSizeK,
    //            unsigned minL, unsigned maxL, unsigned stepSizeL,
    //            unsigned const n, unsigned repeatIter);

    TestParams(unsigned minN, unsigned maxN, unsigned stepSizeN,
               unsigned minK, unsigned maxK, unsigned stepSizeK,
               unsigned minL, unsigned maxL, unsigned stepSizeL,
               unsigned const n, unsigned repeatIter);
};

void getShortestDistances(Graph& G, ResultList& resultsTrue);

void verifyStretchBound(ResultList& resultsTrue, ResultList& resultsOther, 
                        const unsigned stretch, bool print);

void initResultsVector(ResultList* results, unsigned const n);


void printResults(std::string path, 
                      std::vector< std::vector< double> > results,
                      TestParams tstParams);



//output TestParams;
void outputParams(std::string path, TestParams tp);

//outputs "description: p"
template <typename T>
void outputParam(std::string path, std::string description, T p){
    std::ofstream f;
    f.open(path, std::fstream::out | std::fstream::app);
    f << description + " : " << p << std::endl;
    f.close();
}


//expects simple types T. Appends to output file.
template <typename T>
void outputResults(std::string path, std::vector<T>& list){
    std::ofstream f;
    f.open(path, std::fstream::out | std::fstream::app);
    std::string s = "";
    f << "[";
    for(auto a : list){
        s += std::to_string(a) + ", ";
    }
    s = s.substr(0, s.size()-2);
    f << s + "]" << std::endl;
    f.close();
}

template <typename T>
void outputResults(std::string path, std::vector< std::vector< T> >& list){
    std::ofstream f;
    f.open(path, std::fstream::out | std::fstream::app);

    std::string s = "";
    f << "[";
    for(auto a : list){
        s += "[";
        for(auto b : a){
            s += std::to_string(b) + ", ";
        }
        s = s.substr(0, s.size()-2);
        s += "], \n";
    }
    s = s.substr(0, s.size()-3);
    f << s + "]" << std::endl;
    f.close();
}


void outputResultPairs(std::string path, std::vector< std::vector< std::pair<double, double>> >& list);
//     std::ofstream f;
//     f.open(path, std::fstream::out | std::fstream::app);

//     std::string s = "";
//     f << "[";
//     for(auto a : list){
//         s += "[";
//         for(auto b : a){
//             s += "(" + std::to_string(b.first) + ", " 
//                      + std::to_string(b.second) + "), ";
//         }
//         s = s.substr(0, s.size()-2);
//         s += "], \n";
//     }
//     s = s.substr(0, s.size()-3);
//     f << s + "]" << std::endl;
//     f.close();
// }




#endif /* TST_METHODS */