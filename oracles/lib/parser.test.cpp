#include <iostream>
#include <random>
#include <ctime>
#include <string>

#include "parser.hpp"
#include "structs.hpp"



void testLoadGraph(){
    std::string path = "tmp.txt";
    int numVert = 10;
    int edgeCount = numVert*(numVert-1)/2;
    parser::createFullGraph(path, numVert);
    Graph* G = parser::loadGraph(path);

    for(int i = 0; i < numVert; ++i)
        if(G->V->find(i) == G->V->end())
            std::cout << "Error: testLoadGraph(), vertex id " << i << " not found" << std::endl;

    int count = 0;
    for(auto e : (*G->E) ){
        count++;
    }
    if(count != edgeCount)
        std::cout << "Error: testLoadGraph(), edge count = " << count << " but should be " << edgeCount << std::endl;

}



void testParser(){
    std::cout << "Testing parser correctness" << std::endl;
    testLoadGraph();
    std::cout << "Testing parser correctness done" << std::endl;
}