#include <iostream>
#include <random>
#include <ctime>

#include "heap.hpp"
#include "structs.hpp"


void basicTest(){
    std::cout << "Heap basicTest init" << std::endl;
    int size = 1000;
    std::vector<std::pair<Vertex*, double> > inpV; //id,value
    Heap<Vertex > heap = Heap<Vertex >(std::less<double>());

    inpV.reserve(size);

    for(unsigned i = 0; i < size; ++i){
        // int *&a = new int(i);
        inpV.push_back({new Vertex(i),double(i)});
    }
    heap.buildHeap(inpV);

    //test ptrs are set
    
    //// print addresses
    // for(unsigned i = 0; i < size; ++i){
    //     std::cout << printVertex(*(inpV[i].first)) << std::endl;
    // }

    for(unsigned i = 0; i < size; ++i){
        std::pair<Vertex,double> a = heap.pop();
        // if(a.second != i){
        //     std::cout << "Heap correctness test: Error in basicTest " << std::endl;
        // }
        if(a.first.id != i && a.second != double(i)){
            std::cout << "Heap correctness test: Error in basicTest " 
            << a.first.id << " != " << i << std::endl;
        } 
    }

    //cleanup
    for(unsigned i = 0; i < size; ++i){
        // int *&a = new int(i);
        delete inpV[i].first;
    }
    std::cout << "Heap basicTest exit" << std::endl;
    // return 1;
}

void randTest(){
    std::srand(std::time(0)); // use current time as seed for random generator
    std::cout << "Heap randTest init" << std::endl;
    int size = 1000;
    std::vector<std::pair<Vertex*, double> > inpV; //id,value
    Heap<Vertex > heap = Heap<Vertex >(std::less<double>());
    inpV.reserve(size);

    for(unsigned i = 0; i < size; ++i){
        inpV.push_back({new Vertex(i), std::rand() % 100});
    }
    heap.buildHeap(inpV);

    // sort using a custom function object
    struct {
        bool operator()(std::pair<Vertex*, double> a, 
                        std::pair<Vertex*, double> b)
        {   
            return a.second < b.second;
        }   
    } customLess;
    std::sort(inpV.begin(), inpV.end(), customLess);

    for(unsigned i = 0; i < size; ++i){
        std::pair<Vertex,double> a = heap.pop();

        if(a.first.id != inpV[i].first->id && a.second != inpV[i].second){
            std::cout << "Heap correctness test: Error in randTest " 
            << a.second << " != " << inpV[i].second << std::endl;
        } 
    }

    //cleanup
    for(unsigned i = 0; i < size; ++i){
        // int *&a = new int(i);
        delete inpV[i].first;
    }
    // return 1;
    std::cout << "Heap randTest exit" << std::endl;
}


//todo test for changing priority (?)
//todo test for correct handle pointers?



void testHeap(){
    std::cout << "Testing heap correctness" << std::endl;
    basicTest();
    randTest();
    std::cout << "Testing heap correctness done" << std::endl;
}