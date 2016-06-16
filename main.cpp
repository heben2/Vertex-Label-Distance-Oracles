#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include "tests/tests.cpp"
#include "oracles/lib/parser.hpp"
#include "oracles/thorup_zwick_oracle.hpp"
#include "oracles/chechik_oracle.hpp"
#include "oracles/with_dijkstra_oracle.hpp"
#include "oracles/with_thorup_zwick_oracle.hpp"


/**
* @function main
*/
int main(int argc, char** argv) {
    // if (argc != 2) {
    //     cout << "No image provided" << endl;
    //     return -2;
    // }
    /*if (argc != 3) {
        cout << "No image provided" << endl;
        return -2;
    }
*/
    std::cout << "main: init" << std::endl;
    //CORRECTNESS tests
    runTestSuite();
    std::cout << "main: tests are done" << std::endl;

    // testThorupZwickOra();
    // testChechik();
    // testThorupZwickSpanner();
    // testWithDijkstra();
    // testWithThorupZwick();

    // parser::createFullGraph("tmp.txt", 10);


    //Experiments

    return 1;
}