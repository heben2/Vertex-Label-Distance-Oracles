#include <string>
#include <fstream>
#include <vector>
#include "structs.hpp"

#include "tst_methods.hpp"

#ifndef PARSER_H
#define PARSER_H

namespace parser {
    //load graph from given path. Return Graph object.
    Graph* loadGraph(std::string path);

    //construct random full graph and store in file path
    void createFullGraph(std::string path, unsigned numVertices);

    //construct random full graph with random label assignment in file path
    void createFullLabelGraph(std::string path, unsigned numVertices, 
                              unsigned numLabels);
	Graph* returnFullLabelGraph(std::string path, unsigned numVertices,
							    unsigned numLabels);
    void createLabelGraph(std::string path, unsigned numVertices, 
                          unsigned numEdges,unsigned numLabels);
    Graph* returnLabelGraph(std::string path, unsigned numVertices, 
                            unsigned numEdges,unsigned numLabels);

    void constructLabelGraph(unsigned numVertices, unsigned numEdges, 
                             unsigned numLabels, Graph& G);

    //empties/clear the file at path
    void clearFile(std::string path);

    void loadGraphFromFileGr(std::string path, Graph& G, unsigned numLabels);
}

#endif /* PARSER_H */