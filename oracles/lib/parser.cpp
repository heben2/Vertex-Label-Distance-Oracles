#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <random>
#include <set>

#include "parser.hpp"
#include "structs.hpp"
#include "std_methods.hpp" //for some reason split is not found from std_methods.h...

// //copied from http://stackoverflow.com/questions/236129/split-a-string-in-c
// std::vector<std::string> &split(const std::string &s, char delim,
//                                 std::vector<std::string> &elems) {
//     std::stringstream ss(s);
//     std::string item;
//     while (std::getline(ss, item, delim)) {
//         elems.push_back(item);
//     }
//     return elems;
// }


// std::vector<std::string> split(const std::string &s, char delim) {
//     std::vector<std::string> elems;
//     split(s, delim, elems);
//     return elems;
// }


Graph* parser::loadGraph(std::string path){
	std::cout << "parser::loadGraph begin" << std::endl;
    std::ifstream f(path);
    Graph* G = new Graph();
    std::string line;
    while(std::getline(f,line) && line != ""){
        std::vector<std::string> s = split(line, ' ');
        unsigned i = std::stoi(s[0]);
        unsigned l = 0;
        if(s.size() > 1)
            l = std::stoi(s[1]);
        // (*(G->V))[i] = Vertex(i, l);
        G->getV()[i] = Vertex(i, l);
        G->getL().insert(l);
    }
    while(std::getline(f,line)){
        std::vector<std::string> s = split(line, ' ');
        std::string::size_type sz;
        // double d = ,&sz);
        Edge e(Vertex(std::stoi(s[0])), Vertex(std::stoi(s[1])), 
                      std::stod(s[2]));
        
        (*(G->E))[e.getKey()] = e;
    }
	std::cout << "parser::loadGraph begin" << std::endl;
    f.close();
    return G;
}



void parser::createFullGraph(std::string path, unsigned numVertices){
    std::srand(std::time(0));
    std::ofstream f;
    f.open(path);
	std::cout << "parser::createFullGraph begin" << std::endl;
    for(unsigned i=0; i < numVertices; ++i)
        f << std::to_string(i) + "\n";

    f << "\n";

    for (unsigned i = 0; i < numVertices; ++i)
        for (unsigned j = 0; j < i; ++j){
                double dist = (std::rand() % 100)+1;
                f << std::to_string(i) + " " + std::to_string(j) + " " 
                                                + std::to_string(dist) + "\n";
            }
	std::cout << "parser::createFullGraph done" << std::endl;
    f.close();
}


void parser::createFullLabelGraph(std::string path, unsigned numVertices, 
                                  unsigned numLabels){
    std::srand(std::time(0));
    std::random_device rd;
    std::ranlux48_base e2(rd());
    std::uniform_int_distribution<> distLabels(0, numLabels-1);
    std::uniform_int_distribution<> distWeights(1, 100);

    std::ofstream f;
    f.open(path);
	std::cout << "parser::createFullLabelGraph begin" << std::endl;

    //setup vertices
    //make sure each label is assigned at least once.
    for(unsigned i=0; i < numLabels; ++i)
        f << std::to_string(i) + " " << i << "\n";
    for(unsigned i=numLabels; i < numVertices; ++i)
        f << std::to_string(i) + " " << distLabels(e2) << "\n";
    f << "\n";
    //setup edges
    for (unsigned i = 0; i < numVertices; ++i)
        for (unsigned j = 0; j < i; ++j){
                double weigt = distWeights(e2);//(std::rand() % 100)+1;
                f << std::to_string(i) + " " + std::to_string(j) + " " 
                                                + std::to_string(weigt) + "\n";
            }
    f.close();
	std::cout << "parser::createFullLabelGraph end" << std::endl;
}


Graph* parser::returnFullLabelGraph(std::string path, unsigned numVertices,
									unsigned numLabels) {
	std::srand(std::time(0));
	std::random_device rd;
	std::ranlux48_base e2(rd());
	std::uniform_int_distribution<> distLabels(0, numLabels - 1);
	std::uniform_int_distribution<> distWeights(1, 100);
	Graph* G = new Graph();
	std::ofstream f;
	f.open(path);
	std::cout << "parser::returnFullLabelGraph begin" << std::endl;

    VertexMap& V = G->getV();
    EdgeMap& E = G->getE();
    std::unordered_set<unsigned>& L = G->getL();

	//setup vertices
	//make sure each label is assigned at least once.
	for (unsigned i = 0; i < numLabels; ++i) {
		f << std::to_string(i) + " " << i << "\n";
		V[i] = Vertex(i, i);
		L.insert(i);
	}
	for (unsigned i = numLabels; i < numVertices; ++i) {
		unsigned l = distLabels(e2);
		f << std::to_string(i) + " " << l << "\n";
		V[i] = Vertex(i, l);
		L.insert(l);
	}
	std::cout << "parser::returnFullLabelGraph vertices done, now creating " << numVertices*(numVertices-1)/2 << " edges." << std::endl;
	f << "\n";
	//setup edges
    for (unsigned i = 0; i < numVertices; ++i)
        for (unsigned j = 0; j < i; ++j){
				double weight = distWeights(e2);//(std::rand() % 100)+1;
				Edge e = Edge(V[i], V[j], weight);
				f << std::to_string(i) + " " + std::to_string(j) + " "
					+ std::to_string(weight) + "\n";
				E[e.getKey()] = e;
			}
	f.close();
	std::cout << "parser::returnFullLabelGraph end" << std::endl;
	return G;
}



/* Constructs random graph with random edges and random weights.
   Might construct a not fully-connected graph.
*/
void parser::createLabelGraph(std::string path, unsigned numVertices, 
                              unsigned numEdges,unsigned numLabels){
    
    std::srand(std::time(0));
    std::random_device rd;
    std::ranlux48_base e2(rd());
    std::uniform_int_distribution<> distLabels(0, numLabels-1);
    std::uniform_int_distribution<> distWeights(1, 100);

    if(numEdges > numVertices*(numVertices-1)/2)//must not be larger than full graph
        numEdges = numVertices*(numVertices-1)/2;

    std::ofstream f;
    f.open(path);

    //setup vertices
    //make sure each label is assigned at least once.
    for(unsigned i=0; i < numLabels; ++i)
        f << std::to_string(i) + " " << i << "\n";
    for(unsigned i=numLabels; i < numVertices; ++i)
        f << std::to_string(i) + " " << distLabels(e2) << "\n";
    f << "\n";

    //setup edges
    std::vector< Edge > edges;
    for (unsigned i = 0; i < numVertices; ++i)
        for (unsigned j = 0; j < i; ++j)
            edges.push_back(Edge(Vertex(i), Vertex(j), distWeights(e2)));
    
    //randomly select edges.
    if(edges.size() > numEdges){
        std::random_shuffle(edges.begin(), edges.end()); //use default random gen
        //select first numEdges elements.
        edges.resize(numEdges);
    }
    //print edges
    for(auto& e : edges)
        f << std::to_string(e.u.id) + " " + std::to_string(e.v.id) + " " 
                                                + std::to_string(e.weight) + "\n";
    f.close();
}

Graph* parser::returnLabelGraph(std::string path, unsigned numVertices, 
                                unsigned numEdges,unsigned numLabels){
    std::srand(std::time(0));
    std::random_device rd;
    std::ranlux48_base e2(rd());
    std::uniform_int_distribution<> distLabels(0, numLabels-1);
    std::uniform_int_distribution<> distWeights(1, 100);
    Graph* G = new Graph();
    VertexMap& V = G->getV();
    EdgeMap& E = G->getE();
    std::unordered_set<unsigned>& L = G->getL();

    std::ofstream f;
    f.open(path);

    std::cout << "numVertices = " << numVertices << std::endl;
    //setup vertices
    //make sure each label is assigned at least once.
    for(unsigned i=0; i < numLabels; ++i){
        f << std::to_string(i) + " " << i << "\n";
            V[i] = Vertex(i, i);
            L.insert(i);
    }
    for(unsigned i=numLabels; i < numVertices; ++i){
        unsigned l = distLabels(e2);
        f << std::to_string(i) + " " << l << "\n";
        V[i] = Vertex(i, l);
        //L.insert(l);
    }
    f << "\n" << std::endl;

    //setup edges
    std::vector< Edge > edges;
    for (unsigned i = 0; i < numVertices; ++i)
        for (unsigned j = 0; j < i; ++j)
            edges.push_back(Edge(V[i], V[j], distWeights(e2)));
    

    //randomly select edges.
    if(edges.size() > numEdges){
        std::random_shuffle(edges.begin(), edges.end()); //use default random gen
        //select first numEdges elements.
        edges.resize(numEdges);
    }
    //print edges
    for(auto& e : edges){
        E[e.getKey()] = e;
        f << std::to_string(e.u.id) + " " + std::to_string(e.v.id) + " " 
                                                + std::to_string(e.weight) + "\n";
    }
    f.close();
    return G;
}


void parser::constructLabelGraph(unsigned numVertices,unsigned numEdges, unsigned numLabels,
                         Graph& G) {
    std::srand(std::time(0));
    std::random_device rd;
    std::ranlux48_base e2(rd());
    std::uniform_int_distribution<> distLabels(0, numLabels - 1);
    std::uniform_int_distribution<> distWeights(1, 100);
    VertexMap& V = G.getV();
    EdgeMap& E = G.getE();
    std::unordered_set<unsigned>& L = G.getL();

    //setup vertices
    //make sure each label is assigned at least once.
    //for (unsigned i = 0; i < numLabels; ++i) {
    //    V[i] = Vertex(i, i);
    //    L.insert(i);
    //}
    //for (unsigned i = numLabels; i < numVertices; ++i) {
    //    unsigned l = distLabels(e2);
    //    V[i] = Vertex(i, l);
    //    //L.insert(l);
    //}

    for (unsigned i = 0; i < numVertices; ++i) {
        unsigned l = distLabels(e2);
        V[i] = Vertex(i, l);
    }
    for (unsigned i = 0; i < numLabels; ++i) {
        V[i].label = i;
        L.insert(i);
    }

    //setup edges
    std::vector< Edge >* edges = new std::vector<Edge>();
    edges->reserve(numVertices);
    for (unsigned i = 0; i < numVertices; ++i) {
        for (unsigned j = 0; j < i; ++j) {
            edges->push_back(Edge(V[i], V[j], distWeights(e2)));
        }
    }
            
    //randomly select edges.
    if (edges->size() > numEdges) {
        std::random_shuffle(edges->begin(), edges->end()); //use default random gen
        edges->resize(numEdges);
    }
    //print edges
    for (auto& e : *edges) {
        E[e.getKey()] = e;
    }
    delete edges;
}

void parser::clearFile(std::string path){
    std::ofstream f;
    f.open(path);
    f.close();
}



void parser::loadGraphFromFileGr(std::string path, Graph& G, 
                                 unsigned numLabels){
    std::ifstream f(path);
    std::srand(std::time(0));
    std::random_device rd;
    std::ranlux48_base e2(rd());
    std::uniform_int_distribution<> distLabels(0, numLabels - 1);
    VertexMap& V = G.getV();
    EdgeMap& E = G.getE();
    std::unordered_set<unsigned>& L = G.getL();

    std::string line;
    std::vector<std::string> s;
    while(std::getline(f,line) && line[0] != 'p'){}
    s = split(line, ' ');
    unsigned n = std::stoi(s[2]),
             m = std::stoi(s[3]);

    // V.reserve(n);
    E.reserve(m);
    while(std::getline(f,line) && line[0] != 'a'){}

    // unsigned lAssigned = 0;
    do{
        s = split(line, ' ');
        unsigned vid = std::stoi(s[1]),
                 uid = std::stoi(s[2]);
        double weight = std::stoi(s[3]);
        auto it = V.find(vid);
        if(it == V.end()){
            // unsigned l = lAssigned < numLabels? lAssigned++ : distLabels(e2);
            unsigned l = distLabels(e2);
            V[vid] = Vertex(vid, l);
        }
        it = V.find(uid);
        if(it == V.end()){
            unsigned l = distLabels(e2);
            V[vid] = Vertex(uid, l);
        }
        Edge e = Edge(uid, vid, weight);
        E[e.getKey()] = e;
    } while(std::getline(f,line) && line != "");

    std::uniform_int_distribution<> distLabelsVertices(1, n);
    std::set<int> vAlreadySet;
    int i;
    for(unsigned l = 0; l < numLabels; ++l){
        L.insert(l);
        do{//make sure all labels are asigned (at random)
            i = distLabelsVertices(e2);
        }while(vAlreadySet.find(i) != vAlreadySet.end());
        vAlreadySet.insert(i);
        V[i].label = l;
    }

    std::cout << "n = " << G.getVSize() << ", m = " << G.getESize() << ", l = " << G.getLSize() << std::endl;

    f.close();
}