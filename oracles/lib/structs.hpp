#include <exception>
#include <algorithm>
#include <unordered_set>
#include <cstdint>
#include <string>
#include <memory>
#include <unordered_map>
#include <map>


#ifndef STRUCTS_H
#define STRUCTS_H

class Vertex {
    public:
        unsigned id; //name/identifier of vertex
        unsigned label;
        //TODO change ptr to be a view pointer with template type(?) or other smart pointer (which can be empty!!)
        uintptr_t ptr = 0; //castable pointer. Ugly, but works for changeable pointer types
        Vertex* witness = nullptr; //used in trees and as witnesses

        Vertex();
        Vertex(unsigned id, unsigned labelId = 0);
        //sets itself as parent
        Vertex(unsigned id, bool setSelfWitness, unsigned labelId = 0);
        Vertex(Vertex const &v, bool setSelfWitness);
        Vertex(Vertex const &v);
        ~Vertex();
        void setWitnessSelf();
        // unsigned& get_id(){ return id; }
        bool operator< (Vertex const& v) const;
        bool operator== (Vertex const& v) const;
};

std::string printVertex(Vertex& v);


struct EdgeInitException : public virtual std::exception{
    const char * what () const throw ();
};

class Edge {

    public:
        Vertex u;
        Vertex v;
        double weight; //TODO make double or long (sum of )
        Edge(Vertex u, Vertex v, double w);
        Edge();
        Edge& operator=(const Edge& e);

        //how to order in sets,
        bool operator< (const Edge& e) const;

        bool operator== (const Edge& e) const;

        const unsigned get_id_u() const;
        const unsigned get_id_v() const;
        std::string const getKey() const;
};


struct Hash {
    std::size_t operator() (const Vertex &v) const;

    std::size_t operator() (const Edge &e) const;
};




typedef std::unordered_map<std::string, Edge> EdgeMap;
typedef std::map<unsigned, Vertex> VertexMap;

struct Graph{
    std::unique_ptr< VertexMap > V;
    std::unique_ptr< EdgeMap > E;
    std::unique_ptr< std::unordered_set<unsigned> > L;

    Graph();
    ~Graph();
    Graph(Graph &g);
    Graph(Graph &&){};
    Graph(const Graph&);
    Graph(Graph *g);
    Graph(std::map<unsigned, Vertex>* V,
          std::unordered_map<std::string, Edge>* E);
    Graph(std::map<unsigned, Vertex>* V,
          std::unordered_map<std::string, Edge>* E,
          std::unordered_set<unsigned>* L);

    Graph(std::map<unsigned, Vertex>& V,
          std::unordered_map<std::string, Edge>& E,
          std::unordered_set<unsigned>& L);

    Graph& operator=(const Graph&);

    VertexMap& getV();
    EdgeMap& getE();
    std::unordered_set<unsigned>& getL();

    void reset(VertexMap& V, EdgeMap& E, std::unordered_set<unsigned>& L);

    void setV(VertexMap& V);
    void setE(EdgeMap& E);
    void setL(std::unordered_set<unsigned>& L);
    unsigned getLSize();
    unsigned getVSize();
    unsigned getESize();
};

//v.id->dist to v
typedef std::unordered_map< unsigned, double> VertexDistMap;
//v.id -> u.id -> dist(u,v)
typedef std::unordered_map< unsigned, std::unique_ptr< VertexDistMap > > BMap;
//v.id -> p_i(v), dist(v,p_i(v))
typedef std::unordered_map< unsigned, std::pair<Vertex, double> > WitnessMap;

//v.id -> edges (e.g. incident to v)
typedef std::unordered_map< unsigned, 
                            std::unordered_set<Edge, Hash> > AdjacentEdgeMap;


#endif /* STRUCTS_H */