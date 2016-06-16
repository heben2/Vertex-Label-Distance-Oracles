// #include <exception>
// #include <algorithm>
// #include <unordered_set>
#include <iostream>


#include "structs.hpp"

Vertex::Vertex(){}
Vertex::Vertex(unsigned id, unsigned label):id(id), label(label){}
Vertex::Vertex(Vertex const &v):id(v.id), ptr(v.ptr), witness(v.witness), 
    label(v.label){}
Vertex::Vertex(unsigned id, bool setSelfWitness, unsigned label):
    id(id), label(label){
    if(setSelfWitness)
        witness = this;
}
Vertex::Vertex(Vertex const &v, bool setSelfWitness): 
id(v.id), ptr(v.ptr), witness(v.witness), label(v.label){
    if(setSelfWitness)
        witness = this;
}

Vertex::~Vertex(){ }

void Vertex::setWitnessSelf(){
    witness = this;
}
bool Vertex::operator< (const Vertex& v) const {
    return id < v.id;
}
bool Vertex::operator== (const Vertex& v) const {
    return id == v.id;
}

std::string printVertex(Vertex& v){
    return "Vertex.id = "+std::to_string(v.id)+", Vertex.ptr = "
            +std::to_string(v.ptr);
}





const char * EdgeInitException::what () const throw () {
    return "Invalid initialization of Edge class \
    (possibly same vertices given).";
}
Edge::Edge(){}
Edge::Edge(Vertex u, Vertex v, double w): u(u), v(v), weight(w) {}

Edge& Edge::operator=(const Edge& e){
    u = e.u;
    v = e.v;
    weight = e.weight;    
    return *this;
}


bool Edge::operator< (const Edge& e) const {
    return this->getKey() < e.getKey();
}

bool Edge::operator== (const Edge& e) const {
    using std::min;
    using std::max;
    return min(u,v) == min(e.u,e.v) && max(u,v) == max(e.u,e.v);
}

const unsigned Edge::get_id_u() const{
    return u.id;
}
const unsigned Edge::get_id_v() const{
    return v.id;
}

//Edges are undirected: keys are ensured aligned even if v and u are swapped.
std::string const Edge::getKey() const {
    return std::to_string(std::min(u.id,v.id)) + "," 
           + std::to_string(std::max(u.id,v.id));
}


std::size_t Hash::operator() (const Vertex &v) const {
    return v.id;//(hash<unsigned>()(v.id));//v.id; //not used anymore
}

//todo test this
std::size_t Hash::operator() (const Edge &e) const {
    using std::hash;
    std::size_t h1 = hash<std::string>()(std::to_string(std::min(e.u.id, e.v.id)));
    std::size_t h2 = hash<std::string>()(std::to_string(std::max(e.u.id, e.v.id)));
    return h1 ^ (h2 << 1);
}


//TODO does this work?
Graph::Graph():
    V(new std::map<unsigned, Vertex>()),
    E(new std::unordered_map<std::string, Edge>()),
    L(new std::unordered_set<unsigned>()){}


Graph::Graph(Graph& g):
    V(new std::map<unsigned, Vertex>(*(g.V.get())) ),
    E(new std::unordered_map<std::string, Edge>(*(g.E.get())) ),
    L(new std::unordered_set<unsigned>(*(g.L.get())) ){}

Graph::Graph(Graph* g):
    V(new std::map<unsigned, Vertex>(*(g->V.get())) ),
    E(new std::unordered_map<std::string, Edge>(*(g->E.get())) ),
    L(new std::unordered_set<unsigned>(*(g->L.get())) ){}

Graph::Graph(std::map<unsigned, Vertex>* V, 
             std::unordered_map<std::string, Edge>* E): 
    V(V), E(E), L(new std::unordered_set<unsigned>()){}

Graph::Graph(std::map<unsigned, Vertex>* V, 
             std::unordered_map<std::string, Edge>* E,
             std::unordered_set<unsigned>* L): 
    V(V), E(E), L(L){}

Graph::Graph(std::map<unsigned, Vertex>& V,
             std::unordered_map<std::string, Edge>& E,
             std::unordered_set<unsigned>& L):
    V(new std::map<unsigned, Vertex>(V)),
    E(new std::unordered_map<std::string, Edge>(E)),
    L(new std::unordered_set<unsigned>(L)){}

Graph& Graph::operator=(const Graph& g){
    V.reset(new std::map<unsigned, Vertex>(*(g.V.get())) );
    E.reset(new std::unordered_map<std::string, Edge>(*(g.E.get())) );
    L.reset(new std::unordered_set<unsigned>(*(g.L.get())) );
    return *this;
}

Graph::~Graph(){
    V.reset();
    E.reset();
    L.reset();
}

std::map<unsigned, Vertex>& Graph::getV(){
    return *V.get();
}
std::unordered_map<std::string, Edge>& Graph::getE(){
    return *E.get();
}
std::unordered_set<unsigned>& Graph::getL(){
    return *L.get();
}

unsigned Graph::getVSize(){
    return V.get()->size();
}

unsigned Graph::getLSize(){
    return L.get()->size();
}

unsigned Graph::getESize(){
    return E.get()->size();
}

void Graph::setV(VertexMap& VIn){
    V.reset(new VertexMap(VIn));
}
void Graph::setL(std::unordered_set<unsigned>& LIn){
    L.reset(new std::unordered_set<unsigned>(LIn));
}

void Graph::setE(EdgeMap& EIn){
    E.reset(new EdgeMap(EIn));
}

void Graph::reset(VertexMap& V, EdgeMap& E, std::unordered_set<unsigned>& L){
    setV(V);
    setE(E);
    setL(L);
}