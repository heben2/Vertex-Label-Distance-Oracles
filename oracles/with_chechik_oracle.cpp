#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <unordered_map>
#include <float.h>
#include "with_chechik_oracle.hpp"
#include "spanner.hpp"
#include "lib/structs.hpp"
#include "lib/std_methods.hpp"


#include <chrono>


WithChechikOracle::WithChechikOracle() {
    OGS = new ChechikOracle();
    RCO = new RestrictedChechikOracle();
}
WithChechikOracle::~WithChechikOracle() {
    std::cout << "WithChechikOracle::~WithChechikOracle end" << std::endl;
    delete OGS;
    delete RCO;
    std::cout << "WithChechikOracle::~WithChechikOracle end" << std::endl;
}


/*
Make S-set, sampling with probability p from V.
Compute E_S, make graph G_S = (V,E_S).
Compute witnesses and their distances for set S.
Make chechik-ora on G_S.

Make k'-spanner H on G.

Make restricted Chechik's oracle on H.
*/
int WithChechikOracle::preproc(Graph &G, unsigned const k){
std::cout << "WithChechikOracle::preproc begin" << std::endl;
    unsigned n = G.getVSize();
    unsigned l = G.getLSize();
    if(k < 5 || n < 1 || l < 1)
        return 0;
    unsigned checkAllIndices = l > std::pow(n,double(k)/(2.*k-1.));

////////////////// Compute k'=c and k''=kk
    unsigned c, kk; //c=k', kk=k''
    auto boundHoldsLambda = [&]() {return c >= 1 && kk >= 2
                                    && 2 + 4 * (4 * kk - 5)*(2 * c - 1) <= 4 * k - 5; };
    auto tryNewboundHoldsLambda = [&](unsigned kk_in, unsigned c_in) {
        return c_in >= 1 && kk_in >= 2
            && 2 + 4 * (4 * kk_in - 5)*(2 * c_in - 1) <= 4 * k - 5; };

    //fill out for remaining stretch bound if possible
    if(!checkAllIndices){
        std::cout << "WithChechikOracle::preproc Small label size." << std::endl;
        double const x = log(n) / log(l);

        double const f = sqrt(16 *pow(k,4)+8 *pow(k,3) *(104 *x-47)+pow(k,2)
                 *(1600 *pow(x,2)-4656 *x+2425)+54 *k*(40 *x-47)+729);
        double const P = floor(double(-4. * pow(k, 2) + k*(40. * x + 47.) - 27. + f) 
                            / double(16. * (9. * k - 5.)));

        kk = ceil( double(40.*P+4.*k-27.)/double(32.*P-16.) );
        c = floor(double(4 * k - 7) / double(8.*(4.*ceil(double(40.*P + 4.*k - 27.) 
                                                    / double(32. * P - 16.)) - 5.)) +1./2.);
    } else{
        std::cout << "WithChechikOracle::preproc Large label size." << std::endl;
        c = 2; //if possible (k > 5)
        kk = 2; //min bound
        if (!boundHoldsLambda()) { //last resort, should never get here.
            c = 1;
            kk = 2;
        }
        while (tryNewboundHoldsLambda(kk+1, c)) //increase k'' if possible
            ++kk;
        while (tryNewboundHoldsLambda(kk, c + 1)) //else increase c=k' if possible
            ++c;
    }

    std::cout << "WithChechikOracle::preproc kk=k'' = " << kk << ", c=k' = " << c << std::endl;

    
    
    double const p = std::pow(l, 1./k-1./kk-1.);
    std::cout << "WithChechikOracle::preproc kk=k'' = " << kk << ", c=k' = " << c << ", p = " << p << std::endl;
    if(p < 0 || !boundHoldsLambda())
        return 0;

//////////////////////

    
    std::map<unsigned, Vertex>& V = G.getV();
    std::unordered_map<std::string, Edge>& E = G.getE();
    std::unordered_set<unsigned>& L = G.getL();
    

std::cout << "WithChechikOracle::preproc kk=k'' = " << kk << ", c=k' = " << c << ", p = " << p << std::endl;

std::cout << "WithChechikOracle::preproc begin" << std::endl;
    std::map<unsigned, Vertex> S;
   makeRandSubset(p, V, S);


    //get/set witnesses
    witnessDistArraySize = 1;
	witnessDistArray.resize(witnessDistArraySize);
    setWitnesses(G, S, 0);



    //make E_S = UNION_v\in V E_S(v)
    //E_S(v) = the set of edges incident to v with weight < dist_G(v, p_S(v)).
    std::unordered_map<std::string, Edge>* ES = new EdgeMap();
    getEdgeSMap(E, *witnessDistArray[0], *ES);


std::cout << "WithChechikOracle::preproc S.size() = " << S.size() << ", ES.size() = " << ES->size() << std::endl;


    //Make Chechik's oracle on G_S
    Graph* GS = new Graph(V, *ES, L);
    int r = OGS->preproc(*GS, k);
    delete ES;
    delete GS;

    //Make k'-spanner using Baswana and Sen
    Graph* H = new Graph();
    // Graph H; //TODO place in mem?
    auto t1 = std::chrono::high_resolution_clock::now();
    spanner::constructSpanner(G, *H, c);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
    std::cout << "WithChechikOracle::preproc spanner construction time = " << fp_ms.count() << std::endl;

    if (S.size() > 0) {
        r = r*RCO->preproc(G, *H, S, kk, k, *witnessDistArray[0]);
    }
    else
        RCODone = false;
    delete H;
std::cout << "WithChechikOracle::preproc end" << std::endl;
    return r;
}

//d_G(u, p_S(u)) and the distance d_H(p_S(u), λ)
double WithChechikOracle::query(Vertex u, unsigned label){
    double dist1 = OGS->query(u,label);

    double dist2 = RCODone? RCO->query(u,label) : DBL_MAX;

    if(dist1 == DBL_MAX && dist2 == DBL_MAX){
        std::cout << "WithChechikOracle::query u = " << u.id << ", l = " << label << ", dist1 = dist2 = DBL_MAX" << std::endl;
    }

    return std::min(dist1, dist2);
}



size_t WithChechikOracle::getSize(){
    size_t s1 = OGS->getSize(),
           s2 = RCO->getSize();
    return std::max(s1, s2);
}