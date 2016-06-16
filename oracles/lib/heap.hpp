#include <stddef.h>
#include <utility>
#include <vector>
#include <cmath>
#include <float.h>
#include <functional>
#include <iostream>

/*
A min/max (op) binary heap data structure.
Requires template type T to have a
ptr
attribute, to reference local nodes.


TODO change T* (pointers) to T& (references) as null pointers are not handled,
and references are prettier, see http://www.cplusplus.com/articles/z6vU7k9E/

Note that template classes cannot be hidden behind headers in any nice way.
http://stackoverflow.com/questions/1724036/splitting-templated-c-classes-into-hpp-cpp-files-is-it-possible
*/


#ifndef HEAP_H
#define HEAP_H

//Set T as address to Vertex.
template <typename T> 
struct Node {
    T* handle; //cannot be const when using vector::swap!
    double key; //value aka distance
    size_t index; //in heap
    Node(T* k, double v, size_t i): handle(k),key(v),index(i) {}
    ~Node(){
        if(handle)
            handle->ptr = reinterpret_cast<uintptr_t>(nullptr);//nullptr;
    }
};

template <typename T>
class Heap {
private:
    std::vector<Node<T>* > heap;
    std::function<bool(double,double)> Op;    
    size_t parent(size_t index){
        return std::floor(index/2);
    }

    size_t left(size_t index){
        return 2*index;
    }

    size_t right(size_t index){
        return 2*index + 1;
    }

    //heapify heap-vector.
    void heapify(size_t i){
        size_t l = left(i);
        size_t r = right(i);
        //t is either largest or smallest -> depending on Op(erator).
        int t = l <= heap.size()-1 && Op(heap[l]->key, heap[i]->key) ? l : i;
        if(r <= heap.size()-1 && Op(heap[r]->key, heap[t]->key) )
            t = r;
        while(t != i){
            std::swap(heap[i],heap[t]); //swap(i,t);
            heap[i]->index = i;
            heap[t]->index = t;
            i = t;
            //recursion step

            l = left(i);
            r = right(i);
            t = l <= heap.size()-1 && Op(heap[l]->key, heap[i]->key) ? l : i;
            if(r <= heap.size()-1 && Op(heap[r]->key, heap[t]->key))
                t = r;
        }
    }

    void balance(){
        if(heap.size() > 1)
            heapify(0);
    }

    //push node to back.
    //no balancing done
    void push(T* handle, double key){
        Node<T>* n = new Node<T>(handle, key, heap.size());
        handle->ptr = reinterpret_cast<uintptr_t>(n);
        heap.push_back(n);
    }


    //no balancing done
    void deleteHead(){
        if(heap.size() < 1){
            //todo throw error
            return;
        }
        Node<T>* n = heap[0];
        heap[0] = heap[heap.size() - 1];
        heap.pop_back();
        
        if (heap.size() > 0) {
            heap[0]->index = 0;
        }
        delete n;
    }

    //assumes input validation already done
    size_t updateKey(size_t i, double key){
        heap[i]->key = key;
        size_t p = parent(i);
        while(i > 0 && Op(heap[i]->key, heap[p]->key)){
            std::swap(heap[i],heap[p]); //swap(i,p);
            heap[i]->index = i;
            heap[p]->index = p;
            i = p;
            p = parent(i);
        }
        return i;
    }

public:
    Heap(std::function<bool(double,double)> op): Op(op) {}
    ~Heap(){
        for(auto n : heap){
            delete n;
        }
    }
    
    //input pair <handle, key>.
    std::vector<Node<T>* > const * 
    buildHeap(std::vector<std::pair<T*, double> >& inVec){
        //build local data struct
        auto itBegin = inVec.begin();
        for(auto it = itBegin; it != inVec.end(); ++it){
            push(it->first, it->second);
        }
        for(int i=std::floor(heap.size()/2); i >= 0; --i){ //size_t = unsigned; i>= 0 always true
            heapify(i);
        }

        return &heap;
    }

    
    /*************************** inspectors ***************************/
    size_t size(){
        return heap.size();
    }
    bool empty(){
        return heap.size() == 0;
    }

    //get readonly heap
    std::vector<Node<T>* > const * getHeap(){
        return &heap;
    }

    std::pair<T&,double> inspectHead(){
        Node<T>* h = heap[0];
        return {*(h->handle), h->key};
    }

    /*************************** manipulators ***************************/
    void removeHead(){
        deleteHead();
        balance();
    }
    
    std::pair<T&, double> pop(){
        std::pair<T&, double> h = inspectHead();
        removeHead();
        return h;
    }

    //set key if Op(new, old) is true
    //increase/decrease key
    //returns new index
    size_t alterKey(size_t index, double key){
        if (index >= heap.size())
            return 0;
        if(Op(key,heap[index]->key))
            return updateKey(index, key);
    } 

    //insert new key.
    //returns index
    size_t insert(T* handle, double key){
        if(Op(DBL_MAX, DBL_MIN))
            push(handle, DBL_MAX);
        else
            push(handle, DBL_MIN);
        return updateKey(heap.size()-1, key);
    }

    //pop head then insert new handle -- faster than delete-insert, 
    //but cannot return index.
    void replace(T* handle, double key){
        push(handle, key);
        removeHead();
    }
};

#endif /* HEAP_H */