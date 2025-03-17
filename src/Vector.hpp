#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cstring>
#include <initializer_list>

#include "Utilities.hpp"

using std::cout, std::endl;

typedef unsigned int uint;

template<typename T>
class Vector{
    public:
    Vector<T>() : length(0), allocated(1) {
        head = new T[1];
    }

    Vector<T>(uint space) : length(0), allocated(space) {
        head = new T[space];
    }

    Vector<T>(uint quantity, const T& object) : length(quantity), allocated(quantity) {
        head = new T[quantity];
        for(int i=0; i<quantity; ++i)
            head[i] = object;
    }

    Vector<T>(const Vector<T>& rhs) : length(rhs.length), allocated(rhs.allocated) {
        head = new T[rhs.allocated];
        for(int i=0; i<rhs.length; ++i)
            head[i] = rhs[i];
    }

    Vector<T>(const std::initializer_list<T>& list) : length(0), allocated(list.size()) {
        head = new T[list.size()];
        for(const T k : list)
            pushBack(k);
    }

    inline uint size() const{
        return length;
    }
    inline uint capacity() const{
        return allocated;
    }
    inline uint empty() const{
        return !length;
    }
    
    void reserve(const uint newSize);
    void resize(const uint newSize);
    void pushBack(const T& object);
    void popBack();
    
    inline T& operator[](const uint position) const{
        return head[position];
    }

    Vector<T>& operator=(const Vector<T>& rhs){
        delete[] head;
        length = rhs.length;
        allocated = rhs.allocated;

        head = new T[rhs.allocated];
        for(int i=0; i<rhs.length; ++i)
            head[i] = rhs[i];
        return *this;
    }

    template<typename U>
    friend std::ostream& operator<<(std::ostream& out, const Vector<U>& v);

    ~Vector<T>(){
        delete[] head;
    }

    struct iterator{
    public:
        iterator(T* position) : mPtr(position) {}

        iterator& operator++(){
            ++mPtr;
            return *this;
        }
        iterator operator++(int){
            iterator result(mPtr);
            ++mPtr;
            return result;
        }
        T& operator*(){
            return *mPtr;
        }
        T* operator->(){
            return &(*mPtr);
        }
        bool operator!=(const iterator& it) const{
            return mPtr != it.mPtr;
        }
    private:
        T* mPtr;
    };
    iterator begin() { return iterator(head); }
    iterator end() { return iterator(head+length); }

    struct const_iterator{
    public:
        const_iterator(T* position) : mPtr(position) {}

        const_iterator& operator++(){
            ++mPtr;
            return *this;
        }
        const_iterator operator++(int){
            const_iterator result(mPtr);
            ++mPtr;
            return result;
        }
        T operator*() const{
            return *mPtr;
        }
        bool operator!=(const const_iterator& it) const{
            return mPtr != it.mPtr;
        }
    private:
        T* mPtr;
    };
    const_iterator begin() const { return const_iterator(head); }
    const_iterator end() const { return const_iterator(head+length); }

    private:
    T* head;
    uint length;
    uint allocated;
};

template <typename T>
void Vector<T>::reserve(const uint newSize){
    T* newHead = new T[newSize];
    
    if(length)
        std::memcpy(newHead, head, sizeof(T)*length);

    allocated = max(allocated, newSize);

    delete[] head;
    head = newHead;
    
    return;
}

template <typename T>
void Vector<T>::resize(const uint newSize){
    T* newHead = new T[newSize];
    
    if(length)
        std::memcpy(newHead, head, sizeof(T)*min(newSize, length));

    allocated = max(allocated, newSize);
    length = allocated;

    delete[] head;
    head = newHead;
    
    return;
}

template<typename T>
void Vector<T>::pushBack(const T& object){
    if(length == allocated)
        reserve(2*allocated);

    head[length] = object;
    ++length;
    
    return;
}

template<typename T>
void Vector<T>::popBack(){
    --length;            
    return;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const Vector<T>& v){
    for(const T& k : v)
        out << k << ' ';
    return out;
}


#endif
