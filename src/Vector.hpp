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
    Vector() : length(0), allocated(1) {
        head = new T[1];
    }

    Vector(uint space, uint length=0) : length(length), allocated(space) {
        head = new T[space];
    }

    Vector(uint space, uint length, const T& object) : length(length), allocated(space) {
        head = new T[space];
        for(int i=0; i<length; ++i)
            head[i] = object;
    }

    Vector(const Vector<T>& rhs) : length(rhs.length), allocated(rhs.allocated) {
        head = new T[rhs.allocated];
        for(uint i=0; i<rhs.length; ++i)
            head[i] = rhs[i];
    }

    Vector(const std::initializer_list<T>& list) : length(0), allocated(list.size()) {
        head = new T[list.size()];
        for(const T k : list)
            pushBack(k);
    }

    inline uint capacity() const { return allocated; }
    inline const T* data() const { return head; }
    inline uint empty() const { return !length; }
    inline uint size() const { return length; }
    
    void clear();
    void copyBack(const T* const dataPointer, uint numberOfObjectsToCopy);
    void copyBack(const Vector<T>& dataVector, uint begin, uint end);
    Vector<T> slice(uint begin, uint end);
    void popBack();
    void pushBack(const T& object);
    void reserve(uint newSize);
    void resize(uint newSize);

    template<typename U>
    friend bool operator==(const Vector<U>& left, const Vector<U> right);
    
    inline T& operator[](uint position) const{
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

    ~Vector(){
        delete[] head;
        head = nullptr;
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


template<typename T>
void Vector<T>::clear(){
    delete[] head;
    head = new T[1];
    length = 0;
    allocated = 1;
}

template<typename T>
void Vector<T>::copyBack(const T* const dataPointer, uint numberOfObjectsToCopy){
    if(length + numberOfObjectsToCopy > allocated)
        reserve(length + numberOfObjectsToCopy);

    // std::memcpy(head + length, dataPointer, numberOfObjectsToCopy*sizeof(T));
    for(uint i=0; i<numberOfObjectsToCopy; ++i)
        head[length+i] = dataPointer[i];
    length += numberOfObjectsToCopy;
    return;
}

template<typename T>
void Vector<T>::copyBack(const Vector<T>& dataVector, uint begin, uint end){
    if(length + end - begin > allocated)
        reserve(length + end - begin);

    // std::memcpy(head + length, dxataVector.data() + begin, (end-begin)*sizeof(T));
    for(uint i=begin; i<end; ++i)
        head[length+i-begin] = dataVector[i];
    length += end - begin;
    return;
}

template<typename T>
void Vector<T>::popBack(){
    --length;            
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

template <typename T>
void Vector<T>::reserve(uint newSize){
    T* newHead = new T[newSize];
    
    // std::memcpy(newHead, head, length*sizeof(T));
    for(uint i=0; i<length; ++i){
        newHead[i] = head[i];
    }

    allocated = max(allocated, newSize);

    delete[] head;
    head = newHead;
    
    return;
}

template <typename T>
void Vector<T>::resize(uint newSize){
    T* newHead = new T[newSize];
    
    // std::memcpy(newHead, head, min(newSize, length)*sizeof(T));
    for(uint i=0; i<min(length,newSize); ++i){
        newHead[i] = head[i];
    }

    allocated = newSize;
    length = newSize;

    delete[] head;
    head = newHead;
    
    return;
}

template<typename T>
Vector<T> Vector<T>::slice(uint begin, uint end){
    uint size = end - begin;
    Vector<T> subVec(size);
    subVec.copyBack(*this, begin, end);
    return subVec;
}

template<typename T>
bool operator==(const Vector<T>& left, const Vector<T> right){
    if(left.length != right.length)
        return false;
    for(uint i=0; i<left.length; ++i){
        if(left[i] != right[i])
            return false;
    }
    return true;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const Vector<T>& v){
    for(const T& k : v)
        out << k << ' ';
    return out;
}


#endif
