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

    Vector<T>(const uint space, const uint length=0) : length(length), allocated(space) {
        head = new T[space];
    }

    Vector<T>(const uint space, const uint length, const T& object) : length(length), allocated(space) {
        head = new T[space];
        for(int i=0; i<length; ++i)
            head[i] = object;
    }

    Vector<T>(const Vector<T>& rhs) : length(rhs.length), allocated(rhs.allocated) {
        head = new T[rhs.allocated];
        std::memcpy(head, rhs.head, rhs.length*sizeof(T));
    }

    Vector<T>(const std::initializer_list<T>& list) : length(0), allocated(list.size()) {
        head = new T[list.size()];
        for(const T k : list)
            pushBack(k);
    }

    inline uint capacity() const { return allocated; }
    inline const T* data() const { return head; }
    inline uint empty() const { return !length; }
    inline uint size() const { return length; }
    
    void copyBack(const T* const dataPointer, const uint numberOfObjectsToCopy);
    void copyBack(const Vector<T>& dataVector, const uint begin, const uint end);
    Vector<T> slice(const uint begin, const uint end);
    void popBack();
    void pushBack(const T& object);
    void reserve(const uint newSize);
    void resize(const uint newSize);
    
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


template<typename T>
void Vector<T>::copyBack(const T* const dataPointer, const uint numberOfObjectsToCopy){
    if(length + numberOfObjectsToCopy > allocated)
        reserve(length + numberOfObjectsToCopy);

    std::memcpy(head + length, dataPointer, numberOfObjectsToCopy*sizeof(T));
    length += numberOfObjectsToCopy;
    return;
}

template<typename T>
void Vector<T>::copyBack(const Vector<T>& dataVector, const uint begin, const uint end){
    if(length + end - begin > allocated)
        reserve(length + end - begin);

    std::memcpy(head + length, dataVector.data() + begin, (end-begin)*sizeof(T));
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
void Vector<T>::reserve(const uint newSize){
    T* newHead = new T[newSize];
    
    std::memcpy(newHead, head, length*sizeof(T));

    allocated = max(allocated, newSize);

    delete[] head;
    head = newHead;
    
    return;
}

template <typename T>
void Vector<T>::resize(const uint newSize){
    T* newHead = new T[newSize];
    
    std::memcpy(newHead, head, min(newSize, length)*sizeof(T));

    allocated = newSize;
    length = newSize;

    delete[] head;
    head = newHead;
    
    return;
}

template<typename T>
Vector<T> Vector<T>::slice(const uint begin, const uint end){
    uint size = end - begin;
    Vector<T> subVec(size);
    std::memcpy(subVec.head, head, size*sizeof(T));
    subVec.length = size;
    return subVec;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const Vector<T>& v){
    for(const T& k : v)
        out << k << ' ';
    return out;
}


#endif
