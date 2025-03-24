#ifndef UTILITIES_H
#define UTILITIES_H

template<typename T>
const T& max(const T& a, const T& b){
    if(a >= b)
        return a;
    return b; 
}

template <typename T>
T& max(T& a, T& b){
    if(a >= b)
        return a;
    return b;
}

template<typename T>
const T& min(const T& a, const T& b){
    if(a <= b)
        return a;
    return b; 
}

template<typename T>
T& min(T& a, T& b){
    if(a <= b)
        return a;
    return b;
}

template<typename T>
void swap(T& a, T& b){
    T tmp(a);
    a = b;
    b = tmp;
    return;
}

template<typename T>
T floor(const T& dividend, const T& divisor){
    if(dividend < 0 xor divisor < 0)
        return dividend / divisor - 1;
    return dividend / divisor;
}

#endif