#ifndef UTILITIES_H
#define UTILITIES_H

typedef unsigned int uint;

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
    if(dividend == 0)
        return dividend;
    if(dividend < 0 xor divisor < 0 and dividend % divisor != 0)
        return dividend / divisor - 1;
    return dividend / divisor;
}

inline bool floor(const bool dividend, const bool divisor){
    return dividend;
}

template<typename T>
T absoluteValue(const T a){
    if(a >= 0)
        return a;
    return -a;
}

// inline bool abs(const bool value){
//     return value;
// }

template<typename T>
struct OneIndexedValue{
public:
    OneIndexedValue(uint index, T value)
        : index(index), value(value) { }
    
    uint index;
    T value;
};

template<typename T>
struct DoubleIndexedValue{
public:    
    DoubleIndexedValue(uint rowIndex, uint columnIndex, T value)
        : rowIndex(rowIndex), columnIndex(columnIndex), value(value) { }

    uint rowIndex;
    uint columnIndex;
    T value;
};

#endif