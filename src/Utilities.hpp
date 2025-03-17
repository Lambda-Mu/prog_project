#ifndef UTILITIES_H
#define UTILITIES_H

template<typename T>
const T& max(const T& a, const T& b){
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

#endif