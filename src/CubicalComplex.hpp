#ifndef CUBICAL_COMPLEX_H
#define CUBICAL_COMPLEX_H

#include <string>
#include <charconv>
#include <iostream>
#include <set>

#include "Vector.hpp"

using std::string, std::cout;

class Block{
    public:
        Block() {}
        Block(int leftEndpoint, int rightEndpoint) 
            : leftEndpoint(leftEndpoint), rightEndpoint(rightEndpoint) { }
        Block(const string& str){
            int delimiter = str.find(',');
            if(delimiter == string::npos){
                int endpoint;
                std::from_chars(str.data()+1, str.data()+str.size(), endpoint);
                leftEndpoint = endpoint;
                rightEndpoint = endpoint;
            }else{
                int endpoint;
                std::from_chars(str.data()+1, str.data()+delimiter, endpoint);
                leftEndpoint = endpoint;
                rightEndpoint = endpoint+1;
            }
        }

        inline int& leftEnd() { return leftEndpoint; }
        inline int leftEnd() const { return leftEndpoint; }
        inline int& rightEnd() { return rightEndpoint; }
        inline int rightEnd() const { return rightEndpoint; }
        inline bool nonDegenerate() const { return leftEndpoint != rightEndpoint; }

        friend std::ostream& operator<<(std::ostream& out, const Block& block);
        friend bool operator==(const Block& left, const Block& right);
        inline friend bool operator!=(const Block& left, const Block& right) { return !(left==right); }
        friend bool operator<(const Block& left, const Block& right);
        inline friend bool operator>(const Block& left, const Block& right) { return right<left; }
        inline friend bool operator<=(const Block& left, const Block& right) { return (left<right) or (left==right); }
        inline friend bool operator>=(const Block& left, const Block& right) { return (left>right) or (left==right); }

    private:
        int leftEndpoint;
        int rightEndpoint;
    };

class Cube{
public:
    Cube() {}
    Cube(const Vector<Block>& buildingBlocks)
        : buildingBlocks(buildingBlocks) {
        uint nonDegenerateBlocks = 0;
        for(uint i=0; i<dimension(); ++i){
            if((*this)[i].nonDegenerate())
                ++nonDegenerateBlocks;
        }
        embeddingNumber = nonDegenerateBlocks;
    }
    Cube(const string& str){
        uint beginSearch = 0;
        int foundIndex = str.find('x', beginSearch);
        while(foundIndex != string::npos){
            buildingBlocks.pushBack(Block(str.substr(beginSearch, foundIndex-beginSearch)));
            beginSearch = foundIndex+1;
            foundIndex = str.find('x', beginSearch);
        }
        buildingBlocks.pushBack(str.substr(beginSearch, str.size()));
        embeddingNumber = 0;
        for(uint i=0; i<buildingBlocks.size(); ++i){
            if(buildingBlocks[i].nonDegenerate())
                ++embeddingNumber;
        }
    }

    inline uint embeddingDimension() const { return embeddingNumber; }
    inline uint dimension() const { return buildingBlocks.size(); }
    inline Block& operator[](uint position) { return buildingBlocks[position]; }
    inline Block& operator[](uint position) const { return buildingBlocks[position]; }
    Vector<Cube> getPrimaryFaces() const;

    friend std::ostream& operator<<(std::ostream& out, const Cube& cube);
    friend bool operator==(const Cube& left, const Cube& right);
    inline friend bool operator!=(const Cube& left, const Cube& right) { return !(left==right); }
    friend bool operator<(const Cube& left, const Cube& right);
    inline friend bool operator>(const Cube& left, const Cube& right) { return right<left; }
    inline friend bool operator<=(const Cube& left, const Cube& right) { return (left<right) or (left==right); }
    inline friend bool operator>=(const Cube& left, const Cube& right) { return (left>right) or (left==right); }

private:
    uint embeddingNumber;
    Vector<Block> buildingBlocks;
};

class CubicalSet{
public:
    CubicalSet() : cubes(Vector<Cube>()) { }
    CubicalSet(Vector<Cube> cubes) : cubes(cubes) { }

    Cube& operator[](uint position) { return cubes[position]; }
    Cube& operator[](uint position) const { return cubes[position]; }
    void clear() { cubes.clear(); }
    void add(const string& cube) { cubes.pushBack(Cube(cube)); }
    inline bool empty() const { return cubes.empty(); }
    inline uint dimension() const { return cubes[0].dimension(); }
    inline uint numberOfCubes() const { return cubes.size(); }
    static void getFromFile(const string& filename, CubicalSet& cubicalSet);
    std::set<Cube> getAllFaces() const;

    friend std::ostream& operator<<(std::ostream& out, const CubicalSet& cubicalSet);

private:
    Vector<Cube> cubes;
};



#endif