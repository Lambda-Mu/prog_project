#include <fstream>
#include <queue>
#include <set>

#include "CubicalComplex.hpp"
#include "Vector.hpp"

std::ostream& operator<<(std::ostream& out, const Block& block){
    if(block.nonDegenerate())
        out << '[' << block.leftEnd() << ',' << block.rightEnd() << ']';
    else
        out << '[' << block.leftEnd() << ']';
    return out;
}

std::ostream& operator<<(std::ostream& out, const Cube& cube){
    for(uint i=0; i<cube.embeddingDimension()-1; ++i){
        out << cube[i] << 'x';
    }
    out << cube[cube.embeddingDimension()-1];
    return out;
}

std::ostream& operator<<(std::ostream& out, const CubicalSet& cubicalSet){
    for(uint i=0; i<cubicalSet.numberOfCubes()-1; ++i){
        out << cubicalSet[i] << '\n';
    }
    out << cubicalSet[cubicalSet.numberOfCubes()-1];
    return out;
}

bool operator==(const Block& left, const Block& right){
    if(left.leftEnd() != right.leftEnd())
        return false;
    if(left.rightEnd() != right.rightEnd())
        return false;
    return true;
}

bool operator==(const Cube& left, const Cube& right){
    if(left.embeddingDimension() != right.embeddingDimension())
        return false;
    for(uint i=0; i<left.buildingBlocks.size(); ++i){
        if(left.buildingBlocks[i] != right.buildingBlocks[i])
            return false;
    }
    return true;
}

bool operator<(const Block& left, const Block& right){
    if(left.leftEnd() < right.leftEnd())
        return true;
    if(left.leftEnd() > right.rightEnd())
        return false;
    if(left.rightEnd() < right.rightEnd())
        return true;
    return false;
}

bool operator<(const Cube& left, const Cube& right){
    uint leftEmb = left.embeddingDimension();
    uint rightEmb = right.embeddingDimension();
    if(leftEmb < rightEmb)
        return true;
    if(rightEmb < leftEmb)
        return false;
    for(int i=0; i<leftEmb; ++i){
        if(left[i] < right[i])
            return true;
        if(left[i] > right[i])
            return false;
    }
    return false;
}

void CubicalSet::getFromFile(const string& filename, CubicalSet& cubicalSet){
    std::ifstream input(filename);
    if(!input.is_open()){
        std::cerr << "Couldn't open file: \"" << filename << "\"";
        std::cerr << "during execution of CubicalSet::getFromFile. Errors may occur.\n";
        return;
    }
    cubicalSet.clear();
    while(!input.eof()){
        string cube;
        std::getline(input, cube);
        cubicalSet.add(cube);
    }
    input.close();
    return;
}

Vector<Cube> Cube::getPrimaryFaces() const{
    Vector<Cube> primaryFaces(2*dimension());
    for(uint i=0; i<embeddingDimension(); ++i){
        if(!buildingBlocks[i].nonDegenerate())
            continue;
        Cube face(*this);
        --face.dim;
        face[i].rightEnd() = (*this)[i].leftEnd(); 
        primaryFaces.pushBack(face);
        face[i].leftEnd() = (*this)[i].rightEnd();
        face[i].rightEnd() = (*this)[i].rightEnd();
        primaryFaces.pushBack(face);
    }
    return primaryFaces;
}

std::set<Cube> CubicalSet::getAllFaces() const{
    std::queue<Cube> cubes;
    std::set<Cube> faces;
    for(uint i=0; i<numberOfCubes(); ++i){
        cubes.push((*this)[i]);
    }

    while(!cubes.empty()){
        Cube cube = cubes.front();
        faces.insert(cube);
        Vector<Cube> newFaces = cube.getPrimaryFaces();
        for(Cube& q : newFaces)
            cubes.push(q);
        cubes.pop();
    }
    return faces;
}

Vector<Matrix<bool>> CubicalComplexZ2::getMatrixOfBoundaryOperator() const{
    Vector<Matrix<bool>> bd;
    bd.reserve(incidence.size());
    bd.pushBack(Matrix<bool>(
        1, bases[0].size(), 
        Vector<bool>(bases[0].size(), bases[0].size(), 0))
    );
    for(uint i=0; i<incidence.size()-1; ++i){
        Vector<bool> vec(
            bases[i].size()*bases[i+1].size(), 
            bases[i].size()*bases[i+1].size(), false
        );
        bd.pushBack(Matrix<bool>(bases[i].size(), bases[i+1].size(), vec));
        for(uint k=0; k<bases[i+1].size(); ++k){
            for(uint j=0; j<bases[i].size(); ++j)
                if((*this)(bases[i][j], bases[i+1][k])){
                    bd[i+1](j, k) = true;
                }
        }
    }
    return bd;
}