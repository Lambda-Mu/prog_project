#include "LefschetzComplex.hpp"
#include "CubicalComplex.hpp"
#include "Vector.hpp"

#include <queue>
#include <set>
#include <map>

LefschetzComplexZ2<Cube> getLefschetzComplexZ2(const CubicalSet& skeleton){
    std::queue<Cube> cubes;

    uint dim = skeleton.dimension();
    Vector<std::set<Cube>> cells(dim);
    Vector<std::multimap<Cube, Cube>> incidence(dim);

    for(uint i=0; i<skeleton.numberOfCubes(); ++i){
        cubes.push(skeleton[i]);
        cells[skeleton[i].embeddingDimension()].insert(skeleton[i]);
    }

    while(!cubes.empty()){
        Cube cube = cubes.front();
        Vector<Cube> faces = cube.getPrimaryFaces();
        for(Cube& q : faces){
            cubes.push(q);
            incidence[cube.dimension()].insert({cube, q});
        }
        cubes.pop();
    }

    return LefschetzComplexZ2<Cube>(cells, incidence);
}