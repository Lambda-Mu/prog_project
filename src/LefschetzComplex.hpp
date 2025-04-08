#ifndef LEFSCHETZ_COMPLEX_H
#define LEFSCHETZ_COMPLEX_H

#include <map>
#include <set>

#include "Vector.hpp"

using std::multimap, std::set;

template<class Cell>
class LefschetzComplexZ2{
public:
    LefschetzComplexZ2(const Vector<set<Cell>>& cells, const Vector<multimap<Cell, Cell>>& incidenceMap)
        : cells(cells), incidenceMap(incidenceMap) { }

    bool incidence(const Cell& cube, const Cell& face) const;

private:
    Vector<set<Cell>> cells;
    Vector<multimap<Cell, Cell>> incidenceMap;
};

// template<class Cell>
// bool LefschetzComplexZ2<Cell>::incidence(const Cell& cube, const Cell& face) const{
    
// }


#endif
