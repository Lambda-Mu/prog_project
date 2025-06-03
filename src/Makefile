OBJ = main.o Homology.o LefschetzComplex.o CubicalComplex.o Matrix.o Vector.o Utilities.o
CXX = g++
FLAGS = -std=c++20 -O2

main: $(OBJ)
	$(CXX) $(FLAGS) $(OBJ) -o main

main.o: main.cpp LefschetzComplex.hpp CubicalComplex.hpp Matrix.hpp Vector.hpp Utilities.hpp
	$(CXX) $(FLAGS) -c main.cpp

Homology.o: Homology.hpp LefschetzComplex.hpp Matrix.hpp Vector.hpp
	$(CXX) $(FLAGS) -c Homology.cpp

CubicalComplex.o: CubicalComplex.hpp CubicalComplex.cpp Vector.hpp LefschetzComplex.hpp Matrix.hpp
	$(CXX) $(FLAGS) -c CubicalComplex.cpp

LefschetzComplex.o: LefschetzComplex.hpp LefschetzComplex.cpp Matrix.hpp
	$(CXX) $(FLAGS) -c LefschetzComplex.cpp

Matrix.o: Matrix.hpp Matrix.cpp Vector.hpp Utilities.hpp
	$(CXX) $(FLAGS) -c Matrix.cpp

Vector.o: Vector.hpp Vector.cpp Utilities.hpp
	$(CXX) $(FLAGS) -c Vector.cpp

Utilities.o: Utilities.hpp Utilities.cpp
	$(CXX) $(FLAGS) -c Utilities.cpp

clean:
	rm *.o main
