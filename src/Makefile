OBJ = main.o SpecialMatrixForms.o Matrix.o Vector.o Utilities.o
CXX = g++
FLAGS = -std=c++17 -O1

main: $(OBJ)
	$(CXX) $(FLAGS) $(OBJ) -o main

main.o: main.cpp Matrix.hpp Vector.hpp Utilities.hpp
	$(CXX) $(FLAGS) -c main.cpp

SpecialMatrixForms.o: SpecialMatrixForms.hpp SpecialMatrixForms.cpp Matrix.hpp
	$(CXX) $(FLAGS) -c SpecialMatrixForms.cpp

Matrix.o: Matrix.hpp Matrix.cpp Vector.hpp Utilities.hpp
	$(CXX) $(FLAGS) -c Matrix.cpp

Vector.o: Vector.hpp Vector.cpp Utilities.hpp
	$(CXX) $(FLAGS) -c Vector.cpp

Utilities.o: Utilities.hpp Utilities.cpp
	$(CXX) $(FLAGS) -c Utilities.cpp

clean:
	rm *.o main