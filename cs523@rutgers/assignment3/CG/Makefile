ROOT_PATH=/ilab/users/hs806/cpp_workspace/CG
LD_INCLUDE=$(ROOT_PATH)/include 
LD_SRC=$(ROOT_PATH)/src 


CXX = g++ -std=c++11 -DCOMPILE_WITH_DOUBLE_SUPPORT -I$(LD_INCLUDE)

all: main

Krylov_System_Base.o: 
	$(CXX) ./src/Krylov_System_Base.cpp  -c

Conjugate_Gradient.o: 
	$(CXX) ./src/Conjugate_Gradient.cpp -c

main: Krylov_System_Base.o Conjugate_Gradient.o
	$(CXX) ./Krylov_System_Base.o ./Conjugate_Gradient.o ./src/main.cpp -o main

clean:
	rm -f *.o main