DIM = 2
CC = g++

MPICXX_COMPILE_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro -I/usr/include/mpich -L/usr/lib/x86_64-linux-gnu -lmpichcxx -lmpich
MPICXX_LINK_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro -I/usr/include/mpich -L/usr/lib/x86_64-linux-gnu -lmpichcxx -lmpich
MPIFORT_LINK_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro -I/usr/include/mpich -I/usr/include/mpich -L/usr/lib/x86_64-linux-gnu -lmpichfort -lmpich

INCLUDE= -I./src/ -I./src/PFAmr/ -I./src/PFBoundary/ 
CXX_COMPILE_FLAGS = -std=c++11 -DDIM=${DIM} -DBL_SPACEDIM=${DIM}
CXX_LINK_FLAGS = -lamrex -std=c++11
#-lgfortran
#-lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh 


HDR = $(shell find src/ -name *.H)
SRC = $(shell find src/ -mindepth 2  -name "*.cpp" )
OBJ = ${SRC:.cpp=.o}

alamo: ${OBJ} src/main.o
	mpicxx -o bin/alamo $^ ${LIB} ${CXX_LINK_FLAGS} ${MPICXX_LINK_FLAGS} -lgfortran ${MPIFORT_LINK_FLAGS}
fem: ${OBJ} src/fem.o
	mpicxx -o bin/fem $^ ${LIB} ${CXX_LINK_FLAGS} ${MPICXX_LINK_FLAGS} -lgfortran ${MPIFORT_LINK_FLAGS}
%.o: %.cpp ${HDR}
	mpicxx -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} ${MPICXX_COMPILE_FLAGS}

clean:
	find src/ -name "*.o" -exec rm {} \;
	rm -f Backtrace*
