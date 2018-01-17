DIM = 2
CC = g++

MPICXX_COMPILE_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro -I/usr/include/mpich -L/usr/lib/x86_64-linux-gnu -lmpichcxx -lmpich
MPICXX_LINK_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro -I/usr/include/mpich -L/usr/lib/x86_64-linux-gnu -lmpichcxx -lmpich
MPIFORT_LINK_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro -I/usr/include/mpich -I/usr/include/mpich -L/usr/lib/x86_64-linux-gnu -lmpichfort -lmpich

INCLUDE= -I./src/PFAmr/ -I./src/PFBoundary/ 
CXX_COMPILE_FLAGS = -std=c++11 -DDIM=${DIM} -DBL_SPACEDIM=${DIM}
CXX_LINK_FLAGS = -lamrex 
#-lgfortran
#-lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh 


HDR = src/PFAmr/PFAmr.H src/PFBoundary/PFBoundary.H src/PFBoundary/PFBoundarySin.H
SRC = src/main.cpp src/PFAmr/PFAmr.cpp src/PFAmr/PFAmrError.cpp src/PFAmr/PFAmrEvolve.cpp src/PFAmr/PFAmrInit.cpp src/PFAmr/PFAmrIO.cpp 
OBJ = ${SRC:.cpp=.o}

alamo: ${OBJ}
	mkdir bin
	mpicxx -o bin/alamo ${OBJ} ${LIB} ${CXX_LINK_FLAGS} ${MPICXX_LINK_FLAGS} -lgfortran ${MPIFORT_LINK_FLAGS}
%.o: %.cpp ${HDR}
	mpicxx -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} ${MPICXX_COMPILE_FLAGS}

clean:
	find src/ -name "*.o" -exec rm {} \;
	rm -f Backtrace*
