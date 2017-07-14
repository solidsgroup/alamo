DIM = 2
CC = g++

MPICXX_COMPILE_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro -I/usr/include/mpich -L/usr/lib/x86_64-linux-gnu -lmpichcxx -lmpich
MPICXX_LINK_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro -I/usr/include/mpich -L/usr/lib/x86_64-linux-gnu -lmpichcxx -lmpich
MPIFORT_LINK_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro -I/usr/include/mpich -I/usr/include/mpich -L/usr/lib/x86_64-linux-gnu -lmpichfort -lmpich

INCLUDE= -I./PFAmr/ -I./PFBoundary/ 
CXX_COMPILE_FLAGS = -std=c++11 -DDIM=${DIM} -DBL_SPACEDIM=${DIM}
CXX_LINK_FLAGS = -lamrex 
#-lgfortran
#-lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh 


HDR = PFAmr/PFAmr.H PFBoundary/PFBoundary.H PFBoundary/PFBoundarySin.H
SRC = main.cpp PFAmr/PFAmr.cpp PFAmr/PFAmrError.cpp PFAmr/PFAmrEvolve.cpp PFAmr/PFAmrInit.cpp PFAmr/PFAmrIO.cpp 
OBJ = ${SRC:.cpp=.o}

alamo: ${OBJ}
	mpicxx -o alamo ${OBJ} ${LIB} ${CXX_LINK_FLAGS} ${MPICXX_LINK_FLAGS} -lgfortran ${MPIFORT_LINK_FLAGS}
%.o: %.cpp ${HDR}
	${CC} -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} ${MPICXX_COMPILE_FLAGS}

clean:
	find . -name "*.o" -exec rm {} \;
	rm -f Backtrace*
