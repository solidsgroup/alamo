
INCLUDE= -I./PFAmr/ 
COMPILE_FLAGS = -std=c++11 -ggdb -DDIM=2 -BL_SPACEDIM=2
LINK_FLAGS = -lamrex -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran 

HDR = PFAmr/PFAmr.H
SRC = main.cpp PFAmr/PFAmr.cpp PFAmr/PFAmrError.cpp PFAmr/PFAmrEvolve.cpp PFAmr/PFAmrInit.cpp PFAmr/PFAmrIO.cpp 
OBJ = ${SRC:.cpp=.o}

alamo: ${OBJ} ${HDR}
	mpicxx ${OBJ} -std=c++11 ${LIB} ${LINK_FLAGS} -o alamo
%.o: %.cpp ${HDR}
	mpicxx -c $< -o $@ ${INCLUDE} ${COMPILE_FLAGS}

clean:
	find . -name "*.o" -exec rm {} \;
	rm -f Backtrace*
