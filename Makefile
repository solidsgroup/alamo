
INCLUDE=-I/opt/amrex/tmp_install_dir/include/ -I./AmrAdv/ -I/usr/include/vtk-6.2
LIB = /opt/amrex/tmp_install_dir/lib/libamrex.a 
COMPILE_FLAGS = -std=c++11
LINK_FLAGS = -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -L/opt/Trelis-15.1/bin/ 

SRC = main.cpp AmrAdv/AmrAdv.cpp AmrAdv/AmrAdvError.cpp AmrAdv/AmrAdvEvolve.cpp AmrAdv/AmrAdvInit.cpp AmrAdv/AmrAdvIO.cpp 
OBJ = ${SRC:.cpp=.o}

a.out: ${OBJ}
	mpicxx ${OBJ} -std=c++11 ${LIB} ${LINK_FLAGS}
%.o: %.cpp
	mpicxx -c $< -o $@ ${INCLUDE} ${COMPILE_FLAGS}

clean:
	find . -name "*.o" -exec rm {} \;
	rm -f Backtrace*
