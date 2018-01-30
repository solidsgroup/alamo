RESET              = '\033[0m'
B_ON               = '\033[1m'
FG_RED             = '\033[31m'
FG_GREEN           = '\033[32m'
FG_YELLOW          = '\033[33m'
FG_BLUE            = '\033[34m'

DIM = 2

MPICXX_COMPILE_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro 
MPIFORT_COMPILE_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro 
#MPI_LINK_FLAGS =  -Wl,-Bsymbolic-functions -Wl,-z,relro 

CXX_COMPILE_FLAGS = -Wpedantic -Wextra -Wall  -std=c++11 -DDIM=${DIM} -DBL_SPACEDIM=${DIM} 

INCLUDE = -I./src/ -I./src/PFAmr/ -I./src/PFFem/ -I./src/PFFlame/ -I./src/PFBoundary/ -I/opt/amrex/tmp_install_dir/include/
LIB     = -lamrex -lgfortran -lmpichfort -lmpich  



HDR = $(shell find src/ -name *.H)
SRC = $(shell find src/ -mindepth 2  -name "*.cpp" )
SRC_F = $(shell find src/ -mindepth 2  -name "*.F90" )
OBJ = ${SRC:.cpp=.cpp.o}
OBJ_F = ${SRC_F:.F90=.F90.o}

#alamo:bin/alamo

#fem:bin/fem

default: bin/alamo bin/fem bin/flame

bin/alamo: ${OBJ} ${OBJ_F} src/main.cpp.o
	@echo $(B_ON)$(FG_BLUE)"###"
	@echo "### LINKING $@" 
	@echo "###"$(RESET)
	mkdir -p bin/
	mpicxx -o bin/alamo $^ ${LIB} 

bin/fem: ${OBJ} ${OBJ_F} src/fem.cpp.o
	@echo $(B_ON)$(FG_BLUE)"###"
	@echo "### LINKING $@" 
	@echo "###"$(RESET)
	mkdir -p bin/
	mpicxx -o bin/fem $^ ${LIB} 

bin/flame: ${OBJ} ${OBJ_F} src/flame.cpp.o
	@echo $(B_ON)$(FG_BLUE)"###"
	@echo "### LINKING $@" 
	@echo "###"$(RESET)
	mkdir -p bin/
	mpicxx -o bin/flame $^ ${LIB} 

%.cpp.o: %.cpp ${HDR}
	@echo $(B_ON)$(FG_YELLOW)"###"
	@echo "### COMPILING $<" 
	@echo "###"$(RESET)
	mpicxx -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} ${MPICXX_COMPILE_FLAGS}

%.F90.o: %.F90 ${HDR}
	@echo $(B_ON)$(FG_YELLOW)"###"
	@echo "### COMPILING $<" 
	@echo "###"$(RESET)
	mpif90 -c $< -o $@ ${INCLUDE} 

clean:
	find src/ -name "*.o" -exec rm {} \;
	rm -f Backtrace*
