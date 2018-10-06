CC = mpicxx

MPICHFORT ?= mpichfort

COMP ?= GCC
ifeq ($(COMP),INTEL)
 MPI_LIB = -lifcore
else ifeq ($(COMP),GCC)
 MPI_LIB = -lgfortran -l$(MPICHFORT) -lmpich
endif

RESET              = \033[0m
B_ON               = \033[1m
FG_RED             = \033[31m
FG_LIGHTRED        = \033[91m
FG_GREEN           = \033[32m
FG_LIGHTGREEN      = \033[92m
FG_YELLOW          = \033[33m
FG_BLUE            = \033[34m
FG_LIGHTBLUE       = \033[94m
FG_CYAN            = \033[36m



MPICXX_COMPILE_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro
MPIFORT_COMPILE_FLAGS = -Wl,-Bsymbolic-functions -Wl,-z,relro

METADATA_GITHASH  = $(shell git log --pretty=format:'%H' -n 1)
METADATA_USER     = $(shell whoami)
METADATA_PLATFORM = $(shell hostname)
METADATA_COMPILER = $(CC)
METADATA_DATE     = $(shell date +%x)
METADATA_TIME     = $(shell date +%H:%M:%S)


METADATA_FLAGS = -DMETADATA_GITHASH=\"$(METADATA_GITHASH)\" -DMETADATA_USER=\"$(METADATA_USER)\" -DMETADATA_PLATFORM=\"$(METADATA_PLATFORM)\" -DMETADATA_COMPILER=\"$(METADATA_COMPILER)\" -DMETADATA_DATE=\"$(METADATA_DATE)\" -DMETADATA_TIME=\"$(METADATA_TIME)\" 

CXX_COMPILE_FLAGS = -Wpedantic -Wextra -Wall  -std=c++11 $(METADATA_FLAGS) -ggdb


INCLUDE = $(if ${EIGEN}, -I${EIGEN})  $(if ${AMREX}, -I${AMREX}/include/) -I./src/ $(for pth in ${CPLUS_INCLUDE_PATH}; do echo -I"$pth"; done)
LIB     = -L${AMREX}/lib/ -lamrex 

HDR = $(shell find src/ -name *.H)
SRC = $(shell find src/ -mindepth 2  -name "*.cpp" )
SRC_F = $(shell find src/ -mindepth 2  -name "*.F90" )
SRC_MAIN = $(shell find src/ -maxdepth 1  -name "*.cc" )
EXE = $(subst src/,bin/, $(SRC_MAIN:.cc=)) 
OBJ = $(subst src/,obj/, $(SRC:.cpp=.cpp.o)) 
OBJ_MAIN = $(subst src/,obj/, $(SRC_MAIN:.cpp=.cc.o))
OBJ_F = $(subst src/,obj/, $(SRC_F:.F90=.F90.o))

.SECONDARY: 


default: $(EXE)
	@printf "$(B_ON)$(FG_GREEN)###\n"
	@printf "$(B_ON)$(FG_GREEN)### DONE\n" 
	@printf "$(B_ON)$(FG_GREEN)###$(RESET)\n"

bin/%: ${OBJ_F} ${OBJ} obj/%.cc.o
	@printf "$(B_ON)$(FG_BLUE)###\n"
	@printf "$(B_ON)$(FG_BLUE)### LINKING $@\n" 
	@printf "$(B_ON)$(FG_BLUE)###$(RESET)\n"
	mkdir -p bin/
	$(CC) -o $@ $^ ${LIB}  ${MPI_LIB}

obj/%.cc.o: src/%.cc ${HDR}
	@printf "$(B_ON)$(FG_YELLOW)###\n"
	@printf "$(B_ON)$(FG_YELLOW)### COMPILING $<\n" 
	@printf "$(B_ON)$(FG_YELLOW)###$(RESET)\n"
	@mkdir -p $(dir $@)
	$(CC) -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} ${MPICXX_COMPILE_FLAGS}

obj/%.cpp.o: src/%.cpp ${HDR}
	@printf "$(B_ON)$(FG_YELLOW)###\n"
	@printf "$(B_ON)$(FG_YELLOW)### COMPILING $<\n" 
	@printf "$(B_ON)$(FG_YELLOW)###$(RESET)\n"
	@mkdir -p $(dir $@)
	$(CC) -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} ${MPICXX_COMPILE_FLAGS}

obj/IO/WriteMetaData.cpp.o: .FORCE
	@printf "$(B_ON)$(FG_CYAN)###\n"
	@printf "$(B_ON)$(FG_CYAN)### COMPILING $@\n" 
	@printf "$(B_ON)$(FG_CYAN)###$(RESET)\n"
	@mkdir -p $(dir $@)
	$(CC) -c ${subst obj/,src/,${@:.cpp.o=.cpp}} -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} ${MPICXX_COMPILE_FLAGS}
.PHONY: .FORCE



FORT_INCL = $(shell for i in ${CPLUS_INCLUDE_PATH//:/ }; do echo -I'$i'; done)

obj/%.F90.o: src/%.F90 
	@printf "$(B_ON)$(FG_YELLOW)###\n"
	@printf "$(B_ON)$(FG_YELLOW)### COMPILING $<\n" 
	@printf "$(B_ON)$(FG_YELLOW)###$(RESET)\n"
	@mkdir -p $(dir $@)
	mpif90 -c $< -o $@  -I${subst :, -I,$(CPLUS_INCLUDE_PATH)}
	rm *.mod -rf

clean:
	@printf "$(B_ON)$(FG_RED)###\n"
	@printf "$(B_ON)$(FG_RED)### CLEANING\n" 
	@printf "$(B_ON)$(FG_RED)###$(RESET)\n"
	find src/ -name "*.o" -exec rm {} \;
	rm -f bin/*
	rm -rf obj/
	rm -f Backtrace*

%.cpp: $(HDR)

%.cc: $(HDR)

%.F90: $(HDR)

%.H :

help:
	@printf "$(B_ON)$(FG_YELLOW)\n\n============================== ALAMO Makefile help ==============================$(RESET)""\n\n"
	@printf "$(B_ON)Overview: \n$(RESET)"
	@printf "   This makefile automatically compiles all .cpp and .F90 files in \n"
	@printf "   the src directory, and compiles AND LINKS all .cc files into an \n"
	@printf "   executable in the bin directory. \n"
	@printf "   Any modification to a .H file causes everything to recompile. \n"
	@printf "   The file WriteMetaData.cpp recompiles every time to ensure that \n"
	@printf "   all metadata macros are up-to-date. \n"
	@printf "$(B_ON)Usage: $(RESET)\n"
	@printf "$(FG_LIGHTGREEN)   make [exe name] [COMP=INTEL/GCC] [EIGEN=/path/to/eigen]$(RESET)\n"
	@printf "$(FG_LIGHTGREEN)        [ALAMO=/path/to/alamo] [MPICHFORT=mpichfort,mpichf90] [-jNUM]$(RESET) \n"
	@printf "$(B_ON)Examples: $(RESET)\n"
	@printf "$(FG_LIGHTGREEN)   make                              $(RESET) (makes everything using default options)\n"
	@printf "$(FG_LIGHTGREEN)   make bin/alamo                    $(RESET) (makes bin/alamo only)\n"
	@printf "$(FG_LIGHTGREEN)   make COMP=INTEL                   $(RESET) (make using Intel compiler options)\n"
	@printf "$(FG_LIGHTGREEN)   make AMREX=/path/to/amrex         $(RESET) (specify location of AMReX)\n"
	@printf "$(FG_LIGHTGREEN)   make EIGEN=/path/to/eigen         $(RESET) (specify location of Eigen)\n" 
	@printf "$(FG_LIGHTGREEN)   make MPICHFORT=mpichfort,mpichf90 $(RESET) (backwards compatibility for mpichfort)\n" 
	@printf "$(FG_LIGHTGREEN)   make -j8                          $(RESET) (compile in parallel with 8 processors)\n"
	@printf "$(B_ON)Notes: $(RESET)\n"
	@printf "   - Specifying AMREX and EIGEN paths $(FG_LIGHTRED)does not$(RESET) override libraries\n"
	@printf "     that are already loaded in path.   \n"
	@printf "   - The AMREX path must contain directories called $(FG_LIGHTBLUE)lib/ include/$(RESET)   \n"
	@printf "   - The EIGEN path must contain a directory called $(FG_LIGHTBLUE)eigen3$(RESET)   \n"
	@printf "\n"
