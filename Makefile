CC = mpicxx

COMP ?= GCC
ifeq ($(COMP),INTEL)
 MPI_LIB = -lifcore
else ifeq ($(COMP),GCC)
 MPI_LIB = -lgfortran -lmpichfort -lmpich
endif

RESET              = '\033[0m'
B_ON               = '\033[1m'
FG_RED             = '\033[31m'
FG_GREEN           = '\033[32m'
FG_YELLOW          = '\033[33m'
FG_BLUE            = '\033[34m'
FG_CYAN            = '\033[36m'

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


INCLUDE = -I${EIGEN} -I${AMREX}/include/ -I./src/ $(for pth in ${CPLUS_INCLUDE_PATH}; do echo -I"$pth"; done)
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
	@printf $(B_ON)$(FG_GREEN)"###\n"
	@printf "### DONE\n" 
	@printf "###"$(RESET)"\n"

bin/%: ${OBJ_F} ${OBJ} obj/%.cc.o
	@printf $(B_ON)$(FG_BLUE)"###\n"
	@printf "### LINKING $@\n" 
	@printf "###"$(RESET)"\n"
	mkdir -p bin/
	$(CC) -o $@ $^ ${LIB}  ${MPI_LIB}

obj/%.cc.o: src/%.cc ${HDR}
	@printf $(B_ON)$(FG_YELLOW)"###\n"
	@printf "### COMPILING $<\n" 
	@printf "###"$(RESET)"\n"
	@mkdir -p $(dir $@)
	$(CC) -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} ${MPICXX_COMPILE_FLAGS}

obj/%.cpp.o: src/%.cpp ${HDR}
	@printf $(B_ON)$(FG_YELLOW)"###\n"
	@printf "### COMPILING $<\n" 
	@printf "###"$(RESET)"\n"
	@mkdir -p $(dir $@)
	$(CC) -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} ${MPICXX_COMPILE_FLAGS}

obj/IO/WriteMetaData.cpp.o: .FORCE
	@echo $(TEST_VAR)
	@echo $(AMREX)
	@printf $(B_ON)$(FG_CYAN)"###\n"
	@printf "### COMPILING $@\n" 
	@printf "###"$(RESET)"\n"
	@mkdir -p $(dir $@)
	$(CC) -c ${subst obj/,src/,${@:.cpp.o=.cpp}} -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} ${MPICXX_COMPILE_FLAGS}
.PHONY: .FORCE



FORT_INCL = $(shell for i in ${CPLUS_INCLUDE_PATH//:/ }; do echo -I'$i'; done)

obj/%.F90.o: src/%.F90 
	@printf $(B_ON)$(FG_YELLOW)"###\n"
	@printf "### COMPILING $<\n" 
	@printf "###"$(RESET)"\n"
	@mkdir -p $(dir $@)
	mpif90 -c $< -o $@  -I${subst :, -I,$(CPLUS_INCLUDE_PATH)}
	rm *.mod -rf

clean:
	@printf $(B_ON)$(FG_RED)"###\n"
	@printf "### CLEANING\n" 
	@printf "###"$(RESET)"\n"
	find src/ -name "*.o" -exec rm {} \;
	rm -f bin/*
	rm -rf obj/
	rm -f Backtrace*

%.cpp: $(HDR)

%.cc: $(HDR)

%.F90: $(HDR)

%.H :
