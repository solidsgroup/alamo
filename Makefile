
-include Makefile.conf

COMP ?= GCC
ifeq ($(COMP),INTEL)
 CC = mpicxx -cxx=icc
 MPI_LIB = -lifcore
else ifeq ($(COMP),GCC)
 CC = mpicxx -cxx=g++
 MPI_LIB = -lgfortran -lmpich
else ifeq ($(COMP),CLANG)
 CC = mpicxx -cxx=clang++
 MPI_LIB = -lgfortran -lmpich
endif

RESET              = \033[0m
B_ON               = \033[1m
FG_RED             = \033[31m
FG_DIM             = \033[2m
FG_LIGHTRED        = \033[91m
FG_LIGHTGRAY       = \033[90m
FG_GREEN           = \033[32m
FG_LIGHTGREEN      = \033[92m
FG_YELLOW          = \033[33m
FG_LIGHTYELLOW     = \033[93m
FG_BLUE            = \033[34m
FG_LIGHTBLUE       = \033[94m
FG_CYAN            = \033[36m
FG_MAGENTA         = \033[35m




METADATA_GITHASH  = $(shell git log --pretty=format:'%H' -n 1)
METADATA_USER     = $(shell whoami)
METADATA_PLATFORM = $(shell hostname)
METADATA_COMPILER = $(COMP)
METADATA_DATE     = $(shell date +%x)
METADATA_TIME     = $(shell date +%H:%M:%S)
BUILD_DIR         = ${shell pwd}

METADATA_FLAGS = -DMETADATA_GITHASH=\"$(METADATA_GITHASH)\" -DMETADATA_USER=\"$(METADATA_USER)\" -DMETADATA_PLATFORM=\"$(METADATA_PLATFORM)\" -DMETADATA_COMPILER=\"$(METADATA_COMPILER)\" -DMETADATA_DATE=\"$(METADATA_DATE)\" -DMETADATA_TIME=\"$(METADATA_TIME)\" -DBUILD_DIR=\"${BUILD_DIR}\" $(if ${MEME}, -DMEME)


CXX_COMPILE_FLAGS = -fPIC -Winline -Wpedantic -Wextra -Wall  -std=c++11 $(METADATA_FLAGS)
ifeq ($(DEBUG),TRUE)
 CXX_COMPILE_FLAGS += -ggdb -g3
else 
 CXX_COMPILE_FLAGS += -O3
endif

LINKER_FLAGS = -Bsymbolic-functions

INCLUDE = $(if ${EIGEN}, -isystem ${EIGEN})  $(if ${AMREX}, -isystem ${AMREX}/include/) -I./src/ $(for pth in ${CPLUS_INCLUDE_PATH}; do echo -I"$pth"; done) -I/usr/include/python2.7/
LIB     = -L${AMREX}/lib/ -lamrex 

HDR_ALL = $(shell find src/ -name *.H)
HDR_TEST = $(shell find src/ -name *Test.H)
HDR = $(filter-out $(HDR_TEST),$(HDR_ALL))
SRC = $(shell find src/ -mindepth 2  -name "*.cpp" )
SRC_F = $(shell find src/ -mindepth 2  -name "*.F90" )
SRC_MAIN = $(shell find src/ -maxdepth 1  -name "*.cc" )
EXE = $(subst src/,bin/, $(SRC_MAIN:.cc=$(PREFIX))) 
OBJ = $(subst src/,obj/obj$(PREFIX)/, $(SRC:.cpp=.cpp.o)) 
DEP = $(subst src/,obj/obj$(PREFIX)/, $(SRC:.cpp=.cpp.d)) $(subst src/,obj/obj$(PREFIX)/, $(SRC_MAIN:.cc=.cc.d))
OBJ_MAIN = $(subst src/,obj/obj$(PREFIX)/, $(SRC_MAIN:.cpp=.cc.o))
OBJ_F = $(subst src/,obj/obj$(PREFIX)/, $(SRC_F:.F90=.F90.o))

.SECONDARY: 

default: $(DEP) $(EXE)
	@printf "$(B_ON)$(FG_GREEN)DONE $(RESET)\n" 

python: ${OBJ} 
	$(CC) -c src/python.cc -fPIC -o alamo.o ${INCLUDE} ${CXX_COMPILE_FLAGS} 
	g++ -shared -Wl,-soname,alamo.so -o alamo.so  alamo.o ${OBJ} ${LIB} ${MPI_LIB}  /usr/lib/x86_64-linux-gnu/libboost_python-py27.so.1.62.0


clean:
	@printf "$(B_ON)$(FG_RED)CLEANING  $(RESET)\n" 
	find src/ -name "*.o" -exec rm {} \;
	rm -f bin/*
	rm -rf obj
	rm -f Backtrace*
	rm -f Makefile.conf
	rm -rf docs/build docs/doxygen docs/html docs/latex
tidy:
	@printf "$(B_ON)$(FG_RED)TIDYING  $(RESET)\n" 
	rm -f Backtrace*
info:
	@printf "$(B_ON)$(FG_BLUE)Compiler version information$(RESET)\n"
	@$(CC) --version
	@printf "$(B_ON)$(FG_BLUE)MPI Flags$(RESET)\n"
	@$(CC) -show
	@printf "\n"


bin/%: bin/%$(PREFIX) ;

bin/%$(PREFIX): ${OBJ_F} ${OBJ} obj/obj$(PREFIX)/%.cc.o
	@printf "$(B_ON)$(FG_BLUE)LINKING$(RESET)     $@ \n" 
	@mkdir -p bin/
	@$(CC) -o $@ $^ ${LIB}  ${MPI_LIB}  ${LINKER_FLAGS}

obj/obj$(PREFIX)/test.cc.o: src/test.cc
	@printf "$(B_ON)$(FG_YELLOW)COMPILING$(RESET)   $< \n" 
	@mkdir -p $(dir $@)
	@$(CC) -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} 

obj/obj$(PREFIX)/%.cc.o: src/%.cc
	@printf "$(B_ON)$(FG_YELLOW)COMPILING$(RESET)   $< \n" 
	@mkdir -p $(dir $@)
	@$(CC) -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} 

obj/obj$(PREFIX)/%.cpp.o: 
	@printf "$(B_ON)$(FG_YELLOW)COMPILING$(RESET)   $< \n" 
	@mkdir -p $(dir $@)
	@$(CC) -c $< -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} 

obj/obj$(PREFIX)/%.cpp.d: src/%.cpp 
	@printf "$(B_ON)$(FG_LIGHTGRAY)DEPENDENCY$(RESET)  $< \n" 
	@mkdir -p $(dir $@)
	@g++ -I./src/ $< ${INCLUDE} ${CXX_COMPILE_FLAGS} -MM -MT $(@:.cpp.d=.cpp.o) -MF $@

obj/obj$(PREFIX)/%.cc.d: src/%.cc
	@printf "$(B_ON)$(FG_LIGHTGRAY)DEPENDENCY$(RESET)  $< \n" 
	@mkdir -p $(dir $@)
	@g++ -I./src/ $< ${INCLUDE} ${CXX_COMPILE_FLAGS} -MM -MT $(@:.cc.d=.cc.o) -MF $@



obj/obj$(PREFIX)/IO/WriteMetaData.cpp.o: .FORCE
	@printf "$(B_ON)$(FG_LIGHTYELLOW)$(FG_DIM)COMPILING$(RESET)   ${subst obj/obj$(PREFIX)/,src/,${@:.cpp.o=.cpp}} \n" 
	@mkdir -p $(dir $@)
	@$(CC) -c ${subst obj/obj$(PREFIX)/,src/,${@:.cpp.o=.cpp}} -o $@ ${INCLUDE} ${CXX_COMPILE_FLAGS} 

.PHONY: .FORCE



FORT_INCL = $(shell for i in ${CPLUS_INCLUDE_PATH//:/ }; do echo -I'$i'; done)

obj/obj$(PREFIX)/%.F90.o: src/%.F90 
	@printf "$(B_ON)$(FG_YELLOW)COMPILING  $(RESET)$<\n" 
	@mkdir -p $(dir $@)
	mpif90 -c $< -o $@  -I${subst :, -I,$(CPLUS_INCLUDE_PATH)}
	rm *.mod -rf


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
	@printf "$(FG_LIGHTGREEN)        [ALAMO=/path/to/alamo] [-jNUM]$(RESET) \n"
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


docs: docs/doxygen/index.html docs/build/html/index.html .FORCE 
	@printf "$(B_ON)$(FG_MAGENTA)DOCS$(RESET) Done\n" 

docs/doxygen/index.html: $(SRC) $(SRC_F) $(SRC_MAIN) $(HDR_ALL)
	@printf "$(B_ON)$(FG_MAGENTA)DOCS$(RESET) Generating doxygen files\n" 	
	@cd docs && doxygen 
docs/build/html/index.html: $(shell find docs/source/ -type f) Readme.md
	@printf "$(B_ON)$(FG_MAGENTA)DOCS$(RESET) Generating sphinx\n" 	
	@make -C docs html > /dev/null


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),info)
ifneq ($(MAKECMDGOALS),help)
ifneq ($(MAKECMDGOALS),docs)
-include $(DEP)
endif
endif
endif
endif

