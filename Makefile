
-include .make/Makefile.pre.conf

AMREX_TARGET ?= 
CC ?= mpicxx -cxx=g++
MPI_LIB ?= -lmpich

RESET              = \033[0m
B_ON               = \033[1m
FG_RED             = \033[31m
FG_DIM             = \033[2m
FG_LIGHTRED        = \033[91m
FG_LIGHTGRAY       = \033[37m
FG_GRAY            = \033[90m
FG_GREEN           = \033[32m
FG_LIGHTGREEN      = \033[92m
FG_YELLOW          = \033[33m
FG_LIGHTYELLOW     = \033[93m
FG_BLUE            = \033[34m
FG_LIGHTBLUE       = \033[94m
FG_CYAN            = \033[36m
FG_MAGENTA         = \033[35m




QUIET ?= @


METADATA_GITHASH  = $(shell git describe --always --dirty)
METADATA_USER     = $(shell whoami)
METADATA_PLATFORM = $(shell hostname)
METADATA_COMPILER = $(COMP)
METADATA_DATE     = $(shell date +%x)
METADATA_TIME     = $(shell date +%H:%M:%S)
BUILD_DIR         = ${shell pwd}

METADATA_FLAGS = -DMETADATA_GITHASH=\"$(METADATA_GITHASH)\" -DMETADATA_USER=\"$(METADATA_USER)\" -DMETADATA_PLATFORM=\"$(METADATA_PLATFORM)\" -DMETADATA_COMPILER=\"$(METADATA_COMPILER)\" -DMETADATA_DATE=\"$(METADATA_DATE)\" -DMETADATA_TIME=\"$(METADATA_TIME)\" -DBUILD_DIR=\"${BUILD_DIR}\" $(if ${MEME}, -DMEME)


CXX_COMPILE_FLAGS += -Winline -Wextra -Wall -Wno-comment -std=c++17 $(METADATA_FLAGS)

LINKER_FLAGS += -Bsymbolic-functions -lstdc++fs

#CXX_COMPILE_FLAGS += --param inline-unit-growth=100 --param  max-inline-insns-single=1200
#LINKER_FLAGS      += --param inline-unit-growth=100 --param  max-inline-insns-single=1200


ALAMO_INCLUDE += $(if ${EIGEN}, -isystem ${EIGEN})  $(if ${AMREX}, -isystem ${AMREX}/include/) -I./src/ $(for pth in ${CPLUS_INCLUDE_PATH}; do echo -I"$pth"; done)
LIB     += -L${AMREX}/lib/ -lamrex -lpthread

HDR_ALL = $(shell find src/ -name *.H)
HDR_TEST = $(shell find src/ -name *Test.H)
HDR = $(filter-out $(HDR_TEST),$(HDR_ALL))
SRC = $(shell find src/ -mindepth 2  -name "*.cpp" )
SRC_F = $(shell find src/ -mindepth 2  -name "*.F90" )
SRC_MAIN = $(shell find src/ -maxdepth 1  -name "*.cc" )
EXE = $(subst src/,bin/, $(SRC_MAIN:.cc=-$(POSTFIX))) 
OBJ = $(subst src/,obj/obj-$(POSTFIX)/, $(SRC:.cpp=.cpp.o)) 
DEP = $(subst src/,obj/obj-$(POSTFIX)/, $(SRC:.cpp=.cpp.d)) $(subst src/,obj/obj-$(POSTFIX)/, $(SRC_MAIN:.cc=.cc.d))
OBJ_MAIN = $(subst src/,obj/obj-$(POSTFIX)/, $(SRC_MAIN:.cpp=.cc.o))
OBJ_F = $(subst src/,obj/obj-$(POSTFIX)/, $(SRC_F:.F90=.F90.o))

NUM = $(words $(SRC) $(SRC_F) $(SRC_MAIN))
CTR = 0
NUM_DEP = $(words $(DEP))
CTR_DEP = 0
NUM_EXE = $(words $(EXE))
CTR_EXE = 0

.SECONDARY: 





default: $(DEP) $(EXE)
	@printf "$(B_ON)$(FG_GREEN)DONE $(RESET)\n" 

tidy:
	@printf "$(B_ON)$(FG_RED)TIDYING  $(RESET)\n" 
	find src -name "*.orig" -exec rm -rf {} \;
	rm -f Backtrace*
	rm -f amrex.build.log

clean: tidy
	@printf "$(B_ON)$(FG_RED)CLEANING  $(RESET)\n" 
	find src/ -name "*.o" -exec rm {} \;
	rm -rf .diff*
	rm -f bin/*
	rm -rf obj
	rm -f Backtrace*
	rm -rf docs/build docs/doxygen docs/html docs/latex
	rm -f amrex.build.log

realclean: clean
	@printf "$(B_ON)$(FG_RED)CLEANING AMREX $(RESET)\n" 
	-make -C ext/amrex realclean
	git -C ext/amrex reset --hard
	git -C ext/amrex clean -fd
	git -C ext/amrex clean -fx
	rm -rf ext/amrex/1d* ext/amrex/2d* ext/amrex/3d*
	@printf "$(B_ON)$(FG_RED)CLEANING OLD CONFIGURATIONS $(RESET)\n" 
	rm -rf Makefile.conf Makefile.amrex.conf .make


info:
	@printf "$(B_ON)$(FG_BLUE)Compiler version information$(RESET)\n"
	$(CC) --version

-include .make/Makefile.post.conf

bin/%: bin/%-$(POSTFIX) ;

bin/%-$(POSTFIX): ${OBJ_F} ${OBJ} obj/obj-$(POSTFIX)/%.cc.o 
	$(eval CTR_EXE=$(shell echo $$(($(CTR_EXE)+1))))
	@printf "$(B_ON)$(FG_BLUE)LINKING$(RESET)$(FG_LIGHTBLUE)     " 
	@printf '%9s' "($(CTR_EXE)/$(NUM_EXE)) " 
	@printf "$(RESET)$@\n"
	@mkdir -p bin/
	$(QUIET)$(CC) -o $@ $^ ${LIB}  ${MPI_LIB}  ${LINKER_FLAGS}

obj/obj-$(POSTFIX)/test.cc.o: src/test.cc ${AMREX_TARGET}
	$(eval CTR=$(shell echo $$(($(CTR)+1))))
	@printf "$(B_ON)$(FG_YELLOW)COMPILING$(RESET)$(FG_LIGHTYELLOW)   "
	@printf '%9s' "($(CTR)/$(NUM)) " 
	@printf "$(RESET)$<\n"
	@mkdir -p $(dir $@)
	$(QUIET)$(CC) -c $< -o $@ ${ALAMO_INCLUDE} ${CXX_COMPILE_FLAGS} 

obj/obj-$(POSTFIX)/%.cc.o: src/%.cc ${AMREX_TARGET} 
	$(eval CTR=$(shell echo $$(($(CTR)+1))))
	@printf "$(B_ON)$(FG_YELLOW)COMPILING$(RESET)$(FG_LIGHTYELLOW)   "
	@printf '%9s' "($(CTR)/$(NUM)) " 
	@printf "$(RESET)$<\n"
	@mkdir -p $(dir $@)
	$(QUIET)$(CC) -c $< -o $@ ${ALAMO_INCLUDE} ${CXX_COMPILE_FLAGS} 

obj/obj-$(POSTFIX)/%.cpp.o: 
	$(eval CTR=$(shell echo $$(($(CTR)+1))))
	@printf "$(B_ON)$(FG_YELLOW)COMPILING$(RESET)$(FG_LIGHTYELLOW)   "
	@printf '%9s' "($(CTR)/$(NUM)) " 
	@printf "$(RESET)$<\n"
	@mkdir -p $(dir $@)
	$(QUIET)$(CC) -c $< -o $@ ${ALAMO_INCLUDE} ${CXX_COMPILE_FLAGS} 

obj/obj-$(POSTFIX)/%.cpp.d: src/%.cpp  ${AMREX_TARGET}
	$(eval CTR_DEP=$(shell echo $$(($(CTR_DEP)+1))))
	@printf "$(B_ON)$(FG_GRAY)DEPENDENCY$(RESET)$(FG_LIGHTGRAY)  " 
	@printf '%9s' "($(CTR_DEP)/$(NUM)) " 
	@printf "$(RESET)$<\n"
	@mkdir -p $(dir $@)
	$(QUIET)$(CC) -Wno-unused-command-line-argument -I./src/ $< ${ALAMO_INCLUDE} ${CXX_COMPILE_FLAGS}-MM -MT $(@:.cpp.d=.cpp.o) -MF $@

obj/obj-$(POSTFIX)/%.cc.d: src/%.cc ${AMREX_TARGET}
	$(eval CTR_DEP=$(shell echo $$(($(CTR_DEP)+1))))
	@printf "$(B_ON)$(FG_GRAY)DEPENDENCY$(RESET)$(FG_LIGHTGRAY)  " 
	@printf '%9s' "($(CTR_DEP)/$(NUM)) " 
	@printf "$(RESET)$<\n"
	@mkdir -p $(dir $@)
	$(QUIET)$(CC) -Wno-unused-command-line-argument -I./src/ $< ${ALAMO_INCLUDE} ${CXX_COMPILE_FLAGS} -MM -MT $(@:.cc.d=.cc.o) -MF $@

obj/obj-$(POSTFIX)/IO/WriteMetaData.cpp.o: .FORCE ${AMREX_TARGET} ${DEP_DIFF}
	$(eval CTR=$(shell echo $$(($(CTR)+1))))
	@printf "$(B_ON)$(FG_LIGHTYELLOW)COMPILING$(RESET)$(FG_LIGHTYELLOW)   "
	@printf '%9s' "($(CTR)/$(NUM)) " 
	@printf "$(RESET)${subst obj/obj-$(POSTFIX)/,src/,${@:.cpp.o=.cpp}} \n"
	@mkdir -p $(dir $@)
	$(QUIET)$(CC) -c ${subst obj/obj-$(POSTFIX)/,src/,${@:.cpp.o=.cpp}} -o $@ ${ALAMO_INCLUDE} ${CXX_COMPILE_FLAGS} 

.PHONY: .FORCE

FORT_INCL = $(shell for i in ${CPLUS_INCLUDE_PATH//:/ }; do echo -I'$i'; done)

obj/obj-$(POSTFIX)/%.F90.o: src/%.F90 
	@printf "$(B_ON)$(FG_YELLOW)COMPILING  $(RESET)$<\n" 
	@mkdir -p $(dir $@)
	mpif90 -c $< -o $@  -I${subst :, -I,$(CPLUS_INCLUDE_PATH)}
	rm *.mod -rf

docs: docs/build/html/index.html .FORCE 
	@printf "$(B_ON)$(FG_MAGENTA)DOCS$(RESET) Done\n" 

docs/build/html/index.html: $(shell find docs/source/ -type f) README.rst .FORCE
	@printf "$(B_ON)$(FG_MAGENTA)DOCS$(RESET) Generating sphinx\n" 	
	@make -C docs html # > /dev/null


check: .FORCE
	@./scripts/checkdoc.py
	@./.github/workflows/style/check_tabs.py
	@eclint check src

test: .FORCE
	@./.github/workflows/style/check_tabs.py
	@make docs
	@./scripts/runtests.py

GCDA = $(shell mkdir -p obj && find obj/ -name "*.gcda")
GCNO = $(shell mkdir -p obj && find obj/ -name "*.gcno")

GCDA_DIRS  = $(shell mkdir -p obj && find obj/ -maxdepth 1 -name "*coverage*" )
GCDA_DIMS  = $(subst obj-,,$(subst -coverage-g++,,$(notdir $(GCDA_DIRS))))
GCDA_INFOS = $(subst obj-,cov/coverage_,$(subst -coverage-g++,.info,$(notdir $(GCDA_DIRS))))
GCDA_LCOVS = $(subst obj-,--add-tracefile cov/coverage_,$(subst -coverage-g++,.info,$(notdir $(GCDA_DIRS))))

cov-report: cov/index.html
	@echo $(GCDA_LCOVS)
	@echo "Done - output in cov/index.html"

cov-clean: .FORCE
	rm -rf $(GCDA)
	rm -rf ./cov

cov/index.html: cov/coverage_merged.info
	genhtml cov/coverage_merged.info --output-directory cov

cov/coverage_merged.info: $(GCDA_INFOS)
	mkdir -p ./cov/
	lcov --ignore-errors=gcov,source,graph $(GCDA_LCOVS) -o cov/coverage_merged.info  

cov/coverage_%.info: obj/obj-%-coverage-g++/ $(GCDA)
	mkdir -p ./cov/
	geninfo $< -b . -o $@ --exclude "/usr/*" --exclude "ext/*"

githubpages: docs cov-report
	mkdir -p ./githubpages/
	echo "<head><meta http-equiv=\"refresh\" content=\"0; url='docs/index.html\" /></head>" > githubpages/index.html
	cp -rf docs/build/html ./githubpages/docs/
	cp -rf cov/ ./githubpages/cov/

ifneq ($(MAKECMDGOALS),tidy)
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),realclean)
ifneq ($(MAKECMDGOALS),info)
ifneq ($(MAKECMDGOALS),help)
ifneq ($(MAKECMDGOALS),docs)
-include $(DEP)
endif
endif
endif
endif
endif
endif

