
-include .make/Makefile.pre.conf

AMREX_TARGET ?= 
CC ?= mpicxx -cxx=g++

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

FG_ORANGE          = \033[38;5;208m



QUIET ?= @


METADATA_GITHASH  = $(shell git describe --always --dirty)
METADATA_USER     = $(shell whoami)
METADATA_PLATFORM = $(shell hostname)
METADATA_COMPILER = $(COMP)
METADATA_DATE     = $(shell date +%x)
METADATA_TIME     = $(shell date +%H:%M:%S)
BUILD_DIR         = ${shell pwd}

METADATA_FLAGS = -DMETADATA_GITHASH=\"$(METADATA_GITHASH)\" -DMETADATA_USER=\"$(METADATA_USER)\" -DMETADATA_PLATFORM=\"$(METADATA_PLATFORM)\" -DMETADATA_COMPILER=\"$(METADATA_COMPILER)\" -DMETADATA_DATE=\"$(METADATA_DATE)\" -DMETADATA_TIME=\"$(METADATA_TIME)\" -DBUILD_DIR=\"${BUILD_DIR}\" $(if ${MEME}, -DMEME)


CXX_COMPILE_FLAGS += -Winline -Wextra -Wall -Wno-comment -std=c++20 $(METADATA_FLAGS)

# Clang in C++20 adds a warning for failing to explicitly capture "this" in lambda functions.
# Compliance will require changing ParallelFor calls from "[=]" to "[=,this]"
# g++ also issues a warning but has no instrumentation for suppressing this warning only, so we
# will just leave noisy output
ifeq ($(COMP), CLANG)
CXX_COMPILE_FLAGS += -Wno-deprecated-this-capture 
endif

# Skip on macOS: -Bsymbolic-functions is a GNU ld flag not supported by
# Apple's ld64, and -lstdc++fs is a pre-GCC-9 std::filesystem shim that
# doesn't exist in libc++ (used by Apple Clang and Homebrew LLVM).
ifneq ($(shell uname -s),Darwin)
LINKER_FLAGS += -Bsymbolic-functions -lstdc++fs
endif

#CXX_COMPILE_FLAGS += --param inline-unit-growth=100 --param  max-inline-insns-single=1200
#LINKER_FLAGS      += --param inline-unit-growth=100 --param  max-inline-insns-single=1200


ALAMO_INCLUDE += $(if ${EIGEN}, -isystem ${EIGEN})  $(if ${AMREX_TARGET}, -isystem ${AMREX_TARGET}/include/) -I./src/ $(foreach pth,$(subst :, ,${CPLUS_INCLUDE_PATH}),-I$(pth))
LIB     += ${AMREX_TARGET}/lib/libamrex.a -lpthread

HDR_ALL = $(shell find src/ -name *.H)
HDR_TEST = $(shell find src/ -name *Test.H)
HDR = $(filter-out $(HDR_TEST),$(HDR_ALL))
SRC = $(shell find src/ -mindepth 2  -name "*.cpp" )
SRC_MAIN = $(shell find src/ -maxdepth 1  -name "*.cc" )
EXE = $(subst src/,bin/, $(SRC_MAIN:.cc=-$(POSTFIX))) 
# Input-builder generation first ensures executables from the active configuration
# are current, then also includes executable files already present in bin/.
# Keep the complete basename so configurations do not overwrite each other's data.
INPUT_BUILDER_CONFIGURED_EXE = $(EXE)
INPUT_BUILDER_PRESENT_EXE = $(shell find bin -maxdepth 1 -type f -executable -print 2>/dev/null)
INPUT_BUILDER_EXE ?= $(sort $(INPUT_BUILDER_CONFIGURED_EXE) $(INPUT_BUILDER_PRESENT_EXE))
INPUT_BUILDER_EXECUTABLE_NAMES = $(notdir $(INPUT_BUILDER_EXE))
INPUT_BUILDER_ANNOTATE_SCRIPT ?= scripts/annotate_input_schema.py
INPUT_BUILDER_HTML_SCRIPT ?= scripts/make_input_builder.py
INPUT_BUILDER_HTML_TEMPLATE ?= scripts/input_builder_template.html
INPUT_BUILDER_CHECK_SCRIPT ?= scripts/check_input_scrape.py
INPUT_REFERENCE_SCRIPT ?= scripts/make_schema_reference.py
INPUT_BUILDER_CHECK_STRICT ?= 0
INPUT_SCHEMA_GENERATORS = $(INPUT_BUILDER_ANNOTATE_SCRIPT)
INPUT_HTML_GENERATORS = $(INPUT_BUILDER_HTML_SCRIPT) $(INPUT_BUILDER_HTML_TEMPLATE)
INPUT_COVERAGE_GENERATORS = $(INPUT_BUILDER_CHECK_SCRIPT)
DOC_INPUT_SCHEMA_DIR ?= docs/source/_static/input-schemas
DOC_INPUT_BUILDER_DIR ?= docs/source/_static/input-builders
DOC_INPUT_BUILDER_INDEX ?= $(DOC_INPUT_BUILDER_DIR)/index.html
DOC_INPUT_COVERAGE_REPORT ?= $(DOC_INPUT_SCHEMA_DIR)/input-scrape-coverage.json
DOC_INPUT_REFERENCE ?= docs/source/Inputs.generated.rst
DOC_INPUT_REFERENCE_DIR ?= docs/source/InputsReference
DOC_DOXYGEN_INDEX ?= docs/build/html/doxygen/index.html
DOC_DOXYGEN_INPUTS = $(shell find src -type f) README.rst docs/Doxyfile docs/Makefile
DOC_INPUT_SCHEMAS = $(addprefix $(DOC_INPUT_SCHEMA_DIR)/,$(addsuffix .schema.json,$(INPUT_BUILDER_EXECUTABLE_NAMES)))
DOC_INPUT_BUILDERS = $(addprefix $(DOC_INPUT_BUILDER_DIR)/,$(addsuffix .html,$(INPUT_BUILDER_EXECUTABLE_NAMES)))
DOC_SOURCE_FILES = $(shell find docs/source/ -type f ! -path '$(DOC_INPUT_SCHEMA_DIR)/*' ! -path '$(DOC_INPUT_BUILDER_DIR)/*' ! -path '$(DOC_INPUT_REFERENCE_DIR)/*' ! -name '$(notdir $(DOC_INPUT_REFERENCE))')
OBJ = $(subst src/,obj/obj-$(POSTFIX)/, $(SRC:.cpp=.cpp.o)) 
DEP = $(subst src/,obj/obj-$(POSTFIX)/, $(SRC:.cpp=.cpp.d)) $(subst src/,obj/obj-$(POSTFIX)/, $(SRC_MAIN:.cc=.cc.d))
OBJ_MAIN = $(subst src/,obj/obj-$(POSTFIX)/, $(SRC_MAIN:.cpp=.cc.o))

NUM = $(words $(SRC) $(SRC_MAIN))
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
	rm -f profile.prof*

clean: tidy
	@printf "$(B_ON)$(FG_RED)CLEANING  $(RESET)\n" 
	find src/ -name "*.o" -exec rm {} \;
	rm -rf .diff*
	rm -f bin/*
	rm -rf obj
	rm -f Backtrace*
	rm -rf docs/build docs/doxygen docs/html docs/latex
	rm -f amrex.build.log

clean-tests:
	@printf "$(B_ON)$(FG_RED)CLEANING TEST OUTPUT DIRECTORIES $(RESET)\n"
	rm -rf tests/*/output*
	rm -rf report/*
	rm -f report.html

realclean: clean
	@printf "$(B_ON)$(FG_RED)CLEANING AMREX $(RESET)\n" 
	-make -C ${AMREX_ROOT} realclean
	git -C ${AMREX_ROOT} reset --hard
	git -C ${AMREX_ROOT} clean -fd
	git -C ${AMREX_ROOT} clean -fx
	rm -rf ${AMREX_ROOT}/1d* ${AMREX_ROOT}/2d* ${AMREX_ROOT}/3d*
	@printf "$(B_ON)$(FG_RED)CLEANING OLD CONFIGURATIONS $(RESET)\n" 
	rm -rf Makefile.conf Makefile.amrex.conf .make

py: python_ok lib/libalamo-$(POSTFIX).so ${AMREX_TARGET}/lib/libamrex.so
	@python3 ./scripts/make_alamo_package.py --postfix=$(POSTFIX) --amrex=$(AMREX_TARGET)
	@printf "$(B_ON)$(FG_GREEN)DONE $(RESET)\n" 

info:
	@printf "$(B_ON)$(FG_BLUE)Compiler version information$(RESET)\n"
	$(CC) --version

-include .make/Makefile.post.conf

bin/%: bin/%-$(POSTFIX) ;

bin/%-$(POSTFIX): ${OBJ} obj/obj-$(POSTFIX)/%.cc.o
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

obj/obj-$(POSTFIX)/IO/WriteMetaData.cpp.o: src/IO/WriteMetaData.cpp ${AMREX_TARGET} | ${DEP_DIFF}
	$(eval CTR=$(shell echo $$(($(CTR)+1))))
	@printf "$(B_ON)$(FG_LIGHTYELLOW)COMPILING$(RESET)$(FG_LIGHTYELLOW)   "
	@printf '%9s' "($(CTR)/$(NUM)) " 
	@printf "$(RESET)${subst obj/obj-$(POSTFIX)/,src/,${@:.cpp.o=.cpp}} \n"
	@mkdir -p $(dir $@)
	$(QUIET)$(CC) -c ${subst obj/obj-$(POSTFIX)/,src/,${@:.cpp.o=.cpp}} -o $@ ${ALAMO_INCLUDE} ${CXX_COMPILE_FLAGS} 

.PHONY: .FORCE input-builders docs-input-builders input-scrape-check input-reference docs-with-input-builders

docs: docs/build/html/index.html .FORCE
	@printf "$(B_ON)$(FG_MAGENTA)DOCS$(RESET) Done\n" 

docs/build/html/index.html: input-reference $(DOC_SOURCE_FILES) README.rst .FORCE
	@printf "$(B_ON)$(FG_MAGENTA)DOCS$(RESET) Generating sphinx\n" 	
	@make -C docs html SKIP_DOXYGEN=1 # > /dev/null

input-builders: docs-input-builders

docs-input-builders: $(INPUT_BUILDER_CONFIGURED_EXE) $(DOC_INPUT_BUILDER_INDEX)
	@printf "$(B_ON)$(FG_MAGENTA)INPUT BUILDERS$(RESET) Done\n"

input-scrape-check: $(DOC_INPUT_COVERAGE_REPORT)

input-reference: docs-input-builders $(DOC_DOXYGEN_INDEX) $(INPUT_REFERENCE_SCRIPT)
	@printf "$(B_ON)$(FG_MAGENTA)INPUT REFERENCE$(RESET) Sphinx\n"
	@python3 "$(INPUT_REFERENCE_SCRIPT)" \
		--repo-root "." \
		--schema-dir "$(DOC_INPUT_SCHEMA_DIR)" \
		$(foreach schema,$(DOC_INPUT_SCHEMAS),--schema "$(schema)") \
		--output "$(DOC_INPUT_REFERENCE)" \
		--output-dir "$(DOC_INPUT_REFERENCE_DIR)"

$(DOC_DOXYGEN_INDEX): $(DOC_DOXYGEN_INPUTS)
	@$(MAKE) -C docs doxygen

$(DOC_INPUT_SCHEMA_DIR)/%.schema.json: bin/% $(INPUT_SCHEMA_GENERATORS) Makefile
	@printf "$(B_ON)$(FG_MAGENTA)INPUT SCHEMA$(RESET) $*\n"
	@mkdir -p "$(DOC_INPUT_SCHEMA_DIR)"
	@set -e; \
	tmp="$@.tmp"; \
	log=$$(mktemp); \
	trap 'rm -f "$$tmp" "$$log"' EXIT; \
	if ! "$<" --parse-args --parse-args-output "$$tmp" > "$$log" 2>&1; then \
		cat "$$log"; \
		exit 1; \
	fi; \
	python3 "$(INPUT_BUILDER_ANNOTATE_SCRIPT)" --schema "$$tmp" > /dev/null; \
	mv "$$tmp" "$@"; \
	rm -f "$$log"; \
	trap - EXIT

$(DOC_INPUT_BUILDER_DIR)/%.html: $(DOC_INPUT_SCHEMA_DIR)/%.schema.json $(INPUT_HTML_GENERATORS) $(DOC_DOXYGEN_INDEX) Makefile
	@printf "$(B_ON)$(FG_MAGENTA)INPUT BUILDER$(RESET) $*\n"
	@mkdir -p "$(DOC_INPUT_BUILDER_DIR)"
	@python3 "$(INPUT_BUILDER_HTML_SCRIPT)" --schema "$<" --output "$@.tmp" --title "Alamo Input Builder: $*" --current-builder "$*" --doxygen-dir "$(dir $(DOC_DOXYGEN_INDEX))" $(foreach builder,$(INPUT_BUILDER_EXECUTABLE_NAMES),--builder "$(builder)") > /dev/null
	@mv "$@.tmp" "$@"

$(DOC_INPUT_COVERAGE_REPORT): $(DOC_INPUT_SCHEMAS) $(INPUT_COVERAGE_GENERATORS) Makefile
	@printf "$(B_ON)$(FG_MAGENTA)INPUT CHECK$(RESET) Scrape coverage\n"
	@mkdir -p "$(dir $@)"
	@set -e; \
	tmp="$@.tmp"; \
	trap 'rm -f "$$tmp"' EXIT; \
	strict_flags=""; \
	case "$(INPUT_BUILDER_CHECK_STRICT)" in \
		1|yes|YES|true|TRUE) strict_flags="--fail-on-missing";; \
	esac; \
	python3 "$(INPUT_BUILDER_CHECK_SCRIPT)" \
		--repo-root "." \
		--source-dir "src" \
		--schema-dir "$(DOC_INPUT_SCHEMA_DIR)" \
		$(foreach schema,$(DOC_INPUT_SCHEMAS),--schema "$(schema)") \
		--report "$$tmp" \
		--max-list 25 \
		$$strict_flags; \
	mv "$$tmp" "$@"; \
	trap - EXIT

$(DOC_INPUT_BUILDER_INDEX): $(DOC_INPUT_BUILDERS) $(DOC_INPUT_COVERAGE_REPORT) Makefile
	@printf "$(B_ON)$(FG_MAGENTA)INPUT BUILDERS$(RESET) Updating index\n"
	@mkdir -p "$(dir $@)"
	@set -e; \
	tmp="$@.tmp"; \
	{ \
		printf '<!doctype html>\n'; \
		printf '<html lang="en"><head><meta charset="utf-8">\n'; \
		printf '<meta name="viewport" content="width=device-width, initial-scale=1">\n'; \
		printf '<title>Alamo Input Builders</title>\n'; \
		printf '<meta http-equiv="refresh" content="0; url=$(firstword $(INPUT_BUILDER_EXECUTABLE_NAMES)).html">\n'; \
		printf '<style>body{margin:0;padding:32px;font:14px/1.5 system-ui,sans-serif;color:#262626;background:#f7f7f7}main{max-width:720px;margin:auto;background:white;border:1px solid #d9d9d9;padding:24px}h1{margin-top:0;font-size:22px}ul{padding-left:20px}a{color:#1f77b4}</style>\n'; \
		printf '</head><body><main><h1>Alamo Input Builders</h1>\n'; \
		printf '<p>Select the executable and build configuration whose input structure you need.</p>\n'; \
		printf '<ul>\n'; \
	} > "$$tmp"; \
	for executable in $(INPUT_BUILDER_EXECUTABLE_NAMES); do \
		printf '<li><a href="%s.html">%s</a></li>\n' "$$executable" "$$executable" >> "$$tmp"; \
	done; \
	{ \
		printf '</ul>\n'; \
		printf '<p><a href="../input-schemas/$(notdir $(DOC_INPUT_COVERAGE_REPORT))">Input scrape coverage report</a></p>\n'; \
		printf '</main></body></html>\n'; \
	} >> "$$tmp"; \
	mv "$$tmp" "$@"


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

lib/libalamo-$(POSTFIX).so: ${OBJ} 
	@printf "$(B_ON)$(FG_ORANGE)LIBALAMO$(RESET)             $@\n" 	
	$(QUIET)mkdir -p lib
	$(QUIET)$(CC) -shared -fPIC -o $@ $^

${AMREX_TARGET}/lib/libamrex.so : ${AMREX_TARGET}/lib/libamrex.a
	@printf "$(B_ON)$(FG_ORANGE)LIBAMREX$(RESET)             $@\n" 	
	$(QUIET)$(CC) -shared -fPIC -o $@ -Wl,--whole-archive $< -Wl,--no-whole-archive

docs-with-input-builders: docs

githubpages: docs-with-input-builders cov-report
	mkdir -p ./githubpages/
	echo "<head><meta http-equiv=\"refresh\" content=\"0; url='docs/index.html\" /></head>" > githubpages/index.html
	cp -rf docs/build/html ./githubpages/docs/
	cp -rf cov/ ./githubpages/cov/
	cp -rf $(DOC_INPUT_BUILDER_DIR)/ ./githubpages/inputs/
	cp -rf $(DOC_INPUT_SCHEMA_DIR)/ ./githubpages/input-schemas/

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
