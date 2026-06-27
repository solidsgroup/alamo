# CUDA builds compile one GPU-clean integrator closure at a time.
#
# Add a new integrator here only after its entry point and dependency closure are
# nvcc-clean. The top-level Makefile consumes ALAMO_GPU_MAIN and
# ALAMO_GPU_SOURCES; it should not carry integrator-specific source lists.

ALAMO_GPU_SUPPORTED_INTEGRATORS := flame
ALAMO_GPU_INTEGRATOR ?= flame

ifneq ($(filter $(ALAMO_GPU_INTEGRATOR),$(ALAMO_GPU_SUPPORTED_INTEGRATORS)),$(ALAMO_GPU_INTEGRATOR))
$(error Unsupported ALAMO_GPU_INTEGRATOR '$(ALAMO_GPU_INTEGRATOR)'; supported: $(ALAMO_GPU_SUPPORTED_INTEGRATORS))
endif

ALAMO_GPU_MAIN_flame := src/alamo_gpu.cc
ALAMO_GPU_SOURCES_flame := \
    src/BC/BC.cpp \
    src/BC/Constant.cpp \
    src/IO/FileNameParse.cpp \
    src/IO/ParmParse.cpp \
    src/IO/WriteMetaData.cpp \
    src/Integrator/Flame.cpp \
    src/Integrator/Integrator.cpp \
    src/Operator/Elastic.cpp \
    src/Operator/Operator.cpp \
    src/Set/Set.cpp \
    src/Util/Util.cpp

ALAMO_GPU_MAIN := $(ALAMO_GPU_MAIN_$(ALAMO_GPU_INTEGRATOR))
ALAMO_GPU_SOURCES := $(ALAMO_GPU_SOURCES_$(ALAMO_GPU_INTEGRATOR))
