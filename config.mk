CROW_DIR ?= $(ROOT_DIR)/crow

SWIG_VERSION := $(shell swig -version 2>/dev/null)
PYTHON_CONFIG_WHICH := $(shell which python-config 2>/dev/null)
ifeq ($(CROW_USE_PYTHON3),TRUE)
	PYTHON3_CONFIG_WHICH := $(shell which python3-config 2>/dev/null)
endif

UNAME := $(shell uname)

ifneq ($(PYTHON_CONFIG_WHICH),)
	PYTHON2_INCLUDE=$(shell python-config --includes)
	PYTHON2_LIB=$(shell python-config --ldflags)
endif

ifneq ($(PYTHON3_CONFIG_WHICH),)
	PYTHON3_INCLUDE = $(shell python3-config --includes) #-DPy_LIMITED_API
	PYTHON3_LIB = $(shell python3-config --ldflags) #-DPy_LIMITED_API
endif

ifneq ($(PYTHON3_CONFIG_WHICH),)
ifeq ($(findstring SWIG Version 2,$(SWIG_VERSION)),)
	CONTROL_MODULES =
	PYTHON_MODULES =
else
	SWIG_PY_FLAGS=-py3
	CONTROL_MODULES = $(CROW_DIR)/control_modules/_distribution1D.so $(CROW_DIR)/control_modules/_crowtools.so
	PYTHON_MODULES = $(CROW_DIR)/python_modules/_distribution1Dpy2.so $(CROW_DIR)/python_modules/_distribution1Dpy3.so $(CROW_DIR)/python_modules/_interpolationNDpy2.so $(CROW_DIR)/python_modules/_interpolationNDpy3.so
endif #Have SWIG
	PYTHON_INCLUDE=$(PYTHON3_INCLUDE)
	PYTHON_LIB=$(PYTHON3_LIB)

else #no Python3

ifneq ($(PYTHON_CONFIG_WHICH),)
#Python 2 config found but not Python 3
	PYTHON_INCLUDE=$(PYTHON2_INCLUDE)
	PYTHON_LIB=$(PYTHON2_LIB)
	#CONTROL_MODULES=
	SWIG_PY_FLAGS=
	PYTHON_MODULES = $(CROW_DIR)/python_modules/_distribution1Dpy2.so $(CROW_DIR)/python_modules/_interpolationNDpy2.so
	CONTROL_MODULES=$(CROW_DIR)/control_modules/_distribution1D.so $(CROW_DIR)/control_modules/_crowtools.so
else #No python3 config or python2 config
	PYTHON_INCLUDE = -DNO_PYTHON_FOR_YOU
	PYTHON_LIB = -DNO_PYTHON_FOR_YOU
	CONTROL_MODULES =
	PYTHON_MODULES =
endif

endif

CROW_LIB_INCLUDE_DIR := $(CROW_DIR)/contrib/include

ifeq  ($(UNAME),Darwin)
crow_shared_ext := dylib
else
crow_shared_ext := so
endif

HAS_DYNAMIC := $(shell $(libmesh_LIBTOOL) --config | grep build_libtool_libs | cut -d'=' -f2 )

ifeq ($(HAS_DYNAMIC),no)
ifdef CONTROL_MODULES
$(warning CROW modules must be compiled with shared libmesh libraries)
	CONTROL_MODULES =
endif
endif
