
# Set directories

#EXTERNAL_LIBRARY_PATH = /usr/local/lib
EXTERNAL_LIBRARY_PATH = /usr/local/lib64

BIN_PATH = /usr/local/bin

SRC_DIRECTORY = ./

SOURCES = $(wildcard $(SRC_DIRECTORY)/*.f90)

SOURCES_WITH_PREFIX = $(addprefix $(SRC_DIRECTORY)/, $(SOURCES))

OBJECTS := $(patsubst %.f90, %.o, $(wildcard *.f90))

OBJ_DIRECTORY = ../objs

LIB_DIRECTORY = ../lib

OBJECTS_WITH_PREFIX = $(addprefix $(OBJ_DIRECTORY)/, $(OBJECTS))

# Set the compiler options

FC = gfortran

FFLAGS = \
-O3 \
$(LANGUAGE_OPTIONS) \
$(ERROR_WARNING_OPTIONS) \
$(CODE_GENERATION_OPTIONS)

LANGUAGE_OPTIONS = \
-fimplicit-none

ERROR_WARNING_OPTIONS = \
-Wall

CODE_GENERATION_OPTIONS = \
-fcheck=all \
-fbounds-check \
-fno-backtrace \
-fno-unsafe-math-optimizations \
-frounding-math \
-fsignaling-nans

