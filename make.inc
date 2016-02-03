# Set the compiler options

FC = gfortran

FFLAGS = \
-O3 \
$(LANGUAGE_OPTIONS) \
$(ERROR_WARNING_OPTIONS) \
$(CODE_GENERATION_OPTIONS)

LANGUAGE_OPTIONS = \
-fimplicit-none \
-fdefault-real-8 \
-fdefault-double-8

ERROR_WARNING_OPTIONS = \
-Wall

CODE_GENERATION_OPTIONS = \
-fcheck=all \
-fbounds-check \
-fno-backtrace \
-fno-unsafe-math-optimizations \
-frounding-math \
-fsignaling-nans

# Set directories
SRC_DIRECTORY = ./

SOURCES = $(wildcard $(SRC_DIRECTORY)/*.f90)

SOURCES_WITH_PREFIX = $(addprefix $(SRC_DIRECTORY)/, $(SOURCES))

OBJECTS := $(patsubst %.f90, %.o, $(wildcard *.f90))

OBJ_DIRECTORY = ../obj

LIB_DIRECTORY = ../lib

OBJECTS_WITH_PREFIX = $(addprefix $(OBJ_DIRECTORY)/, $(OBJECTS))
