
MAKE = /usr/bin/make

AR = /usr/bin/ar

#EXTERNAL_LIBRARY_PATH = /usr/local/lib
EXTERNAL_LIBRARY_PATH = /usr/local/lib64

BIN_PATH = /usr/local/bin

OBJ_DIRECTORY = ../objs

LIB_DIRECTORY = ../lib

MOD_DIRECTORY = $(LIB_DIRECTORY)

LIB_NAME = spherepack

LIB = $(LIB_DIRECTORY)/lib$(LIB_NAME).a

OBJ := $(patsubst %.f90, %.o, $(wildcard *.f90))

OBJS = $(addprefix $(OBJ_DIRECTORY)/, $(OBJ))

FC = gfortran-6

FFLAGS = \
-O3 \
$(ENFORCE_64_BIT_PRECISION) \
$(LANGUAGE_OPTIONS)

ENFORCE_64_BIT_PRECISION = \
-fdefault-real-8 \
-fdefault-double-8

LANGUAGE_OPTIONS = \
-fimplicit-none

ERROR_WARNING_OPTIONS = \
-Wall

