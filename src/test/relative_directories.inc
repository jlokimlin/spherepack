SRC_DIRECTORY = ./

SOURCES = $(wildcard $(SRC_DIRECTORY)/*.f90)

SOURCES_WITH_PREFIX = $(addprefix $(SRC_DIRECTORY)/, $(SOURCES))

OBJECTS := $(patsubst %.f90, %.o, $(wildcard *.f90))

OBJ_DIRECTORY = ../../obj

LIB_DIRECTORY = ../../lib

OBJECTS_WITH_PREFIX = $(addprefix $(OBJ_DIRECTORY)/, $(OBJECTS))
