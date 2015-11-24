src_directory = ./

sources = $(wildcard $(src_directory)/*.f90)

sources_with_prefix = $(addprefix $(src_directory)/, $(sources))

objects := $(patsubst %.f90, %.o, $(wildcard *.f90))

obj_directory = ../obj

lib_directory = ../lib

objects_with_prefix = $(addprefix $(obj_directory)/, $(objects))
