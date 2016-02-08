
include make.inc

#EXTERNAL_LIBRARY_PATH = /usr/local/lib
EXTERNAL_LIBRARY_PATH = /usr/local/lib64

BIN_PATH = /usr/local/bin

all: lib testlib

lib: 
	mkdir -p ./lib
	mkdir -p ./objs
	cd ./src; make all

testlib:
	cd ./test; make run

install:
	cp ./lib/libspherepack_wrapper.a $(EXTERNAL_LIBRARY_PATH)
	cp -r ../spherepack_wrapper $(BIN_PATH)

.PHONY: all lib testlib install