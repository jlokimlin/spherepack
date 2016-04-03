
include make.inc

all: libspherepack lib testlib

libspherepack:
	mkdir -p ./lib
	mkdir -p ./objs
	( cd ./spherepack_src; $(MAKE) clean; $(MAKE) )

lib: 
	( cd ./src; $(MAKE) all )

testlib:
	( cd ./test; $(MAKE) clean; $(MAKE) run )

install:
	cp ./lib/lib$(LIB_NAME).a $(EXTERNAL_LIBRARY_PATH)
	cp -r ../$(LIB_NAME) $(BIN_PATH)

clean: 
	( cd ./src; $(MAKE) clean; cd ../test; $(MAKE) clean )

.PHONY: all lipspherepack lib testlib install