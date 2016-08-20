
include make.inc

all: lib test examples

lib:
	mkdir -p ./lib
	mkdir -p ./objs
	( cd ./src; $(MAKE) all )

test:
	( cd ./test; $(MAKE) clean; $(MAKE) run )

examples:
	( cd ./examples; $(MAKE) clean; $(MAKE) run )

install:
	cp ./lib/lib$(LIB_NAME).a $(EXTERNAL_LIBRARY_PATH)
	cp -r ../$(LIB_NAME) $(BIN_PATH)

clean: 
	(cd ./src; $(MAKE) clean; cd ../test; $(MAKE) clean)

.PHONY: all lib test examples install
