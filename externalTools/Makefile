rootPath = ../
include ../include.mk

all: 
	cd lastz-distrib-1.03.54 && make install
	mv ./bin/lastz ../bin/cPecanLastz
	mv ./bin/lastz_D ../bin/cPecanLastz_D
	rmdir ./bin

clean: 
	cd lastz-distrib-1.03.54 && make clean 
	rm -rf ${signalAlignBin}/cPecanLastz ${signalAlignBin}/cPecanLastz_D
 
