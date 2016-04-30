rootPath = ../
include ../include.mk

all: 
	cd lastz-distrib-1.03.54 && make install
	mv ./bin/lastz ../bin/cPecanLastz
	mv ./bin/lastz_D ../bin/cPecanLastz_D

clean: 
	cd lastz-distrib-1.03.54 && make clean 
	rm -rf ${binPath}/cPecanLastz ${binPath}/cPecanLastz_D
 
