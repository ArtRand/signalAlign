rootPath = ../
include ../include.mk

all: 
	cd lastz-distrib-1.03.54 && make install
	mv ${binPath}/lastz ${binPath}/cPecanLastz
	mv ${binPath}/lastz_D ${binPath}/cPecanLastz_D

clean: 
	cd lastz-distrib-1.03.54 && make clean 
	rm -rf ${binPath}/cPecanLastz ${binPath}/cPecanLastz_D
 
