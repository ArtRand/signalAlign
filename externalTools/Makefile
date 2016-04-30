rootPath = ../
include ../include.mk

all: 
	cd lastz-distrib-1.03.54 && make install
	mv ${signalAlignBin}/lastz ${signalAlignBin}/cPecanLastz
	mv ${signalAlignBin}/lastz_D ${signalAlignBin}/cPecanLastz_D

clean: 
	cd lastz-distrib-1.03.54 && make clean 
	rm -rf ${binPath}/cPecanLastz ${binPath}/cPecanLastz_D
 
