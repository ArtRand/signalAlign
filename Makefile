rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

signalAlignDependencies =  ${basicLibsDependencies}
signalAlignLib = ${basicLibs}

all : ${libPath}/signalAlignLib.a ${binPath}/signalAlignLibTests ${binPath}/compareDistributions ${binPath}/signalMachine ${binPath}/runSignalAlign ${sonLibrootPath}/signalAlignLib.py
	# disabled right now so that we don't build Lastz every time I do an update
	#cd externalTools && make all
	
clean : 
	rm -f ${binPath}/cPecanRealign ${binPath}/cPecanEm ${binPath}/cPecanLibTests  ${libPath}/cPecanLib.a
	cd externalTools && make clean

${binPath}/compareDistributions : compareDistributions.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/compareDistributions compareDistributions.c ${libPath}/signalAlignLib.a ${signalAlignLib}

${binPath}/signalAlignLibTests : ${libTests} tests/*.h ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${binPath}/signalAlignLibTests ${libTests} ${libPath}/signalAlignLib.a ${signalAlignLib}

${binPath}/signalMachine : signalMachine.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/signalMachine signalMachine.c ${libPath}/signalAlignLib.a ${signalAlignLib}

${binPath}/runSignalAlign : ${rootPath}scripts/runSignalAlign.py
	cp ${rootPath}scripts/runSignalAlign.py ${binPath}/runSignalAlign
	chmod +x ${binPath}/runSignalAlign

${sonLibrootPath}/signalAlignLib.py : ${rootPath}scripts/signalAlignLib.py
	cp ${rootPath}scripts/signalAlignLib.py ${sonLibRootPath}/signalAlignLib.py

${libPath}/signalAlignLib.a : ${libSources} ${libHeaders} ${stBarDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc signalAlignLib.a *.o
	ranlib signalAlignLib.a
	rm *.o
	mv signalAlignLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/
