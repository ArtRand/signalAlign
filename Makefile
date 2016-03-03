rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

signalAlignDependencies =  ${basicLibsDependencies}
signalAlignLib = ${basicLibs}

all : ${libPath}/signalAlignLib.a ${binPath}/signalAlignLibTests ${binPath}/compareDistributions \
      ${binPath}/signalMachine ${binPath}/runSignalAlign ${sonLibrootPath}/signalAlignLib.py \
      ${binPath}/buildHdpUtil ${binPath}/trainModels ${binPath}/hdp_pipeline ${binPath}/testSignalAlign
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

${binPath}/buildHdpUtil : buildHdpUtil.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${binPath}/buildHdpUtil buildHdpUtil.c ${libPath}/signalAlignLib.a ${signalAlignLib}

${binPath}/runSignalAlign : ${rootPath}scripts/runSignalAlign.py
	cp ${rootPath}scripts/runSignalAlign.py ${binPath}/runSignalAlign
	chmod +x ${binPath}/runSignalAlign

${binPath}/trainModels : ${rootPath}scripts/trainModels.py
	cp ${rootPath}scripts/trainModels.py ${binPath}/trainModels
	chmod +x ${binPath}/trainModels

${binPath}/hdp_pipeline : ${rootPath}scripts/hdp_pipeline.py
	cp ${rootPath}scripts/hdp_pipeline.py ${binPath}/hdp_pipeline
	chmod +x ${binPath}/hdp_pipeline

${binPath}/testSignalAlign : ${rootPath}scripts/testSignalAlign.py
	cp ${rootPath}scripts/testSignalAlign.py ${binPath}/testSignalAlign
	chmod +x ${binPath}/testSignalAlign


${sonLibrootPath}/signalAlignLib.py : ${rootPath}scripts/signalAlignLib.py
	cp ${rootPath}scripts/signalAlignLib.py ${sonLibRootPath}/signalAlignLib.py

${libPath}/signalAlignLib.a : ${libSources} ${libHeaders} ${stBarDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc signalAlignLib.a *.o
	ranlib signalAlignLib.a
	rm *.o
	mv signalAlignLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/
