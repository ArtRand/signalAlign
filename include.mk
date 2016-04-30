#Modify this variable to set the location of sonLib
sonLibRootPath=${rootPath}./sonLib
sonLibPath=${sonLibRootPath}/lib
#Use sonLib bin and lib dirs
signalAlignBin = ${rootPath}/bin
binPath=${sonLibRootPath}/bin
libPath=${sonLibPath}

include  ${sonLibRootPath}/include.mk

basicLibs = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a ${dblibs}
basicLibsDependencies = ${sonLibPath}/sonLib.a ${sonLibPath}/cuTest.a 
