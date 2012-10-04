CC=g++
CFLAGS=-Wall -O3 -fopenmp -Wno-non-template-friend -funroll-loops

PYTHONVERSION=2.7

# Generic and Debian GNU/Linux
PYHEADERS=/usr/include/python$(PYTHONVERSION)
NUMPYHEADERS=/usr/lib/python$(PYTHONVERSION)/dist-packages/numpy/core/include
EIGEN3=/usr/include/eigen3
PYLIBS=/usr/lib # location of libboost_python.so
BOOSTPYLIBS=/usr/lib

# Brutus
#module load python
#module load boost
## get NUMPYHEADERS location with command numpy.get_include() 
##NUMPYHEADERS=/cluster/apps/python/2.7.2/x86_64/lib64/python2.7/site-packages/numpy/core/include
#PYHEADERS=/cluster/apps/python/2.7.2/x86_64/include/python2.7/
#EIGEN3=eigen3 # local installation
#BOOSTPYLIBS=/cluster/apps/boost/1_47_0_nompi/x86_64/gcc_4.1.2/lib64
#PYLIBS=/cluster/apps/python/2.7.2/x86_64/lib64

# OS X
#PYHEADERS=/usr/include/python$(PYTHONVERSION)
#PYHEADERS2=/Library/Frameworks/Python.framework/Versions/6.3/include/python2.6
#BOOSTPYLIBS=/usr/lib # location of libboost_python.so
#PYLIBS=/usr/lib # ?? not needed ??

modulename=HyperCubicShape
modulename2=EvaluateBasis

all: python evaluate prog
 
python: ${modulename}.cpp
	$(CC) $(CFLAGS) -c -fPIC -DPYTHONMODULE -I$(PYHEADERS) -I$(EIGEN3) ${modulename}.cpp
	$(CC) $(CFLAGS) -shared -o ${modulename}.so ${modulename}.o $(BOOSTPYLIBS)/libboost_python.so -I$(EIGEN3) -L$(PYLIBS) -lpython$(PYTHONVERSION)

prog: ${modulename}.cpp
	$(CC) $(CFLAGS) -c -I$(PYHEADERS) -I$(EIGEN3) ${modulename}.cpp
	$(CC) $(CFLAGS) -o ${modulename} ${modulename}.o $(BOOSTPYLIBS)/libboost_python.so -I$(EIGEN3) -L$(PYLIBS) -lpython$(PYTHONVERSION)

evaluate: ${modulename2}.cpp
	$(CC) $(CFLAGS) -c -fPIC -DPYTHONMODULE -I$(PYHEADERS) -I$(NUMPYHEADERS) -I$(EIGEN3) ${modulename2}.cpp
	$(CC) $(CFLAGS) -shared -o ${modulename2}.so ${modulename2}.o $(BOOSTPYLIBS)/libboost_python.so -I$(EIGEN3) -L$(PYLIBS) -lpython$(PYTHONVERSION)

install:
	cp ${modulename2}.so ../../src/

clean:
	$(RM) ${modulename} ${modulename}.so *.o *.gch
	$(RM) ${modulename2} ${modulename2}.so *.o
	$(RM) *.pyc a.out


