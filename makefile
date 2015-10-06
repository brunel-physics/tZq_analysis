LIBRARY_SOURCES = $(wildcard src/common/*.cpp)
LIBRARY_OBJECT_FILES = $(patsubst src/common/%.cpp,obj/%.o,${LIBRARY_SOURCES})
LIBRARY = lib/libTQZanalysisTools.so

EXECUTABLE_SOURCES = $(wildcard src/common/*.cxx)
EXECUTABLE_OBJECT_FILES = $(patsubst src/common/%.cxx,obj/%.o,${EXECUTABLE_SOURCES})
EXECUTABLES = $(patsubst src/common/%.cxx,bin/%.exe,${EXECUTABLE_SOURCES})

LIBRARY_PATH = 	-L/home/eepgadm/root/lib \
		-Llib \
		-L/home/eepgadm/lib/local/lib\

LIBRARIES = 	-lCore \
		-lCint  \
		-lRIO  \
		-lNet  \
		-lHist  \
		-lGraf  \
		-lGraf3d  \
		-lGpad  \
		-lTMVA  \
		-lTree  \
		-lRint  \
		-lPostscript  \
		-lMatrix  \
		-lPhysics  \
		-lMathCore  \
		-lThread  \
		-pthread  \
		-lm  \
		-ldl \
		-lconfig++ \
		-lLHAPDF \
		-lz \

INCLUDE_PATH = 	-Iinclude  \
		-I/home/eepgadm/root/include \
		-I/usr/include \
		-I/home/eepgadm/lib/local/include

CFLAGS = -g -O2 -pipe -Wall -W -Woverloaded-virtual -MMD -MP -fPIC -pthread -std=c++0x $(shell root-config --cflags) ${INCLUDE_PATH}

LHAP = -I/cms/cmssw/slc6_amd64_gcc472/external/lhapdf/5.9.1-cms/full/include -L/cms/cmssw/slc6_amd64_gcc472/external/lhapdf/5.9.1-cms/full/lib -lLHAPDF
#LHAPDFLAGS = -I$(shell cd ${CMSSW_BASE}; scram tool tag lhapdffull INCLUDE) -L$(shell cd ${CMSSW_BASE}; scram tool tag lhapdffull LIBDIR) -lLHAPDF -lgfortran -lz

ROOTSYS = /home/eepgadm/root/

LINK_LIBRARY_FLAGS = -shared -Wall -g -O0 -rdynamic ${LIBRARY_PATH} ${LIBRARIES}
LINK_EXECUTABLE_FLAGS = -Wall -g -O0 -rdynamic ${LIBRARY_PATH} ${LIBRARIES} -lTQZanalysisTools

.PHONY: all _all clean _cleanall build _buildall install _installall rpm _rpmall test _testall spec_update

default: build

clean: _cleanall
_cleanall:
	rm -rf obj
	rm -rf bin
	rm -rf lib

all: _all
build: _all
buildall: _all
_all: ${LIBRARY} ${EXECUTABLES}

${LIBRARY}: ${LIBRARY_OBJECT_FILES}
	g++ ${LINK_LIBRARY_FLAGS} ${LIBRARY_OBJECT_FILES} -o $@

${LIBRARY_OBJECT_FILES}: obj/%.o : src/common/%.cpp
	mkdir -p {bin,obj,lib}
	g++ -c ${CFLAGS}  $< -o $@

-include $(LIBRARY_OBJECT_FILES:.o=.d)


${EXECUTABLES}: bin/%.exe: obj/%.o ${EXECUTABLE_OBJECT_FILES}
	g++ ${LINK_EXECUTABLE_FLAGS} $(LHAP) $< -o $@

${EXECUTABLE_OBJECT_FILES}: obj/%.o : src/common/%.cxx
	mkdir -p {bin,obj,lib}
	g++ -c ${CFLAGS} $< -o $@

-include $(EXECUTABLE_OBJECT_FILES:.o=.d)
