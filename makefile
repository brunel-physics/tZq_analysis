LIBRARY_SOURCES = $(wildcard src/common/*.cpp)
LIBRARY_OBJECT_FILES = $(patsubst src/common/%.cpp,obj/%.o,${LIBRARY_SOURCES})
LIBRARY = lib/libTQZanalysisTools.so

EXECUTABLE_SOURCES = $(wildcard src/common/*.cxx)
EXECUTABLE_OBJECT_FILES = $(patsubst src/common/%.cxx,obj/%.o,${EXECUTABLE_SOURCES})
EXECUTABLES = $(patsubst src/common/%.cxx,bin/%,${EXECUTABLE_SOURCES})

LIBRARY_PATH = 	-L$(shell root-config --libdir) \
		-Llib \
		-L${HOME}/usr/lib \
		-L${HOME}/usr/lib64 \

LIBRARIES = $(shell root-config --libs) \
		-lconfig++ \
		-lLHAPDF \
		-lz \

INCLUDE_PATH = 	-Iinclude  \
		-I/usr/include \
		-I${HOME}/usr/include \

CFLAGS = $(shell root-config --cflags) -g -O2 -pipe -Wall -W -Woverloaded-virtual -MMD -MP -fPIC ${INCLUDE_PATH}
LDFLAGS = -Wl,-R${PWD}/lib,-R$(shell root-config --libdir)

#LHAP = -I/cms/cmssw/slc6_amd64_gcc493/external/lhapdf/6.1.5-kpegke3/include/LHAPDF -L/cms/cmssw/slc6_amd64_gcc493/external/lhapdf/6.1.5-kpegke3/lib -lLHAPDF
#LHAPDFLAGS = -I$(shell cd ${CMSSW_BASE}; scram tool tag lhapdffull INCLUDE) -L$(shell cd ${CMSSW_BASE}; scram tool tag lhapdffull LIBDIR) -lLHAPDF -lgfortran -lz

LINK_LIBRARY_FLAGS = -shared -Wall -g -O0 -rdynamic ${LIBRARY_PATH} ${LIBRARIES}
LINK_EXECUTABLE_FLAGS = -Wall -g -O0 -rdynamic ${LIBRARY_PATH} ${LIBRARIES} ${LDFLAGS} -lTQZanalysisTools

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


${EXECUTABLES}: bin/%: obj/%.o ${EXECUTABLE_OBJECT_FILES}
	g++ ${LINK_EXECUTABLE_FLAGS} $(LHAP) $< -o $@

${EXECUTABLE_OBJECT_FILES}: obj/%.o : src/common/%.cxx
	mkdir -p {bin,obj,lib}
	g++ -c ${CFLAGS} $< -o $@

-include $(EXECUTABLE_OBJECT_FILES:.o=.d)
