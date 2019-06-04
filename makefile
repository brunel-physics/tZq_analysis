CXX ?= g++

LIBRARY_SOURCES = $(wildcard src/*.cpp)
LIBRARY_OBJECT_FILES = $(patsubst src/%.cpp,obj/%.o,${LIBRARY_SOURCES})
LIBRARY = lib/libTQZanalysisTools.so

EXECUTABLE_SOURCES = $(wildcard src/*.cxx)
EXECUTABLE_OBJECT_FILES = $(patsubst src/%.cxx,obj/%.o,${EXECUTABLE_SOURCES})
EXECUTABLES = $(patsubst src/%.cxx,bin/%.exe,${EXECUTABLE_SOURCES})

LIBRARY_PATH = 	-L$(shell root-config --libdir) \
		-Llib \
		-L/scratch/shared/sw/yaml-cpp/0.6.2/x86_64-slc6-gcc82-opt/lib \
		-L/cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-slc6-gcc8-opt/lib

LIBRARIES = 	$(shell root-config --libs) \
		-lLHAPDF \
		-lTMVA  \
		-lz \
		-lboost_system \
		-lboost_filesystem \
		-lboost_program_options \
		-lyaml-cpp \

INCLUDE_PATH = 	-Iinclude  \
		-isystem/cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-slc6-gcc8-opt/include \
		-isystem/scratch/shared/sw/yaml-cpp/0.6.2/x86_64-slc6-gcc82-opt/include \
		-isystem/scratch/shared/include \
		-isystem$(shell root-config --incdir) \

CFLAGS = ${INCLUDE_PATH} -std=c++17 -MMD -MP -std=c++17 -march=native \
                         -mtune=native -pipe -O3 -fPIC -m64 -pthread

ifeq ($(CXX),g++)
  CFLAGS += -Wall -Wextra -Wpedantic -Wcast-align -Wcast-qual \
			-Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self \
			-Wlogical-op -Wmissing-declarations -Wmissing-include-dirs \
			-Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls \
			-Wshadow -Wsign-promo -Wstrict-null-sentinel \
			-Wstrict-overflow=5 -Wswitch-default -Wundef \
			-Wuseless-cast -Wzero-as-null-pointer-constant -Wduplicated-cond \
			-Wduplicated-branches -Wrestrict -Wnull-dereference -Wswitch-enum \
			-Wswitch-bool -Wswitch-unreachable -Wno-coverage-mismatch \
			-Wimplicit-fallthrough=5 -Wsync-nand -Wunknown-pragmas \
			-Wstringop-overflow=4 -Wstringop-truncation -Wsuggest-final-types \
			-Wsuggest-final-methods -Wsuggest-override -Walloc-zero -Walloca \
			-Wtrampolines -Wfloat-equal -Wunsafe-loop-optimizations \
			-Wplacement-new=2 -Wunused-macros -Wconditionally-supported \
			-Wsubobject-linkage -Wdate-time -Wextra-semi \
			-Wno-aggressive-loop-optimizations -Wpacked -Winvalid-pch  \
			-Wvector-operation-performance -Wdisabled-optimization \
			-Wstack-protector -Whsa -Wsuggest-attribute=const \
			-Wsuggest-attribute=pure -Wsuggest-attribute=noreturn \
            -Wsuggest-attribute=format -Wsuggest-attribute=cold
  			# -Winline -Wconversion -Wsign-conversion
else ifeq ($(CXX),clang++)
  CFLAGS += -Weverything -Wno-c++98-compat -Wno-double-promotion \
		    -Wno-covered-switch-default
else
  CFLAGS += $(UNKNOWN_CXXFLAGS)
endif

LINK_LIBRARY_FLAGS = -shared -g -O2 -rdynamic ${LIBRARY_PATH} ${LIBRARIES}
LINK_EXECUTABLE_FLAGS = -g -O2 -rdynamic ${LIBRARY_PATH} ${LIBRARIES} \
			-lTQZanalysisTools \
			-Wl,-R/scratch/shared/lib,-Rlib,-R../lib,-R${PWD}/lib,-R/scratch/shared/sw/yaml-cpp/0.6.2/x86_64-slc6-gcc82-opt/lib,--enable-new-dtags

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
	${CXX} ${LINK_LIBRARY_FLAGS} ${LIBRARY_OBJECT_FILES} -o $@

${LIBRARY_OBJECT_FILES}: obj/%.o : src/%.cpp
	mkdir -p {bin,obj,lib}
	${CXX} -c ${CFLAGS}  $< -o $@

-include $(LIBRARY_OBJECT_FILES:.o=.d)


${EXECUTABLES}: bin/%.exe: obj/%.o ${EXECUTABLE_OBJECT_FILES}
	${CXX} ${LINK_EXECUTABLE_FLAGS} $< -o $@

${EXECUTABLE_OBJECT_FILES}: obj/%.o : src/%.cxx
	mkdir -p {bin,obj,lib}
	${CXX} -c ${CFLAGS} $< -o $@

-include $(EXECUTABLE_OBJECT_FILES:.o=.d)
