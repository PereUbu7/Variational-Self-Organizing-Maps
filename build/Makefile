RELEASEDIR_LIB?=/usr/lib
RELEASEDIR_HEADER=/usr/include

CXX = g++ 
CXXFLAGS = -Wall -Wpedantic -Werror -fPIC -std=c++20 -I ../extern -I ../include -I /usr/include/
cc = gcc
CCFLAGS = -Wall -fPIC

LIB = -lpthread -ldl
BIN = sqlite libsom test perfTest VSOM

debug: CXXFLAGS += -g -O0
debug: CCFLAGS += -g -O0
debug: all
release: CXXFLAGS += -O3 -msse2
release: CCFLAGS += -O3
release: all

all : $(BIN)

libsom.a : UMatrix.o SqliteDataLoader.o MnistDataLoader.o DataSet.o Som.o SomIndex.o sqlite3.o
	ar -rcs libsom.a UMatrix.o SqliteDataLoader.o MnistDataLoader.o DataSet.o Som.o SomIndex.o sqlite3.o
libsom.so : UMatrix.o SqliteDataLoader.o MnistDataLoader.o DataSet.o Som.o SomIndex.o sqlite3.o
	$(CXX) -shared -o libsom.so UMatrix.o SqliteDataLoader.o MnistDataLoader.o DataSet.o Som.o SomIndex.o sqlite3.o
SqliteDataLoader.o : ../src/SqliteDataLoader.cpp
	$(CXX) $(CXXFLAGS) -c $^
MnistDataLoader.o : ../src/MnistDataLoader.cpp
	$(CXX) $(CXXFLAGS) -c $^
DataSet.o : ../src/DataSet.cpp
	$(CXX) $(CXXFLAGS) -c $^
Som.o : ../src/Som.cpp
	$(CXX) $(CXXFLAGS) -c $^
SomIndex.o : ../src/SomIndex.cpp
	$(CXX) $(CXXFLAGS) -c $^
UMatrix.o : ../src/UMatrix.cpp
	$(CXX) $(CXXFLAGS) -c $^

sqlite3.o : ../extern/sqlite/sqlite3.c
	$(cc) $(CCFLAGS) -o $@ -c $^

test : ../tests/*.cpp libsom.a
	$(CXX) $(CXXFLAGS) $^ -o test $(LIB)

perfTest : ../tests/performance/*.cpp libsom.a
	$(CXX) $(CXXFLAGS) $^ -o perfTest $(LIB)

VSOM : ../apps/main.cpp libsom.a
	$(CXX) $(CXXFLAGS) $^ -o VSOM $(LIB)

libsom : libsom.a libsom.so
sqlite : sqlite3.o

install : 
	mkdir -p $(RELEASEDIR_HEADER)/libsom/
	mkdir -p $(RELEASEDIR_HEADER)/libsom/sqlite/
	cp libsom.so ${RELEASEDIR_LIB}/
	cp libsom.a ${RELEASEDIR_LIB}/
	cp ../include/*.hpp $(RELEASEDIR_HEADER)/libsom/
	cp ../extern/sqlite/sqlite3.h $(RELEASEDIR_HEADER)/libsom/sqlite/

clean :
	rm -f $(BIN) *.o *.a $(BIN)
	rm $(RELEASEDIR_LIB)/libsom.so
	rm $(RELEASEDIR_LIB)/libsom.a
	rm -r $(RELEASEDIR_HEADER)/libsom/

.PHONY: all, clean, release, debug, VSOM, install