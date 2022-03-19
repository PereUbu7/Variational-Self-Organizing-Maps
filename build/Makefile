CXX = g++ 
CXXFLAGS = -Wall -std=c++20 -I ../extern -I ../include -I /usr/include/
cc = gcc

LIB = -lpthread -ldl
BIN = sqlite libsom test perfTest

debug: CXXFLAGS += -g -O0
debug: all
release: CXXFLAGS += -O3
release: all

all : $(BIN)

libsom.a : DataBase.o DataSet.o Som.o SomIndex.o sqlite3.o
	ar -rcs libsom.a DataBase.o DataSet.o Som.o SomIndex.o sqlite3.o
DataBase.o : ../src/DataBase.cpp
	$(CXX) $(CXXFLAGS) -c $^
DataSet.o : ../src/DataSet.cpp
	$(CXX) $(CXXFLAGS) -c $^
Som.o : ../src/Som.cpp
	$(CXX) $(CXXFLAGS) -c $^
SomIndex.o : ../src/SomIndex.cpp
	$(CXX) $(CXXFLAGS) -c $^

sqlite3.o : ../extern/sqlite/sqlite3.c
	$(cc) -o $@ -c $^

test.a : ../tests/*.cpp libsom.a
	$(CXX) $(CXXFLAGS) $^ -o test $(LIB)

perfTest.a : ../tests/performance/*.cpp libsom.a
	$(CXX) $(CXXFLAGS) $^ -o perfTest $(LIB)

libsom : libsom.a
sqlite : sqlite3.o
test : test.a
perfTest : perfTest.a

clean :
	rm -f $(BIN) *.o *.a $(BIN)

.PHONY: all, clean, release, debug