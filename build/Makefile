CXX = g++ 
CXXFLAGS = -Wall -Wpedantic -Werror -std=c++20 -I ../extern -I ../include -I /usr/include/
cc = gcc
CCFLAGS = -Wall

LIB = -lpthread -ldl
BIN = sqlite libsom test perfTest VSOM

debug: CXXFLAGS += -g -O0
debug: CCFLAGS += -g -O0
debug: all
release: CXXFLAGS += -O3
release: CCFLAGS += -O3
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

test : ../tests/*.cpp libsom.a
	$(CXX) $(CXXFLAGS) $^ -o test $(LIB)

perfTest : ../tests/performance/*.cpp libsom.a
	$(CXX) $(CXXFLAGS) $^ -o perfTest $(LIB)

VSOM : ../apps/main.cpp libsom.a
	$(CXX) $(CXXFLAGS) $^ -o VSOM $(LIB)

libsom : libsom.a
sqlite : sqlite3.o

clean :
	rm -f $(BIN) *.o *.a $(BIN)

.PHONY: all, clean, release, debug, VSOM