CXX = g++ -I ../extern -I ../include -I /usr/include/
cc = gcc

LIB = -lpthread -ldl
BIN = sqlite libsom test

all : $(BIN)

libsom.a : DataBase.o DataSet.o Som.o SomIndex.o sqlite3.o
	ar -rcs libsom.a DataBase.o DataSet.o Som.o SomIndex.o sqlite3.o
DataBase.o : ../src/DataBase.cpp
	$(CXX) -c $^
DataSet.o : ../src/DataSet.cpp
	$(CXX) -c $^
Som.o : ../src/Som.cpp
	$(CXX) -c $^
SomIndex.o : ../src/SomIndex.cpp
	$(CXX) -c $^

sqlite3.o : ../extern/sqlite/sqlite3.c
	$(cc) -o $@ -c $^

test.a : ../tests/*.cpp libsom.a
	$(CXX) $^ -o test $(LIB)

libsom : libsom.a
sqlite : sqlite3.o
test : test.a

clean :
	rm -f $(BIN) *.o *.a $(BIN)

.PHONY: all, clean