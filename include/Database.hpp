#pragma once

#include <string>
#include "sqlite/sqlite3.h"

class DataBase
{
	protected:
		sqlite3 *db;
        bool _verbose;
	public:
		DataBase(bool verbose = false) : _verbose{verbose} {};
		~DataBase();
		bool open(const char*);
		void close();
		int getElement(char*, int, std::string);
		int getMax(char*, std::string);
		int getRecord(char*[], int, std::string) const;
		int rows();
		int minId();
		int maxId();
		int doesExist(int);
		int startTransaction();
		int endTransaction();
};