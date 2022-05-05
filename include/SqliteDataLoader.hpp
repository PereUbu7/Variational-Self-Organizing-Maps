#pragma once

#include "IDataLoader.hpp"

#include <string>
#include "sqlite/sqlite3.h"

class SqliteDataLoader : public IDataLoader
{
protected:
	sqlite3 *db;
	bool _verbose;

	void close();
	int getMax(char *, std::string);
	int getRecord(char *[], int, std::string) const;

public:
	SqliteDataLoader(bool verbose = false) : _verbose{verbose} {};
	~SqliteDataLoader() override;
	int getElement(char *, int, std::string) override;
	int rows() override;
	bool open(const char *);
	int startTransaction() override;
	int endTransaction() override;
	int minId() override;
	int maxId() override;
	int doesExist(int) override;
};