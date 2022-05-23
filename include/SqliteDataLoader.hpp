#pragma once

#include "IDataLoader.hpp"

#include <string>
#include <vector>
#include "sqlite/sqlite3.h"


class SqliteDataLoader : public IDataLoader
{
protected:
	std::vector<std::string> _columnNames;
	std::vector<float> _weight;
	std::vector<int> _isBinary;
	std::vector<int> _isContinuous;
	sqlite3 *db;
	bool _verbose;
	bool hasOpenTransaction, hasOpenDatabase;
	size_t totalNumberOfRows;
	size_t vectorLength;

	void close();
	int getMax(char *, std::string);
	int getRecord(char *[], int, std::string) const;
	int startTransaction();
	int endTransaction();
	int minId();
	int maxId();
	int doesExist(int);
	int getElement(char *, int, std::string);
	size_t numberOfRows();
	void loadColumnSpecData(const std::string &path);

public:
	SqliteDataLoader(std::vector<std::string> columnNames, std::vector<float> weights, std::vector<int> isBinary, bool verbose = false) : 
		_columnNames{columnNames},
		_weight{weights},
		_isBinary{isBinary},
		// TODO: invert _isContinuous
		_isContinuous{isBinary},
		_verbose{verbose}, 
		hasOpenTransaction{false},
		hasOpenDatabase{false},
		totalNumberOfRows{},
		vectorLength{columnNames.size()}
	{};
	SqliteDataLoader(const std::string &specFilePath, bool verbose = false) :
		_columnNames{},
		_weight{},
		_isBinary{},
		_verbose{verbose}, 
		hasOpenTransaction{false},
		hasOpenDatabase{false},
		totalNumberOfRows{},
		vectorLength{0}
	{
		loadColumnSpecData(specFilePath);
	};
	~SqliteDataLoader() override;
	size_t load(std::optional<size_t> maxCount) override;
	bool open(const char *dbPath) override;
	const std::vector<float> &getWeights() const noexcept override;
	const std::vector<int> &getBinary() const noexcept override;
	const std::vector<int> &getContinuous() const noexcept override;
	std::string getName(size_t index) const noexcept;
	const std::vector<std::string> &getNames() const noexcept override;
	size_t getDepth() const noexcept;
};