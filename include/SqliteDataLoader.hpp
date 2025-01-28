#pragma once

#include "IDataLoader.hpp"
#include "ColumnSpec.hpp"

#include <string>
#include <vector>
#include "sqlite/sqlite3.h"


class SqliteDataLoader : public IDataLoader
{
public:
	SqliteDataLoader(const std::string &specFilePath, std::optional<size_t> maxLoadCount = std::nullopt, bool verbose = false) :
		IDataLoader{maxLoadCount},
		_columnSpec{},
		_tableNames{},
		_columnNames{},
		_weight{},
		_isBinary{},
		_tableName{},
		_verbose{verbose}, 
		hasOpenTransaction{false},
		hasOpenDatabase{false},
		totalNumberOfRows{},
		vectorLength{0u},
		currentLoadId{std::nullopt}
	{
		loadColumnSpecData(specFilePath);
	};
	SqliteDataLoader(bool db, const std::string &dbPath, std::optional<size_t> maxLoadCount = std::nullopt) :
		IDataLoader{maxLoadCount},
		_columnSpec{},
		_tableNames{},
		_columnNames{},
		_weight{},
		_isBinary{},
		_tableName{},
		_dbPath{dbPath},
		_verbose{false}, 
		hasOpenTransaction{false},
		hasOpenDatabase{false},
		totalNumberOfRows{},
		vectorLength{0u},
		currentLoadId{std::nullopt} {}
		
	~SqliteDataLoader() override;
	void setTable(const std::string& name) noexcept;
	void setColumnSpec(const std::vector<ColumnSpec> columnSpec) noexcept override;
	const std::vector<ColumnSpec> getColumnSpec() noexcept override;
	size_t load() override;
	std::vector<RowData> getPreview(size_t count) override;
	bool open(const char *dbPath) override;
	bool open();
	std::vector<std::string> findAllColumns() override;
	std::vector<std::string> findAllTables();

	void addColumnName(std::string columnName) noexcept;
	void addTableName(std::string tableName) noexcept;

	const std::vector<float> getWeights() const noexcept override;
	float getWeight(size_t index) override;
	const std::vector<int> &getBinary() const noexcept override;
	const std::vector<int> &getContinuous() const noexcept override;
	std::string getName(size_t index) const noexcept;
	const std::vector<std::string> getNames() const noexcept override;
	size_t getDepth() const noexcept;
	bool isAtStartOfDataStream() const noexcept override 
	{ 
		return !currentLoadId.has_value();
	}

protected:
	std::vector<ColumnSpec> _columnSpec;
	std::vector<std::string> _tableNames;
	std::vector<std::string> _columnNames;
	std::vector<float> _weight;
	std::vector<int> _isBinary;
	std::vector<int> _isContinuous;
	std::string _tableName;
	std::string _dbPath;
	sqlite3 *db;
	bool _verbose;
	bool hasOpenTransaction, hasOpenDatabase;
	size_t totalNumberOfRows;
	size_t vectorLength;
	std::optional<int> currentLoadId;

	void close();
	int getMax(char *, std::string);
	int getRecord(char *[], int, std::string) const;
	int startTransaction();
	int endTransaction();
	size_t minId();
	size_t maxId();
	bool doesExist(size_t);
	int getElement(char *, size_t, std::string);
	void executeStringStatement(std::stringstream& statementStream, int(func)(void*, int, char**, char**));
	size_t numberOfRowsToLoad(size_t minId, size_t maxId);
	void loadColumnSpecData(const std::string &path);
	std::tuple<std::vector<RowData>, std::optional<int>, int> fetchData(std::optional<size_t> startId, std::optional<size_t> maxCount);
	std::tuple<std::vector<RowData>, std::optional<int>, int> fetchData2(std::optional<size_t> startId, std::optional<size_t> maxCount);

};