#include "SqliteDataLoader.hpp"

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <algorithm>

// Deprecated
void SqliteDataLoader::loadColumnSpecData(const std::string &path)
{
	std::cerr << "Use of deprecated method SqliteDataLoader::loadColumnSpecData(const std::string &path)\n";

	auto inFile = std::ifstream();
	inFile.open(path);

	std::string tempLine;
	while (std::getline(inFile, tempLine))
	{
		int currentColumn = 0;
		std::istringstream iline(tempLine);

		std::string token;
		std::string name{};
		float parsedValue = 1.0f;
		bool isBinary = false;
		while (std::getline(iline, token, '\t'))
		{
			// Extract first column; database column name
			if (std::isalnum(token[0]) && currentColumn == 0)
				name = token;
			// Extract second column; database column weight
			else if (currentColumn > 0)
			{
				char *ptr;
				parsedValue = strtof(token.c_str(), &ptr);

				// If not a number, set weight to 1
				if (*ptr)
				{
					parsedValue = 1.0;
					if (token.compare("binary") == 0)
					{
						isBinary = true;
					}
				}
			}
			++currentColumn;
		}
		_columnNames.push_back(name);
		_weight.push_back(parsedValue);
		_isContinuous.push_back(!isBinary);
		_isBinary.push_back(isBinary);

		_columnSpec.emplace_back(_columnNames.back(), _weight.back(), _isBinary.back());

		parsedValue = 1;
		if (_verbose)
			std::cout << _columnNames.back() << "\tWeight:" << _weight.back() << "\tBinary:" << _isBinary.back() << "\n";
	}
	vectorLength = _columnNames.size();

	inFile.close();

	_tableName = "ican";
}

void SqliteDataLoader::setTable(const std::string& name) noexcept
{
	_tableName = name;
}

void SqliteDataLoader::setColumnSpec(const std::vector<ColumnSpec> columnSpec) noexcept
{
	_columnNames.clear();
	_weight.clear();
	_isBinary.clear();
	_isContinuous.clear();
	_columnSpec.clear();

	auto numberOfElements = columnSpec.size();

	_columnNames.reserve(numberOfElements);
	_weight.reserve(numberOfElements);
	_isBinary.reserve(numberOfElements);
	_isContinuous.reserve(numberOfElements);
	_columnSpec.reserve(numberOfElements);

	for (auto spec : columnSpec)
	{
		_columnNames.emplace_back(spec.name);
		_weight.emplace_back(spec.weight);
		_isBinary.emplace_back(spec.isBinary);
		_isContinuous.emplace_back(!spec.isBinary);

		_columnSpec.emplace_back(_columnNames.back(), _weight.back(), _isBinary.back());
	}

	vectorLength = _columnNames.size();
}

const std::vector<ColumnSpec> SqliteDataLoader::getColumnSpec() noexcept
{
	return _columnSpec;
}

int saveElement(void *buff, int argc, char **argv, char **azColName)
{
	char *ptr = (char *)buff;

	// If value string is not NULL
	if (argv[0])
	{
		strcpy(ptr, argv[0]);
	}
	else
	{
		strcpy(ptr, "NULL");
	}
	return 0;
}

SqliteDataLoader::~SqliteDataLoader()
{
	sqlite3_close(db);
}

bool SqliteDataLoader::open()
{
	return open(_dbPath.c_str());
}

bool SqliteDataLoader::open(const char *fileName)
{
	int rc = sqlite3_open_v2(fileName, &db, SQLITE_OPEN_READONLY, NULL);

	if (_verbose)
		std::cout << "Opening database: " << fileName << "\n";

	// If error
	if (rc != SQLITE_OK)
	{

		std::cout << "Cannot open database: " << sqlite3_errmsg(db) << "\n";
		sqlite3_close(db);

		return false;
	}

	hasOpenDatabase = true;

	return true;
}

void SqliteDataLoader::close()
{
	sqlite3_close(db);

	hasOpenDatabase = false;
	hasOpenTransaction = false;
}

int SqliteDataLoader::getRecord(char *buff[], int row, std::string icanName) const
{
	return 0;
}

const std::vector<int> &SqliteDataLoader::getBinary() const noexcept
{
	return _isBinary;
}

const std::vector<int> &SqliteDataLoader::getContinuous() const noexcept
{
	return _isContinuous;
}

const std::vector<float> SqliteDataLoader::getWeights() const noexcept
{
	std::vector<float> res{};
	std::transform(_columnSpec.begin(), _columnSpec.end(), std::back_inserter(res), [](auto spec) { return spec.weight; });
	return res;
}

float SqliteDataLoader::getWeight(size_t index)
{
	return _columnSpec.at(index).weight;
}

std::string SqliteDataLoader::getName(size_t index) const noexcept
{
	return _columnSpec.at(index).name;
}

const std::vector<std::string> SqliteDataLoader::getNames() const noexcept
{
	std::vector<std::string> res{};
	std::transform(_columnSpec.begin(), _columnSpec.end(), std::back_inserter(res), [](auto spec) { return spec.name; });
	return res;
}

size_t SqliteDataLoader::getDepth() const noexcept
{
	return vectorLength;
}

int SqliteDataLoader::startTransaction()
{
	char *err_msg;
	int retValue{};
	retValue = sqlite3_exec(db, "BEGIN;", NULL, NULL, &err_msg);
	if (err_msg)
	{
		fprintf(stderr, "Failed to BEGIN transaction: Error message:%s\n", err_msg);
		sqlite3_free(err_msg);
		return (-1);
	}

	hasOpenTransaction = true;

	return retValue;
}

int SqliteDataLoader::endTransaction()
{
	char *err_msg;
	int retValue{};
	retValue = sqlite3_exec(db, "END;", NULL, NULL, &err_msg);
	if (err_msg)
	{
		fprintf(stderr, "Failed to END transaction: Error message:%s\n", err_msg);
		sqlite3_free(err_msg);
		return (-1);
	}

	hasOpenTransaction = false;

	return retValue;
}

int saveColumnName(void *buff, int argc, char **argv, char **azColName)
{
	auto instance = static_cast<SqliteDataLoader *>(buff);

	if (argc == 1)
	{
		instance->addColumnName(argv[0]);
	}

	return 0;
}

int saveTableName(void *buff, int argc, char **argv, char **azColName)
{
	auto instance = static_cast<SqliteDataLoader *>(buff);

	if (argc == 1)
	{
		instance->addTableName(argv[0]);
	}

	return 0;
}

void SqliteDataLoader::addColumnName(std::string columnName) noexcept
{
	_columnNames.push_back(columnName);
}

void SqliteDataLoader::addTableName(std::string tableName) noexcept
{
	_tableNames.push_back(tableName);
}

std::vector<std::string> SqliteDataLoader::findAllColumns()
{
	if (!hasOpenDatabase)
	{
		std::cerr << "No open database - cannot parse columns\n";
		return std::vector<std::string>{};
	}

	_columnNames.clear();

	std::stringstream ss;

	ss << "select name from pragma_table_info('" << _tableName << "') as tblInfo;";

	executeStringStatement(ss, saveColumnName);

	return _columnNames;
}

std::vector<std::string> SqliteDataLoader::findAllTables()
{
	if (!hasOpenDatabase)
	{
		std::cerr << "No open database - cannot parse columns\n";
		return std::vector<std::string>{};
	}

	_tableNames.clear();

	std::stringstream ss;

	ss << "SELECT name FROM sqlite_master WHERE type = 'table' AND name NOT LIKE 'sqlite_%';";

	executeStringStatement(ss, saveTableName);

	return _tableNames;
}

void SqliteDataLoader::executeStringStatement(std::stringstream& statementStream, int(func)(void*, int, char**, char**))
{
	const std::string &tmp = statementStream.str();
	const char *sql = tmp.c_str();

	char *err_msg;
	sqlite3_exec(db, sql, func, static_cast<void *>(this), &err_msg);

	if (err_msg)
	{
		sqlite3_free(err_msg);
		throw std::runtime_error(std::string("Failed to execute statement: ") + sql + "\nError message:" + err_msg);
	}
}

int SqliteDataLoader::getElement(char *buff, size_t row, std::string icanName)
{
	char sql[200];
	char *err_msg;

	sprintf(sql, "SELECT %s FROM %s WHERE Id = '%ld';", icanName.c_str(), _tableName.c_str(), row);

	// Exekvera SQL statement
	sqlite3_exec(db, sql, saveElement, buff, &err_msg);

	if (err_msg)
	{
		fprintf(stderr, "Failed to execute statement: %s\nError message:%s\n", sql, err_msg);
		sqlite3_free(err_msg);
		return (-1);
	}

	return (0);
}

size_t SqliteDataLoader::numberOfRowsToLoad(size_t minId, size_t maxId)
{
	char *err_msg;
	char sql[200];
	sprintf(sql, "SELECT COUNT(*) FROM %s WHERE Id >= \'%ld\' AND Id <= \'%ld\';", _tableName.c_str(), minId, maxId);

	char buff[50];
	unsigned long value;

	// Exekvera SQL statement
	sqlite3_exec(db, sql, saveElement, buff, &err_msg);

	if (err_msg)
	{
		fprintf(stderr, "Failed to execute statement: %s\nError message:%s\n", sql, err_msg);
		sqlite3_free(err_msg);
		return (-1);
	}

	value = strtoul(buff, NULL, 10);

	return ((int)value);
}

size_t SqliteDataLoader::minId()
{
	char *err_msg;
	char sql[200];
	sprintf(sql, "SELECT MIN(Id) FROM %s;", _tableName.c_str());
	char buff[50];
	size_t value;

	// Exekvera SQL statement
	sqlite3_exec(db, sql, saveElement, buff, &err_msg);

	if (err_msg)
	{
		fprintf(stderr, "Failed to execute statement: %s\nError message:%s\n", sql, err_msg);
		sqlite3_free(err_msg);
		return (-1);
	}

	value = strtoul(buff, NULL, 10);

	return value;
}

size_t SqliteDataLoader::maxId()
{
	char *err_msg;
	char sql[200];
	sprintf(sql, "SELECT MAX(Id) FROM %s;", _tableName.c_str());
	char buff[50];
	size_t value;

	// Exekvera SQL statement
	sqlite3_exec(db, sql, saveElement, buff, &err_msg);

	if (err_msg)
	{
		fprintf(stderr, "Failed to execute statement: %s\nError message:%s\n", sql, err_msg);
		sqlite3_free(err_msg);
		return (-1);
	}

	value = strtoul(buff, NULL, 10);

	return value;
}

bool SqliteDataLoader::doesExist(size_t row)
{
	char *err_msg;
	char sql[200];
	char buff[50];
	size_t value;

	sprintf(sql, "SELECT COUNT(*) FROM %s WHERE Id = \'%ld\';", _tableName.c_str(), row);

	// Exekvera SQL statement
	sqlite3_exec(db, sql, saveElement, buff, &err_msg);

	if (err_msg)
	{
		fprintf(stderr, "Failed to execute statement: %s\nError message:%s\n", sql, err_msg);
		sqlite3_free(err_msg);
		return (-1);
	}

	value = strtoul(buff, NULL, 10);

	// if( !value ) printf("\nDoes not exist\n");

	return value == 1;
}

int SqliteDataLoader::getMax(char *buff, std::string icanName)
{
	char sql[200];
	char *err_msg;

	sprintf(sql, "SELECT MAX(CAST(%s AS DECIMAL)) FROM %s;", icanName.c_str(), _tableName.c_str());

	// Exekvera SQL statement
	sqlite3_exec(db, sql, saveElement, buff, &err_msg);

	if (err_msg)
	{
		fprintf(stderr, "Failed to execute statement: %s\nError message:%s\n", sql, err_msg);
		sqlite3_free(err_msg);
		return (-1);
	}

	return (0);
}

size_t SqliteDataLoader::load()
{
	auto [loadedData, currentId, maxId] = fetchData2(currentLoadId, m_maxLoadCount);

	data = loadedData;
	currentLoadId = currentId;

	if (_verbose)
		std::cout << "\rLoading database:100%\n";

	if (currentLoadId > maxId)
		currentLoadId.reset();

	return data.size();
}

std::tuple<std::vector<RowData>, std::optional<int>, int> SqliteDataLoader::fetchData2(std::optional<size_t> currentId, std::optional<size_t> maxCount)
{
	if(_tableName.empty())
	{
		std::cerr << "Cannot read database - No table name specified\n";
		return std::make_tuple(std::vector<RowData>{}, std::nullopt, 0);
	}

	auto depth = getDepth();

	startTransaction();

	auto maxIdValue = maxId();

	currentId = currentId.has_value() ? currentId.value() : minId();
	auto maxLoadId = maxCount.has_value() ? (currentId.value() + static_cast<long>(maxCount.value())) > maxIdValue ? maxIdValue : currentId.value() + maxCount.value() - 1
										  : maxIdValue;

	totalNumberOfRows = numberOfRowsToLoad(currentId.value(), maxLoadId);

	auto newData = std::vector<RowData>();
	newData.reserve(totalNumberOfRows);

	std::stringstream ss;
	ss << "SELECT ";
	for(auto column : getNames())
		ss << column << ',';
	ss << "Id FROM " << _tableName << " WHERE Id>=" << currentId.value() << " AND Id <=" << maxLoadId << ';';
	auto query = ss.str();

	sqlite3_stmt *res;
	const char *tail;

	if(sqlite3_prepare_v2(db, query.c_str(), static_cast<int>(query.length()), &res, &tail) != SQLITE_OK)
	{
		std::cerr << "Cannot fetch data: " << sqlite3_errmsg(db) << '\n';
		return std::make_tuple(std::vector<RowData>{}, std::nullopt, 0);
	}

	while(sqlite3_step(res) == SQLITE_ROW)
	{
		auto currentRowData = RowData{Eigen::VectorXf(depth), std::vector<int>(depth)};

		for(int columnIndex{}; columnIndex < static_cast<int>(depth); ++columnIndex)
		{
			/* TODO: Check for value type somehow */
			currentRowData.values[columnIndex] = static_cast<float>(sqlite3_column_double(res, columnIndex));
			currentRowData.valid[columnIndex] = 1;
		}

		currentId = sqlite3_column_int64(res, static_cast<int>(depth));

		newData.push_back(std::move(currentRowData));
	}

	endTransaction();

	/* Free res */
	sqlite3_finalize(res);

	std::cout << "Fetched whole batch\n";

	return std::make_tuple(newData, currentId.emplace(currentId.value() + 1), maxIdValue);
}

std::tuple<std::vector<RowData>, std::optional<int>, int> SqliteDataLoader::fetchData(std::optional<size_t> currentId, std::optional<size_t> maxCount)
{
	const auto depth = getDepth();

	if(_tableName.empty())
	{
		std::cerr << "Cannot read database - No table name specified\n";
		return std::make_tuple(std::vector<RowData>{}, std::nullopt, 0);
	}

	startTransaction();

	auto maxIdValue = maxId();

	currentId = currentId.has_value() ? currentId.value() : minId();
	auto maxLoadId = maxCount.has_value() ? (currentId.value() + static_cast<long>(maxCount.value())) > maxIdValue ? maxIdValue : currentId.value() + maxCount.value() - 1
										  : maxIdValue;

	totalNumberOfRows = numberOfRowsToLoad(currentId.value(), maxLoadId);

	auto newData = std::vector<RowData>();
	newData.reserve(totalNumberOfRows);


	while (currentId <= maxLoadId)
	{
		if(currentId.value() % 100 == 0)
			std::cout << "\rFetching data with id = [" << currentId.value() << ", " << maxLoadId << "]";

		auto currentRowData = RowData{Eigen::VectorXf(depth), std::vector<int>(depth)};
		// Loop over all columns of record
		for (size_t i = 0; i < depth; i++)
		{
			char valueBuffer[20];

			if (doesExist(currentId.value()))
				getElement(valueBuffer, currentId.value(), _columnNames[i]);
			else
			{
				--i;
				continue;
			}

			if (valueBuffer[0] == '\0' || strcmp(valueBuffer, "NULL") == 0)
			{
				currentRowData.values[i] = 0.0f;
				currentRowData.valid[i] = 0;
			}
			else
			{
				char *ptr;
				currentRowData.values[i] = strtof(valueBuffer, &ptr);
				currentRowData.valid[i] = 1;
			}
		}

		newData.push_back(std::move(currentRowData));

		currentId.emplace(currentId.value() + 1);
	}

	endTransaction();

	std::cout << "Fetched whole batch\n";

	return std::make_tuple(newData, currentId, maxIdValue);
}

std::vector<RowData> SqliteDataLoader::getPreview(size_t count)
{
	auto [loadedData, currentId, maxId] = fetchData2(std::nullopt, 100);
	return loadedData;
}

#include "doctest/doctest.h"

TEST_CASE("testing file stuff")
{
	auto sut = SqliteDataLoader(true, "../tests/performance/data/testDb.sq3");
	sut.open();

	SUBCASE("print all tables in database")
	{
		auto tables = sut.findAllTables();

		for (auto table : tables)
		{
			std::cout << table << '\n';
		}
	}

	SUBCASE("print all columns of table")
	{
		std::cout << "Choosing table 'ican'\n";
		sut.setTable("ican");

		auto columns = sut.findAllColumns();

		for (auto col : columns)
		{
			std::cout << col << '\n';
		}
	}
	SUBCASE("compareLoadingSpecfileWithDynamicSpecs")
	{
		/* First remove, then write new columnSpecFile */
		std::remove("tempColumnSpec.txt");
		std::ofstream ofs;
		ofs.open("tempColumnSpec.txt", std::ios::out);

		ofs << 
"A	1\n\
B	1\n\
C	1\n\
D	1\n\
E	1	binary\n\
F	1\n\
G	1\n\
H	0.42\n\
I	1";

		ofs.close();

		auto fileSut = SqliteDataLoader("tempColumnSpec.txt");
		fileSut.open();
		auto specFileColumnNames = fileSut.getNames();

		std::string colA{"A"};
		std::string colB{"B"};
		std::string colC{"C"};
		std::string colD{"D"};
		std::string colE{"E"};
		std::string colF{"F"};
		std::string colG{"G"};
		std::string colH{"H"};
		std::string colI{"I"};
		sut.setColumnSpec(std::vector<ColumnSpec>{
			ColumnSpec(colA, 1, false),
			ColumnSpec(colB, 1, false),
			ColumnSpec(colC, 1, false),
			ColumnSpec(colD, 1, false),
			ColumnSpec(colE, 1, true),
			ColumnSpec(colF, 1, false),
			ColumnSpec(colG, 1, false),
			ColumnSpec(colH, 0.42f, false),
			ColumnSpec(colI, 1, false),
		});

		/* Check that column names are equal */
		auto dynamicColumnNames = sut.getNames();

		CHECK(specFileColumnNames.size() == dynamicColumnNames.size());

		for (size_t index{0}; index < specFileColumnNames.size(); ++index)
		{
			REQUIRE(specFileColumnNames.at(index) == dynamicColumnNames.at(index));
		}

		/* Check that column binary is equal */
		auto specFileColumnBinary = fileSut.getBinary();
		auto dynamicColumnBinary = sut.getBinary();

		CHECK(specFileColumnBinary.size() == dynamicColumnBinary.size());

		for (size_t index{0}; index < specFileColumnBinary.size(); ++index)
		{
			REQUIRE(specFileColumnBinary.at(index) == dynamicColumnBinary.at(index));
		}

		/* Check that column weights are equal */
		auto specFileColumnWeights = fileSut.getWeights();
		auto dynamicColumnWeights = sut.getWeights();

		CHECK(specFileColumnWeights.size() == dynamicColumnWeights.size());

		for (size_t index{0}; index < specFileColumnWeights.size(); ++index)
		{
			REQUIRE(specFileColumnWeights.at(index) == dynamicColumnWeights.at(index));
		}
	}
}
