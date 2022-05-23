#include "SqliteDataLoader.hpp"

#include <iostream>
#include <fstream>
#include <cstring>
#include <stdexcept>

void SqliteDataLoader::loadColumnSpecData(const std::string &path)
{
	auto inFile = std::ifstream();
	inFile.open(path);

	std::string tempLine;
	while (std::getline(inFile, tempLine))
	{
		int currentColumn = 0;
		std::istringstream iline(tempLine);

		std::string token;
		std::string name{};
		double parsedValue = 1.0;
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
				parsedValue = strtod(token.c_str(), &ptr);

				// If not a number, set weight to 1
				if (*ptr)
				{
					parsedValue = 1.0;
					if (token.compare("binary") == 0)
					{
						// std::cout << name << " is treated as a binary variable.\n";
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
		parsedValue = 1;
		if (_verbose)
			std::cout << _columnNames.back() << "\tWeight:" << _weight.back() << "\tBinary:" << _isBinary.back() << "\n";
	}
	vectorLength = _columnNames.size();

	inFile.close();
}

int saveElement(void *buff, int argc, char **argv, char **azColName)
{
	char *ptr = (char *)buff;

	// If value string is not NULL
	if (argv[0])
	{
		strcpy(ptr, argv[0]);
		// std::cout << azColName[0] << ":\t" << ptr << "\n";
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

const std::vector<float> &SqliteDataLoader::getWeights() const noexcept
{
	return _weight;
}

std::string SqliteDataLoader::getName(size_t index) const noexcept
{
	assert(_columnNames.size() > index);
	return _columnNames[index];
}

const std::vector<std::string> &SqliteDataLoader::getNames() const noexcept
{
	return _columnNames;
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

int SqliteDataLoader::getElement(char *buff, int row, std::string icanName)
{
	char sql[200];
	char *err_msg;

	sprintf(sql, "SELECT %s FROM Ican WHERE Id = '%d';", icanName.c_str(), row);

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

size_t SqliteDataLoader::numberOfRows()
{
	char *err_msg;
	char sql[] = "SELECT COUNT(*) FROM Ican;";
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

int SqliteDataLoader::minId()
{
	char *err_msg;
	char sql[] = "SELECT MIN(Id) FROM Ican;";
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

int SqliteDataLoader::maxId()
{
	char *err_msg;
	char sql[] = "SELECT MAX(Id) FROM Ican;";
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

int SqliteDataLoader::doesExist(int row)
{
	char *err_msg;
	char sql[200];
	char buff[50];
	unsigned long value;

	sprintf(sql, "SELECT COUNT(*) FROM Ican WHERE Id = \'%d\';", row);

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

	return ((int)value);
}

int SqliteDataLoader::getMax(char *buff, std::string icanName)
{
	char sql[200];
	char *err_msg;

	sprintf(sql, "SELECT MAX(CAST(%s AS DECIMAL)) FROM Ican;", icanName.c_str());

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

size_t SqliteDataLoader::load(std::optional<size_t> maxCount)
{
	totalNumberOfRows = numberOfRows();

	auto numberOfRowsToLoad = maxCount.value_or(totalNumberOfRows) >= totalNumberOfRows ? totalNumberOfRows : maxCount.value();

	data = std::vector<RowData>();
	data.reserve(numberOfRowsToLoad);

	const auto depthVar = getDepth();
	int currentRow = 0;

	startTransaction();
	// Loop over all records in the entire database
	for (int r = minId(); r < maxId() + 1; ++r)
	{
		auto currentRowData = RowData{Eigen::VectorXf(depthVar), std::vector<int>(depthVar)};
		// Loop over all columns of record
		for (size_t i = 0; i < depthVar; i++)
		{
			char valueBuffer[20];

			if (doesExist(r))
				getElement(valueBuffer, r, _columnNames[i]);
			else
			{
				--i;
				++r;
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

		data.push_back(std::move(currentRowData));

		++currentRow;

		if ((((maxId() - minId()) > 100 && (r % (int)((maxId() - minId()) / 100)) == 0) || (maxId() - minId()) < 100) && _verbose)
		{
			if (_verbose)
				std::cout << "\rLoading database:" << 100 * r / (maxId() - minId()) << "%";
		}
	}
	if (_verbose)
		std::cout << "\rLoading database:100%\n";

	endTransaction();

	return currentRow;
}