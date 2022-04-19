#include "DataSet.hpp"

#include <iostream>
#include <fstream>
#include <ctype.h>

DataSet::DataSet(const char *fileName, bool verbose)
{
	_verbose = verbose;

	if (_verbose)
		std::cout << "Opening: " << fileName << "\n";

	std::ifstream inFile;
	inFile.open(fileName);

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
		variableNames.push_back(name);
		binary.push_back(isBinary);
		continuous.push_back(!isBinary);
		weight.push_back(parsedValue);
		parsedValue = 1;
		if(_verbose) std::cout << variableNames.back() << "\tWeight:" << weight.back() << "\tBinary:" << binary.back() << "\n";
	}

	inFile.close();

	depth = variableNames.size();

	if (_verbose)
		std::cout << "Dataset depth: " << depth << "\n";
}

Eigen::VectorXf DataSet::getData(size_t index) const
{
	assert(n > index);
	return data[index];
}

const Eigen::VectorXi DataSet::getValidity(size_t index) const
{
	assert(n > index);
	return Eigen::Map<const Eigen::VectorXi>(valid[index].data(), valid[index].size());
}

const Eigen::ArrayXi DataSet::getBinary() const
{
	return Eigen::Map<const Eigen::ArrayXi>(binary.data(), binary.size());
}

const Eigen::ArrayXi DataSet::getContinuous() const
{
	return Eigen::Map<const Eigen::ArrayXi>(continuous.data(), continuous.size());
}

size_t &DataSet::getLastBMU(size_t index)
{
	assert(n > index);
	return lastBMU[index];
}

const std::vector<std::string> &DataSet::getNames() const noexcept
{
	return variableNames;
}

std::string DataSet::getName(size_t index) const
{
	assert(depth > index);
	return variableNames[index];
}

const Eigen::VectorXf DataSet::getWeights() const
{
	return Eigen::Map<const Eigen::VectorXf>(weight.data(), weight.size());
}

size_t DataSet::size() const
{
	return n;
}

void DataSet::addVector(Eigen::VectorXf v)
{
	if (static_cast<size_t>(v.rows()) == depth)
	{
		data.push_back(v);
		n += 1;
	}
	else
		std::cout << "Added vector size does not correspond to data set depth!\n";
}

void DataSet::addSet(const DataSet &d)
{
	if (d.vectorLength() == depth)
	{
		data.reserve(n + d.size());
		for (size_t i = 0; i < d.size(); ++i)
		{
			data.push_back(d.getData(i));
			n += 1;
		}
	}
	else
		std::cout << "Data set depths do not correspond!\n";
}

void DataSet::loadDataBase(DataBase *db)
{
	const int numberOfRows = db->rows();
	
	// Create template vector for whole data set
	const Eigen::VectorXf initV(depth);

	// Initialize new element of outer list
	valid.resize(numberOfRows, std::vector<int>(depth, 0));

	data.resize(numberOfRows, initV);

	lastBMU.resize(numberOfRows, 0);

	// Shuffle indices for data and valid maps
	index.resize(numberOfRows);
	for (size_t k = 0; k < index.size(); ++k)
		index[k] = k;
	shuffle();

	n = numberOfRows;

	db->startTransaction();

	if (_verbose)
		std::cout << "Min:" << db->minId() << "\tMax:" << db->maxId() << "\n";

	// Loop over all records in the entire database
	for (int currentRow = 0, r = db->minId(); r < db->maxId() + 1; ++r)
	{
		// Loop over all columns of record
		for (size_t i = 0; i < depth; i++)
		{
			char valueBuffer[20];

			if (db->doesExist(r))
				db->getElement(valueBuffer, r, variableNames[i]);
			else
			{
				--i;
				++r;
				continue;
			}
			// Save dataset with shuffled indices in data and valid
			if (valueBuffer[0] == '\0' || strcmp(valueBuffer, "NULL") == 0)
			{
				data[index[currentRow]](i) = 0;
				valid[index[currentRow]][i] = 0;
			}
			else
			{
				// Save dataset with shuffled indices in data and valid
				char *ptr;
				data[index[currentRow]](i) = strtod(valueBuffer, &ptr);
				valid[index[currentRow]][i] = 1;
			}
		}
		++currentRow;

		if ((((db->maxId() - db->minId()) > 100 
			&& (r % (int)((db->maxId() - db->minId()) / 100)) == 0) 
			|| (db->maxId() - db->minId()) < 100) 
			&& _verbose)
		{
			if (_verbose)
				std::cout << "\rLoading database:" << 100 * r / (db->maxId() - db->minId()) << "%";
		}
	}
	if (_verbose)
		std::cout << "\rLoading database:100%\n";

	db->endTransaction();
}

// This function should load a csv file into a DataSet.
// yet to make functional
void DataSet::loadTextFile(const char *fileName)
{
	FILE *fp;
	std::vector<char> lineString;
	Eigen::VectorXf tempData(depth);
	char ch,
		*ptr,
		*nextColumn;
	int numberOfNewlines = 0,
		numberOfColumns = 0,
		maxNumberOfColumns = 0,
		minNumberOfColumns = 1000000,
		i = 0,
		icanIndex = 0;
	size_t vectorElement = 0;
	double tempValue;

	if ((fp = fopen(fileName, "r")) == NULL)
	{
		std::cout << "Couldn't open file " << fileName << ". Quitting...\n";
		exit(EXIT_FAILURE);
	}

	// Load complete file
	while ((ch = fgetc(fp)) != EOF)
		lineString.push_back(ch);

	// Count number of lines
	for (unsigned int i = 0; i < lineString.size(); i++)
	{
		if (lineString[i] == '\n')
		{
			numberOfNewlines += 1;
			maxNumberOfColumns = numberOfColumns > maxNumberOfColumns ? numberOfColumns : maxNumberOfColumns;
			minNumberOfColumns = numberOfColumns < minNumberOfColumns ? numberOfColumns : minNumberOfColumns;
			numberOfColumns = 0;
		}
		else if (lineString[i] == '\t')
			numberOfColumns++;
	}

	if (_verbose)
		std::cout << fileName << ": File size: " << lineString.size() << "\nNumber of lines: " << numberOfNewlines << "\nBiggest depth: " << maxNumberOfColumns << "\nSmallest depth: " << minNumberOfColumns << "\n";

	nextColumn = &lineString[0];

	// Parse loaded file
	while (i < numberOfNewlines)
	{
		// Skip first line if commented with leading '#'
		if (lineString[0] == '#' && i == 0)
		{
			int j = 1;
			i++;

			do
			{
				j++;
			} while (lineString[j] != '\n');

			nextColumn = &lineString[j];
		}
		tempValue = strtod(nextColumn, &ptr);

		nextColumn = ptr;

		// std::cout << "vector element: " << vectorElement << "\ticanIndex: " << icanIndex << "\tptr: " << (char)*(ptr+1) << " value: " << tempValue  << "\n";

		/* HÄR ÄR JAG OCH KRÅNGLAR JUST NU!!!!
		 * ----------------------------------*/
		if ((std::isdigit(*(ptr + 1)) || *(ptr + 1) == '-') && vectorElement < depth /*&& filter[icanIndex]*/)
		{
			tempData[vectorElement] = tempValue;
			vectorElement++;
			icanIndex++;
		}
		else if ((std::isdigit(*(ptr + 1)) || *(ptr + 1) == '-') && vectorElement < depth)
		{
			icanIndex++;
		}
		else if (*(ptr + 1) == '\n')
		{
			this->addVector(tempData);
			// std::cout << "vector size: " << vectorElement << "\n";
			vectorElement = 0;
			icanIndex = 0;
			i++;
		}
		else
		{
			// std::cout << "Omitting ican index " << icanIndex << "\n";
			icanIndex++;
			// std::cout << "File " << fileName << " and data set depths do not correspond. Quitting...\n";
			// exit(EXIT_FAILURE);
		}
	}
}

void DataSet::display() const
{
	std::cout << "Number of samples: " << n << "\nVector length: " << depth << "\n";
}

size_t DataSet::vectorLength() const
{
	return depth;
}

// Shuffle data set record indices. I.e. shuffle references, but keep order in memory
void DataSet::shuffle()
{
	std::random_shuffle(index.begin(), index.end());
}