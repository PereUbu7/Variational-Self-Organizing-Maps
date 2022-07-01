#include "DataSet.hpp"
#include "SqliteDataLoader.hpp"

#include <iostream>
#include <fstream>
#include <ctype.h>

const std::vector<DataSet::DataRow> DataSet::getAll() const
{
	return allData;
}

std::vector<DataSet::DataRow> DataSet::getAll()
{
	return allData;
}

std::vector<Eigen::VectorXf> DataSet::getPreviewData(size_t count) const
{
	auto data = _loader.getPreview(count);

	auto previewData = std::vector<Eigen::VectorXf>();
	previewData.reserve(data.size());

	for(auto item : data)
	{
		previewData.emplace_back(
			item.values
		);
	}
	return previewData;
}

Eigen::VectorXf DataSet::getData(size_t index) const
{
	if(_loader.data.size() > index) 
		return data[index];
	else 
		return Eigen::VectorXf::Zero(_loader.getDepth());
}

const Eigen::VectorXi DataSet::getValidity(size_t index) const
{
	if(n > index)
		return Eigen::Map<const Eigen::VectorXi>(valid[index].data(), valid[index].size());
	else
		return Eigen::VectorXi::Zero(_loader.getDepth());
}

const Eigen::ArrayXi DataSet::getBinary() const
{
	return Eigen::Map<const Eigen::ArrayXi>(_loader.getBinary().data(), _loader.getBinary().size());
}

const Eigen::ArrayXi DataSet::getContinuous() const
{
	return Eigen::Map<const Eigen::ArrayXi>(_loader.getContinuous().data(), _loader.getContinuous().size());
}

const std::vector<size_t> &DataSet::getLastBMU() const noexcept
{
	return lastBMU;
}

size_t &DataSet::getLastBMU(size_t index)
{
	assert(n > index);
	return lastBMU[index];
}

const std::vector<std::string> &DataSet::getNames() const noexcept
{
	return _loader.getNames();
}

std::string DataSet::getName(size_t index) const
{
	return _loader.getName(index);
}

const Eigen::VectorXf DataSet::getWeights() const
{
	return Eigen::Map<const Eigen::VectorXf>(_loader.getWeights().data(), _loader.getWeights().size());
}

float &DataSet::getWeight(size_t index)
{
	return _loader.getWeight(index);
}

size_t DataSet::size() const
{
	return n;
}

// Used by csv loader - Should become its own data loader
void DataSet::addVector(Eigen::VectorXf v)
{
	if (static_cast<size_t>(v.rows()) == _loader.getDepth())
	{
		data.push_back(v);
		n += 1;
	}
	else
		std::cout << "Added vector size does not correspond to data set depth!\n";
}

void DataSet::resetStreamLoadPosition() noexcept
{
	loadedNumberOfChunks = 0;
}

bool DataSet::hasReadWholeDataStream() const noexcept
{
	return loadedNumberOfChunks > 0 && _loader.isAtStartOfDataStream();
}

void DataSet::loadNextDataFromStream()
{
	if(_loader.isAtStartOfDataStream()) loadedNumberOfChunks = 0;

	const size_t numberOfRows = _loader.load();

	n = numberOfRows;
	
	// Create template vector for whole data set
	const Eigen::VectorXf initV(_loader.getDepth());

	// Initialize new element of outer list
	valid = std::vector<std::vector<int>>();
	valid.reserve(numberOfRows);

	data = std::vector<Eigen::VectorXf>();
	data.reserve(numberOfRows);

	lastBMU = std::vector<size_t>();
	lastBMU.resize(numberOfRows, 0);

	allData = std::vector<DataSet::DataRow>();
	allData.reserve(numberOfRows);

	// Shuffle indices for data and valid maps
	index.resize(numberOfRows);
	for (size_t k = 0; k < index.size(); ++k)
		index[k] = k;
	shuffle();
	
	for(size_t currentIndex{0}; auto &[record, validity] : _loader.data)
	{
		data.push_back(record);
		valid.push_back(validity);
		allData.emplace_back(
			&(data.back()), 
			&(valid.back()), 
			&(lastBMU[currentIndex++])
			);
	}
	++loadedNumberOfChunks;
	std::cout << "Loaded " << loadedNumberOfChunks << " number of chunks\n";
}

// This function should load a csv file into a DataSet.
// yet to make functional
// TODO: Make into DataLoader
void DataSet::loadTextFile(const char *fileName)
{
	FILE *fp;
	std::vector<char> lineString;
	Eigen::VectorXf tempData(_loader.getDepth());
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
		if ((std::isdigit(*(ptr + 1)) || *(ptr + 1) == '-') && vectorElement < _loader.getDepth() /*&& filter[icanIndex]*/)
		{
			tempData[vectorElement] = tempValue;
			vectorElement++;
			icanIndex++;
		}
		else if ((std::isdigit(*(ptr + 1)) || *(ptr + 1) == '-') && vectorElement < _loader.getDepth())
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
	std::cout << "Number of samples: " << _loader.data.size() << "\nVector length: " << _loader.getDepth() << "\n";
}

size_t DataSet::vectorLength() const
{
	return _loader.getDepth();
}

// Shuffle data set record indices. I.e. shuffle references, but keep order in memory
void DataSet::shuffle()
{
	std::random_shuffle(index.begin(), index.end());
}