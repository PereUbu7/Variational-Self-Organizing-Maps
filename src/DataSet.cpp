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

const std::vector<std::string> DataSet::getNames() const noexcept
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

float DataSet::getWeight(size_t index)
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