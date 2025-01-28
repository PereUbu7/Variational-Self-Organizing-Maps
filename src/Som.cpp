#include "SOM.hpp"

#include <iostream>
#include <fstream>
#include <random>
#include <execution>
#include <ranges>
#include <array>

// Initialize an empty SOM
void Som::Construct(size_t inWidth, size_t inHeight, size_t inDepth, std::vector<std::string> names)
{
	// Create template vector for whole map
	Eigen::VectorXf initV(inDepth);

	// Set map to correct size with every Vector equal to initV
	map.resize(inWidth * inHeight, initV);
	sigmaMap.resize(inWidth * inHeight, initV);
	weightMap.resize(inWidth * inHeight);
	SMap.resize(inWidth * inHeight, initV);

	// Set BMU hits map to zero
	bmuHits.resize(inWidth * inHeight, 0u);

	// Set U-matrix to zero
	uMatrix.resize(inWidth * inHeight, 0);

	for (size_t i = 0; i < inWidth * inHeight; ++i)
	{
		// Initialize each element in each neuron to 0
		for (int n = 0; n < sigmaMap[i].rows(); ++n)
		{
			map[i](n) = 0.0f;
			sigmaMap[i](n) = 0.0f;
			SMap[i](n) = 0.0f;
		}
		weightMap[i] = 0.0f;
		bmuHits[i] = 0u;
		uMatrix[i] = 0.0;
	}

	// Save size of map to object
	width = inWidth;
	height = inHeight;
	depth = inDepth;

	transform.names = names;
}

// Initialize SOM from trained data
Som::Som(const char *SOMFileName)
{
	// Create template vector for whole map
	Eigen::VectorXf initV(0);
	initV = getSizeFromFile(SOMFileName);

	depth = initV.size();

	// Set map to correct size with every Vector equal to size returned by getSizeFromFile()
	map.resize(width * height, initV);
	sigmaMap.resize(width * height, initV);
	weightMap.resize(width * height);
	SMap.resize(width * height, initV);

	for (size_t i = 0; i < width * height; ++i)
	{
		// Initialize each element in each neuron to 0
		for (int n = 0; n < sigmaMap[i].size(); ++n)
		{
			sigmaMap[i](n) = 0;
			SMap[i](n) = 0;
		}
		weightMap[i] = 0;
	}

	// Set BMU hits map to zero
	bmuHits.resize(width * height, 0uz);

	// Set U-matrix to zero
	uMatrix.resize(width * height, 0);

	this->load(SOMFileName);
}

void Som::display() const
{
	std::cout << "Map size: " << map.size() << "\n";
	std::cout << "Map width: " << width << "\n";
	std::cout << "Map height: " << height << "\n";
	std::cout << "M size: " << map[0].rows() << "\n\n";

	// Print map
	/*for(unsigned int i = 0; i < map.size(); i++)
	{
		std::cout << "M" << i << "\t";
		for(int j = 0; j < map[0].rows(); j++)
			std::cout << map[i][j] << "\t";

		std::cout << "\n";
	}
	*/
	// Print Bmu hits
	for (size_t i = 0; i < this->getHeight(); ++i)
	{
		for (size_t j = 0; j < this->getWidth(); ++j)
			std::cout << bmuHits[i * width + j] << "\t";

		std::cout << "\n";
	}
}

// Returns the Euclidian squared distance between neuron at position pos and vector v.
// Considers only dimensions where valid is true.
// All dimensions are weighed individually by weights vector
double Som::euclidianWeightedDist(
	const SomIndex &pos, const Eigen::VectorXf &v,
	const Eigen::VectorXf &valid, const Eigen::VectorXf &weights) const
{
	auto intPos = pos.getY() * width + pos.getX();

	return euclidianWeightedDist(intPos, v, valid, weights);
}

double Som::euclidianWeightedDist(
	const size_t &pos, const Eigen::VectorXf &v,
	const Eigen::VectorXf &valid, const Eigen::VectorXf &weights) const
{
 	Eigen::ArrayXf sM(sigmaMap[pos].size());
	Eigen::VectorXf validEigen;

	sM = (sigmaMap[pos].array() < 0.00001f).select(0.00001f, sigmaMap[pos]);

	// validEigen is an Eigen vector with weight if valid and 0 if not
	validEigen = valid.array() * weights.array();

	auto comparer = transform.Comparer(v, this->getNeuron(pos), sM, validEigen);

	// (pos-v)*(pos-v)*weight*valid/sigma2
	// return (((comparer).array() / sM).matrix().dot(((comparer).array() * validEigen.array() / sM).matrix()));
	return (comparer).matrix().dot(comparer.matrix());
}

double Som::euclidianWeightedDistRaw(
	const size_t &pos, const Eigen::VectorXf &v,
	const Eigen::VectorXf &valid, const Eigen::VectorXf &weights) const
{
 	Eigen::ArrayXf sM(sigmaMap[pos].size());
	Eigen::VectorXf validEigen;

	sM = (sigmaMap[pos].array() < 0.00001f).select(0.00001f, sigmaMap[pos]);

	// validEigen is an Eigen vector with weight if valid and 0 if not
	validEigen = valid.array() * weights.array();

	// (pos-v)*(pos-v)*weight*valid/sigma2
	return (((this->getNeuron(pos) - v).array() / sM).matrix().dot(((this->getNeuron(pos) - v).array() * validEigen.array() / sM).matrix()));
}

UMatrix Som::getUMatrix() const noexcept
{
	return UMatrix{uMatrix, width, height};
}

Eigen::VectorXf Som::getWeigthMap() const noexcept
{
	return weightMap;
}

std::vector<size_t> Som::getBmuHits() const noexcept
{
	return bmuHits;
}

size_t Som::getHeight() const noexcept
{
	return height;
}

size_t Som::getWidth() const noexcept
{
	return width;
}

size_t Som::getDepth() const noexcept
{
	return depth;
}

size_t Som::getIndex(SomIndex i) const noexcept
{
	return i.getY() * width + i.getX();
}

Eigen::VectorXf Som::getNeuron(SomIndex i) const noexcept
{
	return map[i.getY() * width + i.getX()];
}

Eigen::VectorXf Som::getNeuron(size_t i) const noexcept
{
	return map[i];
}

Eigen::VectorXf Som::getSigmaNeuron(SomIndex i) const noexcept
{
	return sigmaMap[i.getY() * width + i.getX()];
}

Eigen::VectorXf Som::getSigmaNeuron(size_t i) const noexcept
{
	return sigmaMap[i];
}

std::vector<std::string> Som::getNeuronStrings(SomIndex index) const noexcept
{
	return transform.Displayer(map[index.getY() * width + index.getX()]);
}

std::vector<std::string> Som::getSigmaNeuronStrings(SomIndex index) const noexcept
{
	return transform.Displayer(sigmaMap[index.getY() * width + index.getX()]);
}

float Som::getMaxValueOfFeature(size_t modelVectorIndex) const
{
	assert(modelVectorIndex < static_cast<size_t>(map.at(0).size()));

	return (*std::max_element(map.begin(), map.end(),
							  [modelVectorIndex](const auto &a, const auto &b)
							  {
								  return a[modelVectorIndex] < b[modelVectorIndex];
							  }))[modelVectorIndex];
}

float Som::getMinValueOfFeature(size_t modelVectorIndex) const
{
	assert(modelVectorIndex < static_cast<size_t>(map.at(0).size()));

	return (*std::min_element(map.begin(), map.end(),
							  [modelVectorIndex](const auto &a, const auto &b)
							  {
								  return a[modelVectorIndex] < b[modelVectorIndex];
							  }))[modelVectorIndex];
}

float Som::getMaxSigmaOfFeature(size_t modelVectorIndex) const
{
	assert(modelVectorIndex < static_cast<size_t>(sigmaMap.at(0).size()));

	return (*std::max_element(sigmaMap.begin(), sigmaMap.end(),
							  [modelVectorIndex](const auto &a, const auto &b)
							  {
								  return a[modelVectorIndex] < b[modelVectorIndex];
							  }))[modelVectorIndex];
}

float Som::getMinSigmaOfFeature(size_t modelVectorIndex) const
{
	assert(modelVectorIndex < static_cast<size_t>(sigmaMap.at(0).size()));

	return (*std::min_element(sigmaMap.begin(), sigmaMap.end(),
							  [modelVectorIndex](const auto &a, const auto &b)
							  {
								  return a[modelVectorIndex] < b[modelVectorIndex];
							  }))[modelVectorIndex];
}

Som::Metrics Som::getMetrics() const noexcept
{
	return metrics;
}

bool Som::isTraining() const noexcept
{
	return _isTraining;
}

bool Som::isCompatibleWithData(DataSet &data) const noexcept
{
	return transform.Length(data.vectorLength()) == depth;
}

SomIndex Som::findBmu(const Eigen::VectorXf &v) const
{
	auto valid = Eigen::VectorXf::Ones(v.size());
	auto weights = Eigen::VectorXf::Ones(v.size());
	
	return findBmu(v, valid, weights);
}

SomIndex Som::findBmu(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, const Eigen::VectorXf &weights) const
{
	auto minDist = this->euclidianWeightedDist(0, v, valid, weights);
	size_t minIndex = 0;

	for (size_t i = 0; i < height * width; ++i)
	{
		double currentDist;
		if ((currentDist = this->euclidianWeightedDist(i, v, valid, weights)) < minDist)
		{
			minDist = currentDist;
			minIndex = i;
		}
	}

	SomIndex returnIndex(minIndex % width, minIndex / width);

	return returnIndex;
}

// Search through whole map where bmuHits > minBmuHits and return index of neuron that is closest to v.
// Only considers valid dimensions with weights
SomIndex Som::findRestrictedBmu(const Eigen::VectorXf &v, const Eigen::VectorXf &valid,
								const size_t minBmuHits, const Eigen::VectorXf &weights) const
{
	auto minDist = this->euclidianWeightedDist(0, v, valid, weights);
	size_t minIndex = 0;

	for (size_t i = 0; i < height * width; ++i)
	{
		double currentDist;
		if ((currentDist = this->euclidianWeightedDist(i, v, valid, weights)) < minDist && bmuHits[i] >= minBmuHits)
		{
			minDist = currentDist;
			minIndex = i;
		}
	}

	SomIndex returnIndex(minIndex % width, minIndex / width);

	return returnIndex;
}

// Searches the neighbourhood of map starting from neuron lastBMU and returns index of neuron closest to v
SomIndex Som::findLocalBmu(const Eigen::VectorXf &v, const Eigen::VectorXf &valid,
						   const size_t &lastBMUref, const Eigen::VectorXf &weights) const
{
	size_t lastBMU = lastBMUref;
	double minDist = this->euclidianWeightedDist(lastBMU, v, valid, weights);
	size_t minIndex = lastBMU;
	constexpr std::array<size_t, 8> firstSearchX{-1uz, 0uz, 1uz, 1uz, 1uz, 0uz, -1uz, -1uz};
	constexpr std::array<size_t, 8> firstSearchY{1uz, 1uz, 1uz, 0uz, -1uz, -1uz, -1uz, 0uz};

	size_t lastMeasured = lastBMU;

	size_t lastMeasuredX, lastMeasuredY, lastBMUX, lastBMUY;

	size_t currentX, currentY;
	size_t startX, endX;

	bool success = false;

	while (!success)
	{
		lastMeasuredX = lastMeasured % width;
		lastMeasuredY = lastMeasured / width;

		lastBMUX = lastBMU % width;
		lastBMUY = lastBMU / width;

		// If first try
		if (lastMeasured == lastBMU)
		{
			for (size_t i{0uz}; i < 8uz; ++i)
			{
				currentX = std::max(std::min(lastMeasuredX + firstSearchX[i], width - 1uz), 0uz);
				currentY = std::max(std::min(lastMeasuredY + firstSearchY[i], height - 1uz), 0uz);
				
				if (double currentDist; (currentDist = this->euclidianWeightedDist(currentY * width + currentX, v, valid, weights)) < minDist)
				{
					minDist = currentDist;
					minIndex = currentY * width + currentX;
				}
			}

			// If BMU has not changed. Do no more
			if (minIndex == lastBMU)
			{
				return SomIndex{minIndex % width, minIndex / width};
			}
			else
			{
				lastMeasured = minIndex;
			}
		}
		// If not first try
		else
		{
			// If we are moving in X direction
			if (lastMeasuredX - lastBMUX)
			{
				for (int i = -1; i < 2; ++i)
				{
					currentX = std::max(std::min(lastMeasuredX + lastMeasuredX - lastBMUX, width - 1uz), 0uz);
					currentY = std::max(std::min((lastMeasuredY + i), height - 1uz), 0uz);
					
					if (double currentDist; (currentDist = this->euclidianWeightedDist(currentY * width + currentX, v, valid, weights)) < minDist)
					{
						minDist = currentDist;
						minIndex = currentY * width + currentX;
					}
				}
			}

			// If we are moving in Y direction
			if (lastMeasuredY - lastBMUY)
			{
				// When moving in positive Y direction, we have already measured distance from X-point + 1
				if (lastMeasuredX - lastBMUX > 0)
				{
					startX = -1uz;
					endX = 0uz;
				}
				// When moving in negative Y direction, we have already measured distance from X-point - 1
				else if (lastMeasuredX - lastBMUX < 0)
				{
					startX = 0uz;
					endX = 1uz;
				}
				// When not moving at all in Y direction, we have to measure all three points in X range from -1 to 1
				else
				{
					startX = -1uz;
					endX = 1uz;
				}
				for (size_t i = startX; i < (endX + 1uz); ++i)
				{
					currentX = std::max(std::min((lastMeasuredX + i), width - 1uz), 0uz);
					currentY = std::max(std::min(lastMeasuredY + lastMeasuredY - lastBMUY, height - 1uz), 0uz);
					
					if (double currentDist; (currentDist = this->euclidianWeightedDist(currentY * width + currentX, v, valid, weights)) < minDist)
					{
						minDist = currentDist;
						minIndex = currentY * width + currentX;
					}
				}
			}

			// If BMU did not change since last try, do no more
			if (minIndex == lastMeasured)
			{
				return SomIndex{minIndex % width, minIndex / width};
			}
			// Otherwise update BMUs
			else
			{
				lastBMU = lastMeasured;
				lastMeasured = minIndex;
			}
		}
	}

	return SomIndex{minIndex % width, minIndex / width};
}

/*				find restricted Best Matching Distribution				*/
std::vector<double> Som::findRestrictedBmd(const Eigen::VectorXf &v, const Eigen::VectorXf &valid,
										   size_t minBmuHits, const Eigen::VectorXf &weights) const
{
	std::vector<double> dist(height * width, -1);

	// Normalization constant
	double C = 0;

	for (size_t i = 0; i < height * width; ++i)
	{
		// Restriction
		if (bmuHits[i] >= minBmuHits)
		{
			dist[i] = this->euclidianWeightedDist(i, v, valid, weights);

			// Transform euclidian distance into a probability measure
			dist[i] = (double)std::exp(-dist[i] * dist[i] / 2);

			// Integrate all probabilities
			C += dist[i];
		}
		else
			dist[i] = 0;
	}

	// Normalize with respect to C in order to get a distribution
	for (size_t i = 0; i < height * width; ++i)
		dist[i] /= C;

	return dist;
}

// Evaluates mean error of a map with dataset data
double Som::evaluate(const DataSet &data) const
{
	double error = 0;

	Eigen::VectorXf binaryError;
	Eigen::ArrayXf ones(map[0].size(), 1);

	auto continuous = data.getContinuous();

	auto binary = data.getBinary();

	auto binaryEigen = binary.cast<float>();

	for (size_t i = 0; i < data.size(); i++)
	{
		auto val = (data.getValidity(i).array() * continuous).cast<float>();

		SomIndex bmu = this->findBmu(data.getData(i), val.cast<float>(), data.getWeights());

		binaryError = (((this->getNeuron(bmu).array().log() * data.getData(i).array()) + ((ones - this->getNeuron(bmu).array()).log() * (ones - data.getData(i).array())))).array();

		// Replace NaNs and Infs with something big
		binaryError = (binaryError.array().isNaN() || binaryError.array().isInf()).select(-99999, binaryError);

		binaryError = (binaryError.array() * binaryEigen.array() * val.array()).matrix();
		/*	Ref.
		 * Incremental calculation of weighed mean and variance. Tony Finch. University of Cambrige Computing Service. February 2009.
		 * Equation 4
		 * Mean of binaryError + continious error */
		error += 1. / (static_cast<double>(i) + 1.0) * (this->euclidianWeightedDist(bmu, data.getData(i), val.cast<float>(), data.getWeights()) + std::sqrt(binaryError.dot(binaryError)) - error);
	}

	return error;
}

size_t Som::variationalAutoEncoder(const DataSet *data, size_t minBmuHits) const
{
	std::random_device rd;
	std::mt19937 gen(rd());

	size_t modelVector{0};

	for (size_t i = 0; i < data->size(); ++i)
	{

		// Extract sample vector
		Eigen::VectorXf v = data->getData(i);

		Eigen::VectorXf val = data->getValidity(i).cast<float>();

		auto w = data->getWeights();

		// Find distribution for sample vector
		auto probability = findRestrictedBmd(v, val, minBmuHits, w);

		std::discrete_distribution<size_t> d(probability.begin(), probability.end());

		// TODO: This only outputs a model vector for the last row of the dataset
		modelVector = d(gen);

		// std::cout << "Model vector: "<< modelVector << " is chosen with probability: " << probability[modelVector] << "\n";

		// Print in an Octave-friendly format
		// std::cout << "prob" << i << "= [";
		// for(unsigned int c = 0; c < width; c++)
		//{
		//	for(unsigned int d = 0; d < height; d++)
		//	{
		//		std::cout << probability[d*width + c] << " ";
		//	}
		//	if(c < width -1) std::cout << ";";
		//}
		// std::cout << "]\n";
	}

	return modelVector;
}

int Som::autoEncoder(const DataSet *data, size_t minBmuHits) const
{
	bool success = true;

	std::srand((unsigned)(time(NULL) + clock()));

	for (size_t i = 0; i < data->size(); i++)
	{

		// Extract sample vector
		Eigen::VectorXf v = data->getData(i);

		// std::vector<int> val = data->getValidity(i);

		// std::vector<float> w = data->getWeights();

		// Now we're just getting the bmu to sample from that univariate normal distribution.
		// What we could do is to take a distribution over model vectors (use variationalAutoEncoder function),
		// sample from that distribution to get a stocastic model vector and then use that to sample
		// from its univariate normal distribution

		// Fixed

		auto bmuInt = variationalAutoEncoder(data, minBmuHits);

		SomIndex bmu(bmuInt % width, bmuInt / width);

		// Find best matching unit (bmu)
		// SomIndex bmu = findRestrictedBmu(v, val, minBmuHits, w);

		for (int n = 0; n < v.size(); n++)
		{
			std::cout << v(n) << "\n";

			{ // This value should be sampled from a normal distribution with mean this->getNeuron(bmu)[n] and standard deviation this->getSigmaNeuron(bmu)[n] if continuous
				// and from a probability of (1-n)*this->getNeuron(bmu)[n] + n*(1-this->getNeuron(bmu)[n]) if binary
				//
				// Approximate normal distribution by sampling from L(0,1) and calculate logit(L)/1.6*sigma + mean = log(L/(1-L))/1.6*sigma + mean
				//

				double L = (double)(std::rand() % (int)(1000)) / 1000;
				double N = std::log(L / (1 - L)) / 1.6 * this->getSigmaNeuron(bmu)[n] + this->getNeuron(bmu)[n];

				std::cout << data->getName(n) << "\t" << N << "\t"
						  << "\n";
				//
				// As for now, we're just printing a deterministic value of this->getNeuron(bmu)[n]
				// std::cout << data->getName(n) << "\t" << this->getNeuron(bmu)[n] << "\n";
			}
		}

		std::cout << "\n";
	}

	return success;
}

// Find BMU of each record in DataSet data and calculate normalized distance between them.
// In non-verbose mode. Only print distances between record and BMU with biggest distance in any dimension
// COLUMN_NAME \t Relative distance \t Record value \t Min of approved interval \t Max of approved interval
//
// In verbose mode. Print distances of all columns of all records with format:
// COLUMN_NAME \t Record value \t Min of approved interval \t Max of approved interval
int Som::measureSimilarity(const DataSet *data, int numOfSigmas, size_t minBmuHits) const
{
	bool success = true;
	float maxValue{-99999999.f};
	size_t maxValueDataSetRow{0uz};
	bool last{false};

	auto binaryEigen = data->getBinary();
	auto contEigen = data->getContinuous();

	for (size_t i = 0; i < data->size() + 1; i++)
	{
		if (i == data->size())
		{
			i = maxValueDataSetRow;
			last = true;
		}

		// Extract sample vector
		Eigen::VectorXf v = data->getData(i);

		// Find best matching unit (bmu)
		SomIndex bmu = findRestrictedBmu(v, data->getValidity(i).cast<float>(), minBmuHits, data->getWeights());

		auto pos = width * bmu.getY() + bmu.getX();

		// Modified sigma vector
		Eigen::VectorXf sM = (sigmaMap[pos].array() > 0.00001f).select(0.00001f, sigmaMap[pos]);

		// Make map binary in order to get a proper interval (min - max) when comparing binaries
		Eigen::ArrayXf binaryMap = (binaryEigen.array() > 0).select(0.0f, map[pos]);
		binaryMap = (binaryEigen.array() > 0 && this->getNeuron(bmu).array() >= 0.5)
						.select(Eigen::VectorXf::Constant(binaryMap.size(), 1.0f), 0.0f);

		auto validEigen = data->getValidity(i);

		// Calculate delta vector
		// Eigen::VectorXf delta = (validEigen.array()*((contEigen.array()*(v - this->getNeuron(bmu)).array()/sM.array()/numOfSigmas) +
		//(binaryEigen.array()*(v.array()*(this->getNeuron(bmu).array().log())+(ones-v.array())*((ones-this->getNeuron(bmu).array()).log())))) ).matrix();

		Eigen::ArrayXf delta = /*contEigen.array()*/ (v - this->getNeuron(bmu)).array() / sM.array() / numOfSigmas;
		// Eigen::ArrayXf binOnLoss = binaryEigen.array()*( v.array()*(this->getNeuron(bmu).array().log()) );
		// Eigen::ArrayXf binOffLoss = binaryEigen.array()*( (v-this->getNeuron(bmu)).array()*((ones - this->getNeuron(bmu).array()).log()) );

		Eigen::ArrayXf min = ((/*contEigen.array()*/ ((this->getNeuron(bmu) - sM * numOfSigmas).array()))); // +
		//(binaryEigen.array()*(binaryMap - 0.0001*ones)));
		Eigen::ArrayXf max = ((/*contEigen.array()*/ ((this->getNeuron(bmu) + sM * numOfSigmas).array()))); // +
		//(binaryEigen.array()*(binaryMap + 0.0001*ones)));

		// std::cout << "# BMU hits: " << bmuHits[pos] << "\tX: " << bmu.getX() << "\tY: " << bmu.getY() << "\n";
		// if( ( ( (data->size()) > 100 && (i % (int)((data->size())/100)) == 0 ) || (data->size()) < 100 ) && _verbose )
		//	std::cout << "\rMeasuring dataset:" << 100*i/(data->size()) << "%";

		for (int n = 0; n < v.size(); n++)
		{
			if (delta(n) > maxValue)
			{
				maxValue = static_cast<float>(fabs(delta(n)));
				maxValueDataSetRow = i;
			}

			// Replace NaN and Inf with zero (0)
			delta(n) = std::isnan(delta(n)) || std::isinf(delta(n)) ? 0 : delta(n);
			// binOnLoss(n) = std::isnan(binOnLoss(n)) ? 0 : binOnLoss(n);
			// binOffLoss(n) = std::isnan(binOffLoss(n)) ? 0 : binOffLoss(n);

			if (validEigen(n) && last)
			{
				// TODO: Make this some kind of return instead
				// std::cout << data->getName(n) << "\t" << delta(n)/*+binOnLoss(n)+binOffLoss(n)*/ << "\t" << v(n) << "\t" << min(n) << "\t" << max(n) << "\n";

				if (((v(n) < min(n) || v(n) > max(n)) && validEigen(n)))
					success = false;
			}
		}

		if (last)
			break;

		// std::cout << "\n";
	}

	return success;
}

void Som::trainBatchSom(DataSet &data, size_t numberOfEpochs, double sigma0, double sigmaDecay, bool updateUMatrixAfterEpoch)
{
	/* Reset metrics */
	metrics = Som::Metrics(numberOfEpochs);

	auto weights = data.getWeights();

	for (size_t i = 0; i < numberOfEpochs; ++i)
	{
		std::cout << "Training VSOM epoch " << i << "/" << numberOfEpochs << '\n';
		
		auto sigma = sigma0 * std::exp(-sigmaDecay * static_cast<double>(i));

		if (sigma < 1.0)
			return;


		auto meanSquareError = float{0.0f};
		auto countDataChunks = size_t{0};
		while(!data.hasReadWholeDataStream())
		{
			data.loadNextDataFromStream();
			meanSquareError += trainBatchSomEpoch(data, sigma, i == 0);

			++countDataChunks;
		}

		meanSquareError /= static_cast<float>(countDataChunks);
		{
			const std::lock_guard<std::mutex> lock(metricsMutex);
			metrics.MeanSquaredError[i] = meanSquareError;
		}

		data.resetStreamLoadPosition();

		if (updateUMatrixAfterEpoch)
			updateUMatrix(data.getWeights());
	}
}

float Som::trainBatchSomEpoch(DataSet &dataset, double currentSigma, bool isFirst)
{
	auto meanSquareError = std::atomic<float>{0.0f};
	auto wholeDataset = dataset.getAll();
	auto epochSize = wholeDataset.size();
	// Assign a som neuron to each data row
	if (isFirst)
	{
		std::for_each(
			std::execution::par_unseq,
			wholeDataset.begin(),
			wholeDataset.end(),
			[this, &dataset, &meanSquareError, epochSize](auto data)
			{
				auto valid = Eigen::Map<const Eigen::VectorXi>(data.valid->data(), data.valid->size()).cast<float>();
				auto index = findBmu(
								 *data.data,
								 valid,
								 dataset.getWeights())
								 .getSomIndex(*this);

				*data.lastBMU = index;
				bmuHits[index] += 1uz;

				auto residual = this->transform.Comparer(*data.data, this->getNeuron(index), this->SMap[index], valid);
				meanSquareError.fetch_add(residual.squaredNorm() / static_cast<float>(epochSize));
			});
	}
	else
	{
		std::for_each(
			std::execution::par_unseq,
			wholeDataset.begin(),
			wholeDataset.end(),
			[this, &dataset, &meanSquareError, epochSize](auto data)
			{
				auto valid = Eigen::Map<const Eigen::VectorXi>((*data.valid).data(), (*data.valid).size()).cast<float>();
				auto index = this->findLocalBmu(
									 *data.data,
									 valid,
									 *data.lastBMU,
									 dataset.getWeights())
								 .getSomIndex(*this);

				*data.lastBMU = index;
				bmuHits[index] += 1uz;

				auto residual = this->transform.Comparer(*data.data, this->getNeuron(index), this->SMap[index], valid);
				meanSquareError.fetch_add(residual.squaredNorm() / static_cast<float>(epochSize));
			});
	}

	// Calculate all neurons' new features based on their respective data points
	std::for_each(
		std::execution::par_unseq,
		map.begin(),
		map.end(),
		[this, &dataset, currentSigma](auto &neuron)
		{
			int index = static_cast<int>(&neuron - map.data());
			auto somIndex = SomIndex(*this, index);
			auto data = dataset.getAll();

			auto currentX = somIndex.getX();
			auto currentY = somIndex.getY();

			auto bmus = dataset.getLastBMU();

			/* Calcualte new model vector m_i as mean of all training data x_j
			 *	weighted by their respective neighborhood weight h_ij
			 *	(that's dependent on their respective BMU)
			 *
			 *		   S_j(h_ij * x_j)
			 *	m_i = -----------------
			 *			   S_j(h_ij)
			 *
			 * Self-Organizing Maps T. Kohonen ISBN 3-540-67921-9 3rd Edition
			 * Eq. 3.29
			 */

			/*	Ref.
			 * Incremental calculation of weighed mean and variance. Tony Finch. University of Cambrige Computing Service. February 2009.
			 * */

			auto sumOfWeights = float{0.f};
			auto currentModel = Eigen::VectorXf(getDepth());
			auto currentModelSigma = Eigen::ArrayXf(getDepth());
			currentModel.setZero();
			currentModelSigma.setZero();
			for (const auto &datapoint : data)
			{
				auto bmuIndex = SomIndex(*this, *datapoint.lastBMU);
				auto bmuX = bmuIndex.getX();
				auto bmuY = bmuIndex.getY();

				float currentWeight = static_cast<float>(this->calculateNeighbourhoodWeight(
					currentX, currentY, bmuX, bmuY, currentSigma));

				auto valid = Eigen::Map<const Eigen::VectorXi>((*datapoint.valid).data(), (*datapoint.valid).size()).cast<float>();
				
				/* Eq. 47 */
				sumOfWeights += currentWeight;

				auto lastModel = currentModel;

				Eigen::VectorXf currentDelta = this->transform.Stepper(*(datapoint.data), currentModel, valid);

				/* Eq. 53 */
				currentModel = currentModel + currentWeight / sumOfWeights * currentDelta;

				/* Eq. 68 */
				currentModelSigma = currentModelSigma + currentWeight * (this->transform.Stepper(*(datapoint.data), lastModel, valid)).array() * (currentDelta.array());
			}

			neuron = currentModel;
	
			/* Eq. 69 */
			sigmaMap[index] = (currentModelSigma / sumOfWeights).sqrt().matrix();
			
			weightMap[index] = sumOfWeights;
		});
	
	return meanSquareError;
}

// Trains on a single data point, i.e. a single record vector v is used to update
// the map weights with parameters valid, weights, eta, sigma and weightDecayFunction
// lastBMU is the index of the BMU that this specific record found last epoch.
// It's used as a starting point if we're looking for a local BMU.
Som::TrainingReturnValue Som::trainSingle(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, const Eigen::VectorXf &weights,
										  const double eta, const double sigma, size_t &lastBMU, const WeigthDecayFunction weightDecayFunction)
{
	// Find best matching unit (bmu)
	const auto bmu = [this, sigma, v, valid, weights, lastBMU]()
	{
		return sigma > SIGMA_SWITCH_TO_LOCAL ? findBmu(v, valid, weights) : findLocalBmu(v, valid, lastBMU, weights);
	}();

	// Save last BMU to dataset
	lastBMU = bmu.getY() * width + bmu.getX();

	// Size of neighbourhood
	// Only calculate new values within +- 2.5 standard deviations of neighbourhood
	size_t startX = static_cast<size_t>(std::max((static_cast<double>(bmu.getX()) - 2.5 * sigma), 0.));
	size_t startY = static_cast<size_t>(std::max((static_cast<double>(bmu.getY()) - 2.5 * sigma), 0.));

	size_t endX = static_cast<size_t>(std::min((static_cast<double>(bmu.getX()) + 2.5 * sigma), static_cast<double>(width)));
	size_t endY = static_cast<size_t>(std::min((static_cast<double>(bmu.getY()) + 2.5 * sigma), static_cast<double>(height)));

	Eigen::VectorXf totalWeight = valid.array() * weights.array();

	for (size_t j = startY; j < endY; j++)
	{
		for (size_t i = startX; i < endX; i++)
		{
			// Calculate delta vector
			const Eigen::VectorXf delta = transform.Stepper(v, map[j * width + i], totalWeight).matrix();

			// Strength of neighbourhood. Calculated for each neuron
			auto neighbourhoodWeight = calculateNeighbourhoodWeight(i, j, bmu.getX(), bmu.getY(), sigma);

			/*	Ref.
			 * Incremental calculation of weighed mean and variance. Tony Finch. University of Cambrige Computing Service. February 2009.
			 * */

			// Exponential weight decay function
			if (weightDecayFunction == WeigthDecayFunction::Exponential)
			{
				weightMap[j * width + i] += static_cast<float>(neighbourhoodWeight * eta);
				map[j * width + i] += neighbourhoodWeight * eta * delta;
			}
			// Inverse proportional weight decay function
			else
			{
				weightMap[j * width + i] += static_cast<float>(neighbourhoodWeight); // Equation 47

				// Protection against division by zero
				auto tempWeight = weightMap[j * width + i] == 0 ? 1.0 : neighbourhoodWeight / weightMap[j * width + i];

				map[j * width + i] += tempWeight * transform.Stepper(v, map[j * width + i], totalWeight); // Equation 53
			}

			// Protection against division by zero
			auto tempWeight = weightMap[j * width + i] == 0 ? 0.000001 : weightMap[j * width + i];

			SMap[j * width + i] += neighbourhoodWeight * (delta.array() * transform.Stepper(v, map[j * width + i], totalWeight).array()).matrix(); // Equation 68
			sigmaMap[j * width + i] = (SMap[j * width + i].array() / tempWeight).abs().sqrt().matrix();				  // Equation 69
		}
	}
	// Return best matching unit
	return TrainingReturnValue{bmu, transform.Comparer(v, getNeuron(bmu), getSigmaNeuron(bmu), totalWeight), static_cast<float>(euclidianWeightedDist(bmu, v, valid, weights))};
}

double Som::calculateNeighbourhoodWeight(
	const size_t &currentX, const size_t &currentY,
	const size_t &bmuX, const size_t &bmuY,
	const double &currentSigma)
{
	// Calculate neighbourhood fundction
	if (currentSigma > 1.0)
	{
		double currentX_d = static_cast<double>(currentX);
		double currentY_d = static_cast<double>(currentY);
		double bmuX_d = static_cast<double>(bmuX);
		double bmuY_d = static_cast<double>(bmuY);

		return std::exp(-((currentX_d - bmuX_d) * (currentX_d - bmuX_d) / 2.0 / currentSigma / currentSigma +
						  (currentY_d - bmuY_d) * (currentY_d - bmuY_d) / 2.0 / currentSigma / currentSigma));
	}

	// If sigma is equal or smaller than 1, only update bmu with weight = 1 from neighbourhood function
	else if (currentX == bmuX && currentY == bmuY)
	{
		return 1.0;
	}
	else
	{
		return 0.0;
	}
}

void Som::randomInitialize(int seed, float sigma)
{
	std::srand(seed);

	metrics = Metrics{depth};

	for (size_t i = 0; i < width * height; i++)
	{
		// Initialize each element in each neuron to a value between -sigma and sigma
		for (int n = 0; n < map[i].size(); n++)
		{
			map[i](n) = (static_cast<float>(std::rand() % static_cast<int>((2000 * sigma))) - (1000.f * sigma)) / 1000.f;
			sigmaMap[i](n) = 0.0f;
			SMap[i](n) = 0.0f;
		}
		weightMap[i] = 0.0f;
		bmuHits[i] = 0uz;
		uMatrix[i] = 0.0;
		// std::cout << "InitV[" << i << "]: " << map[i] << "\n";
	}
}

void Som::updateUMatrix(const Eigen::VectorXf &argWeights)
{
	std::vector<double> U(width * height);
	auto val = Eigen::VectorXf::Constant(depth, 1);
	auto weights = Eigen::VectorXf::Constant(depth, 1);
	double diagonalFactor = 0.3;
	double min = 10000000,
		   max = 0;
	// span;

	// For each vector in the map, calculate mean distance to closest neighbouring vectors.
	// Put these mean distances in the U matrix
	for (size_t i = 0; i < height; ++i)
	{
		for (size_t j = 0; j < width; ++j)
		{
			if (j > 0 && i > 0 && j < (width - 1) && i < (height - 1)) // Middle part of map
			{
				U[i * width + j] = (this->euclidianWeightedDistRaw(i * width + j, map[(i + 0) * width + j - 1], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 0) * width + j + 1], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j - 0], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j - 0], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j - 1], val, weights) * diagonalFactor +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j - 1], val, weights) * diagonalFactor +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j + 1], val, weights) * diagonalFactor +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j + 1], val, weights) * diagonalFactor) /
								   8;
			}
			else if (i == 0 && j > 0 && j < (width - 1)) // Upper edge
			{
				U[i * width + j] = (this->euclidianWeightedDistRaw(i * width + j, map[(i + 0) * width + j - 1], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 0) * width + j + 1], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j - 0], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j - 1], val, weights) * diagonalFactor +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j + 1], val, weights) * diagonalFactor) /
								   5;
			}
			else if (i == (height - 1) && j > 0 && j < (width - 1)) // Lower edge
			{
				U[i * width + j] = (this->euclidianWeightedDistRaw(i * width + j, map[(i + 0) * width + j - 1], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 0) * width + j + 1], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j - 0], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j - 1], val, weights) * diagonalFactor +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j + 1], val, weights) * diagonalFactor) /
								   5;
			}
			else if (j == 0 && i > 0 && i < (height - 1)) // Left edge
			{
				U[i * width + j] = (this->euclidianWeightedDistRaw(i * width + j, map[(i + 0) * width + j + 1], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j - 0], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j - 0], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j + 1], val, weights) * diagonalFactor +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j + 1], val, weights) * diagonalFactor) /
								   5;
			}
			else if (j == (width - 1) && i > 0 && i < (height - 1)) // Right edge
			{
				U[i * width + j] = (this->euclidianWeightedDistRaw(i * width + j, map[(i + 0) * width + j - 1], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j - 0], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j - 0], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j - 1], val, weights) * diagonalFactor +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j - 1], val, weights) * diagonalFactor) /
								   5;
			}
			else if (j == 0 && i == 0) // Top left corner
			{
				U[i * width + j] = (this->euclidianWeightedDistRaw(i * width + j, map[(i + 0) * width + j + 1], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j - 0], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j + 1], val, weights) * diagonalFactor) /
								   3;
			}
			else if (j == (width - 1) && i == 0) // Top rigth corner
			{
				U[i * width + j] = (this->euclidianWeightedDistRaw(i * width + j, map[(i + 0) * width + j - 1], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j - 0], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i + 1) * width + j - 1], val, weights) * diagonalFactor) /
								   3;
			}
			else if (j == 0 && i == (height - 1)) // Bottom left corner
			{
				U[i * width + j] = (this->euclidianWeightedDistRaw(i * width + j, map[(i + 0) * width + j + 1], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j - 0], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j + 1], val, weights) * diagonalFactor) /
								   3;
			}
			else if (j == (width - 1) && i == (height - 1)) // Bottom right corner
			{
				U[i * width + j] = (this->euclidianWeightedDistRaw(i * width + j, map[(i + 0) * width + j - 1], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j - 0], val, weights) +
									this->euclidianWeightedDistRaw(i * width + j, map[(i - 1) * width + j - 1], val, weights) * diagonalFactor) /
								   3;
			}
			else
			{
				U[i * width + j] = 0;
			}

			// Find upper and lower bounds in U-matrix
			if (U[i * width + j] > max)
				max = U[i * width + j];
			if (U[i * width + j] < min)
				min = U[i * width + j];
		}
	}

	// Normalize to 0-255 interval
	// span = max - min;
	// if(_verbose) std::cout << "Span: " << span << "\nMin: " << min << "\nMax: " << max << "\n";
	for (size_t i = 0; i < uMatrix.size(); ++i)
	{
		uMatrix[i] = U[i];
	}
}

void Som::train(DataSet &data, size_t numberOfEpochs, double eta0, double etaDecay, double sigma0, double sigmaDecay, WeigthDecayFunction weightDecayFunction, bool updateUMatrixAfterEpoch)
{
	_isTraining = true;
	try
	{
		switch (weightDecayFunction)
		{
		case WeigthDecayFunction::BatchMap:
			trainBatchSom(data, numberOfEpochs, sigma0, sigmaDecay, updateUMatrixAfterEpoch);
			break;
		default:
			trainBasicSom(data, numberOfEpochs, eta0, etaDecay, sigma0, sigmaDecay, weightDecayFunction, updateUMatrixAfterEpoch);
		}
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}
	_isTraining = false;
}

// Loops through the data set for numberOfEpochs turns while changing eta and sigma accordingly and updates map weights on each turn by calling trainSingle
void Som::trainBasicSom(DataSet &data, size_t numberOfEpochs, double eta0, double etaDecay, double sigma0, double sigmaDecay, WeigthDecayFunction weightDecayFunction, bool updateUMatrixAfterEpoch)
{

	/* Reset metrics */
	metrics = Som::Metrics(numberOfEpochs);

	auto weights = data.getWeights();

	for (size_t i = 0; i < numberOfEpochs; ++i)
	{
		auto eta = eta0 * std::exp(-etaDecay * static_cast<double>(i));
		auto sigma = sigma0 * std::exp(-sigmaDecay * static_cast<double>(i));

		if (sigma < 1.0)
			sigma = 1.0;

		std::cout << "Epoch: " << i + 1 << "/" << numberOfEpochs << "\teta: " << eta << "\tsigma: " << sigma << "\n";

		float meanSquareError{0.0};
		size_t countDataChunks{0};
		while(!data.hasReadWholeDataStream())
		{
			data.loadNextDataFromStream();

			auto epochSize = data.size();

			for (size_t j = 0; j < epochSize; ++j)
			{
				auto [pos, residual, distanceError] = trainSingle(data.getData(j), data.getValidity(j).cast<float>(), weights, eta, sigma, data.getLastBMU(j), weightDecayFunction);

				addBmu(pos);

				meanSquareError += residual.squaredNorm() / static_cast<float>(epochSize);

				if ((data.size() > 100 && (j % (int)(data.size() / 100)) == 0) || data.size() < 100)
					std::cout << "\rTraining SOM:" << 100 * j / data.size() << "%";
			}

			++countDataChunks;
		}
		meanSquareError /= static_cast<float>(countDataChunks);
		{
			const std::lock_guard<std::mutex> lock(metricsMutex);
			metrics.MeanSquaredError[i] = meanSquareError;
		}

		data.resetStreamLoadPosition();

		if (updateUMatrixAfterEpoch)
			updateUMatrix(data.getWeights());
	}
	std::cout << "\rTraining SOM:100%\n";
}

void Som::addBmu(SomIndex pos)
{
	bmuHits[this->getIndex(pos)] += 1uz;
}

void Som::displayUMatrix() const
{
	std::cout << "\nU-matrix:\n[";
	for (unsigned int i = 0; i < height; i++)
	{
		for (unsigned int j = 0; j < width; j++)
		{
			std::cout << (double)uMatrix[i * width + j] << " ";
		}
		std::cout << "; ";
	}
	std::cout << "]\n";
}

// The SOM is saved in an octave/matlab friendly file that holds all variables needed to analyze what's been learned in octave/matlab
void Som::save(const char *fileName) const
{
	FILE *fp;

	if ((fp = fopen(fileName, "w")) == NULL)
	{
		std::cout << "Could not open file " << fileName << " for writing. Quitting...\n";
		exit(EXIT_FAILURE);
	}

	fprintf(fp, "# This file is generated by Som.exe\n# It contains all data needed to evaluate a sample according to the SOM trained by Som.exe\n");
	fprintf(fp, "# name: som\n");
	fprintf(fp, "# type: matrix\n");
	fprintf(fp, "# ndims: 3\n");
	fprintf(fp, " %lu %lu %lu\n", this->getWidth(), this->getHeight(), (long unsigned int)map[0].rows());
	// fprintf(fp, "Vector length:%lu\n", (long unsigned int)map[0].rows());
	// fprintf(fp, "Som width:%u\nSom height:%u\n", this->getWidth(), this->getHeight());
	// fprintf(fp, "Som data:\n");

	// Print map data
	for (int i = 0; i < map[0].size(); i++)
	{
		for (unsigned int j = 0; j < height * width; j++)
		{
			fprintf(fp, " %f\n", map[j](i));
		}
		// fprintf(fp, "\n");
	}

	fprintf(fp, "\n\n# name: sigmaSom\n");
	fprintf(fp, "# type: matrix\n");
	fprintf(fp, "# ndims: 3\n");
	fprintf(fp, " %lu %lu %lu\n", this->getWidth(), this->getHeight(), (long unsigned int)sigmaMap[0].size());
	// fprintf(fp, "Sigma SOM data:\n");
	//  Print map data
	for (int i = 0; i < sigmaMap[0].size(); i++)
	{
		for (unsigned int j = 0; j < height * width; j++)
		{
			fprintf(fp, " %f\n", sigmaMap[j](i));
		}
		// fprintf(fp, "\n");
	}

	fprintf(fp, "\n\n# name: weightMap\n");
	fprintf(fp, "# type: matrix\n");
	fprintf(fp, "# rows: %lu\n", this->getWidth());
	fprintf(fp, "# columns: %lu\n", this->getHeight());
	// fprintf(fp, "Weight som data:\n");
	//  Print bmuHits data
	for (unsigned int j = 0; j < width; j++)
	{
		for (unsigned int i = 0; i < height; i++)
			fprintf(fp, "%f\t", weightMap[j * height + i]);

		fprintf(fp, "\n");
	}

	fprintf(fp, "\n\n# name: bmuHits\n");
	fprintf(fp, "# type: matrix\n");
	fprintf(fp, "# rows: %lu\n", this->getWidth());
	fprintf(fp, "# columns: %lu\n", this->getHeight());
	// fprintf(fp, "\nBmuHits data:\n");
	//  Print bmuHits data
	for (unsigned int j = 0; j < width; j++)
	{
		for (unsigned int i = 0; i < height; i++)
			fprintf(fp, "%lu\t", bmuHits[j * height + i]);
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n\n# name: U\n");
	fprintf(fp, "# type: matrix\n");
	fprintf(fp, "# rows: %lu\n", this->getWidth());
	fprintf(fp, "# columns: %lu\n", this->getHeight());
	// fprintf(fp, "\nU-matrix data:\n");
	//  Print U-matrix data
	for (unsigned int j = 0; j < width; j++)
	{
		for (unsigned int i = 0; i < height; i++)
			fprintf(fp, "%f\t", uMatrix[j * height + i]);
		fprintf(fp, "\n");
	}

	fclose(fp);
}

Eigen::VectorXf Som::getSizeFromFile(const char *fileName)
{
	std::ifstream file(fileName);
	std::string line;
	std::string::size_type ptr;
	std::size_t found;
	unsigned long loadedVectorLength;
	Eigen::VectorXf initV(1);

	if (!file.is_open())
	{
		std::cout << "Could not open file " << fileName << " for reading. Quitting...\n";
		exit(EXIT_FAILURE);
	}

	while (std::getline(file, line))
	{
		// Check vector length
		if (line.compare(0, 9, "# ndims: ") == 0)
		{
			std::getline(file, line);
			found = line.find_last_of(" ");
			loadedVectorLength = std::stoul(line.substr(found + 1), &ptr, 10);
			initV.resize(loadedVectorLength);
			// std::cout << "Found vector length:" << loadedVectorLength << "\n";
		}

		// Find width of som
		else if (line.compare(0, 8, "# rows: ") == 0)
		{
			height = std::stoul(&(line[8]), &ptr, 10);
			continue;
		}

		// Find height of som
		else if (line.compare(0, 11, "# columns: ") == 0)
		{
			width = std::stoul(&(line[11]), &ptr, 10);
			continue;
		}
	}

	file.close();

	return initV;
}

void Som::load(const char *fileName)
{
	std::ifstream file(fileName);
	std::string line;
	std::string::size_type ptr;
	int readSomData = 0;
	int readSigmaSomData = 0;
	int readWeightMapData = 0;
	int readBmuHitsData = 0;
	int readUMatrixData = 0;
	// unsigned long dimension = 0;

	char *columnMarker;

	unsigned int i = 0, d = 0;

	if (!file.is_open() || !file.good())
	{
		std::cout << "Could not open file " << fileName << " for reading. Quitting...\n";
		exit(EXIT_FAILURE);
	}

	while (!file.eof())
	{

		getline(file, line);
		// std::cout << "|" << line << "|\n";

		if (line[0] == '#')
		{
			// std::cout << line << "\n";
			//  Check vector name
			if (line.compare(0, 8, "# name: ") == 0)
			{
				// std::cout << "Found name line: ";
				if (line.compare(8, 20, "som") == 0)
				{
					readSomData = 1;
					readSigmaSomData = 0;
					readWeightMapData = 0;
					readBmuHitsData = 0;
					readUMatrixData = 0;
					// std::cout << "SOM\n";
					continue;
				}
				else if (line.compare(8, 20, "sigmaSom") == 0)
				{
					readSomData = 0;
					readSigmaSomData = 1;
					readWeightMapData = 0;
					readBmuHitsData = 0;
					readUMatrixData = 0;
					// std::cout << "sigmaSom\n";
					continue;
				}
				else if (line.compare(8, 20, "weightMap") == 0)
				{
					readSomData = 0;
					readSigmaSomData = 0;
					readWeightMapData = 1;
					readBmuHitsData = 0;
					readUMatrixData = 0;
					// std::cout << "weightMap\n";
					continue;
				}
				else if (line.compare(8, 20, "bmuHits") == 0)
				{
					readSomData = 0;
					readSigmaSomData = 0;
					readWeightMapData = 0;
					readBmuHitsData = 1;
					readUMatrixData = 0;
					// std::cout << "bmuHits\n";
					continue;
				}
				else if (line.compare(8, 20, "U") == 0)
				{
					readSomData = 0;
					readSigmaSomData = 0;
					readWeightMapData = 0;
					readBmuHitsData = 0;
					readUMatrixData = 1;
					// std::cout << "U\n";
					continue;
				}
				else
				{
					std::cout << "No name for vector in SOM file!\n";
					exit(EXIT_FAILURE);
				}

				i = 0;
				d = 0;
			}

			// Find vector type
			else if (line.compare(0, 8, "# type: ") == 0)
			{
				// std::cout << "Found type line: ";
				if (line.compare(8, 20, "matrix") != 0)
				{
					std::cout << "Incorrect type in SOM file!\n";
					exit(EXIT_FAILURE);
				}
				// std::cout << "matrix\n";
			}

			// Find vector dimension
			else if (line.compare(0, 9, "# ndims: ") == 0)
			{
				// dimension = std::stoul(&(line[9]), &ptr, 10);
				// std::cout << "Found ndim line: " << dimension << "\n";

				getline(file, line);

				// std::cout << "|" << line << "\n";

				// To go: Do something with dimension length

				continue;
			}

			// Find 2D matrix height
			else if (line.compare(0, 8, "# rows: ") == 0)
			{
				width = std::stoul(&(line[8]), &ptr, 10);
				// std::cout << "Found rows line: " << height << "\n";
				continue;
			}

			// Find 2D matrix width
			else if (line.compare(0, 11, "# columns: ") == 0)
			{
				height = std::stoul(&(line[11]), &ptr, 10);
				// std::cout << "Found columns line: " << width << "\n";
				continue;
			}
		}

		// Save SOM data to Eigen vector
		else if (readSomData && (std::isdigit(line[1]) || line[1] == '-'))
		{
			map[d](i) = static_cast<float>(std::stod(line, &ptr));

			// if( i == 0 )
			//	std::cout << "Found som value:" << map[d](i) << " for element:" << i << "\tat pos:" << d << "\n";

			d++;

			if (d >= height * width)
			{
				d = 0;
				i++;
				if (i >= map[0].size())
				{
					i = 0;
					readSomData = 0;
				}
			}

			continue;
		}

		// Save sigmaSOM data to Eigen vector
		else if (readSigmaSomData && (std::isdigit(line[1]) || line[1] == '-'))
		{
			sigmaMap[d](i) = static_cast<float>(std::stod(line, &ptr));
			// std::cout << "Found sigmaSom value:" << sigmaMap[i](d) << " for element:" << d << "\n";
			d++;

			if (d >= height * width)
			{
				d = 0;
				i++;
				if (i >= sigmaMap[0].size())
				{
					i = 0;
					readSigmaSomData = 0;
				}
			}

			continue;
		}

		// Save  weightMap data to Eigen vector
		else if (readWeightMapData && std::isdigit(line.c_str()[0]))
		{
			for (unsigned int l = 0; l < width; l++)
			{
				columnMarker = (char *)line.c_str();

				for (unsigned int k = 0; k < height; k++)
				{
					// printf("|%s\n", columnMarker);
					weightMap(l * height + k) = static_cast<float>(std::strtod(columnMarker, &columnMarker));
					// columnMarker = (char*)ptr;
					// std::cout << weightMap(l*width + k) << "\t";
				}
				// std::cout << "\n";
				getline(file, line);
			}

			continue;
		}

		// Save bmu hits data to vector
		else if (readBmuHitsData && std::isdigit(line.c_str()[0]))
		{
			for (unsigned int l = 0; l < width; l++)
			{
				columnMarker = (char *)line.c_str();

				for (unsigned int k = 0; k < height; k++)
				{
					// printf("|%s\n", columnMarker);
					bmuHits[l * height + k] = std::strtoul(columnMarker, &columnMarker, 10);
					// columnMarker = (char*)ptr;
					// std::cout << bmuHits[l*height + k] << "\t";
				}
				// std::cout << "\n";
				getline(file, line);
			}

			continue;
		}

		// Save U-matrix data to vector
		else if (readUMatrixData && std::isdigit(line.c_str()[0]))
		{
			for (unsigned int l = 0; l < width; l++)
			{
				columnMarker = (char *)line.c_str();

				for (unsigned int k = 0; k < height; k++)
				{
					// printf("|%s\n", columnMarker);
					uMatrix[l * height + k] = std::strtod(columnMarker, &columnMarker);
					// columnMarker = (char*)ptr;
					// std::cout << uMatrix[l*width + k] << "\t";
				}
				// std::cout << "\n";
				getline(file, line);
			}

			continue;
		}
		else
		{
			;
		}

		// std::cout << line << '\n';
	}

	file.close();
}