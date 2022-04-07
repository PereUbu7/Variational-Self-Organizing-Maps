#include "SOM.hpp"

#include <iostream>
#include <fstream>
#include <random>

// Initialize an empty SOM
Som::Som(size_t inWidth = 1, size_t inHeight = 1, size_t inDepth = 1, bool verbose)
{
	_verbose = verbose;
	// Create template vector for whole map
	Eigen::VectorXf initV(inDepth);
	
	// Set map to correct size with every Vector equal to initV
	map.resize(inWidth*inHeight, initV);
	sigmaMap.resize(inWidth*inHeight, initV);
	weightMap.resize(inWidth*inHeight);
	SMap.resize(inWidth*inHeight, initV);
	
	for(size_t i = 0; i < inWidth*inHeight; ++i)
	{
		// Initialize each element in each neuron to 0
		for(int n = 0; n < sigmaMap[i].rows(); ++n)
		{
			map[i](n) = 0;
			sigmaMap[i](n) = 0;
			SMap[i](n) = 0;
		}
		weightMap[i] = 0;
	}
	
	// Set BMU hits map to zero
	bmuHits.resize(inWidth*inHeight, 0);
	
	// Set U-matrix to zero
	uMatrix.resize(inWidth*inHeight, 0);
	
	// Save size of map to object
	width = inWidth;
	height = inHeight;

}

// Initialize SOM from trained data
Som::Som(const char *SOMFileName, bool verbose)
{
	_verbose = verbose;
	// Create template vector for whole map
	Eigen::VectorXf initV(0);
	initV = getSizeFromFile(SOMFileName);
	
	if(_verbose) std::cout << "width: " << width << "\theight: " << height << "\tlength: " << initV.size() << "\n";
	
	// Set map to correct size with every Vector equal to size returned by getSizeFromFile()
	map.resize(width*height, initV);
	sigmaMap.resize(width*height, initV);
	weightMap.resize(width*height);
	SMap.resize(width*height, initV);
	
	for(size_t i = 0; i < width*height; ++i)
	{
		// Initialize each element in each neuron to 0
		for(int n = 0; n < sigmaMap[i].size(); ++n)
		{
			sigmaMap[i](n) = 0;
			SMap[i](n) = 0;
		}
		weightMap[i] = 0;
	}
	
	// Set BMU hits map to zero
	bmuHits.resize(width*height, 0);
	
	// Set U-matrix to zero
	uMatrix.resize(width*height, 0);
	
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
	//Print Bmu hits
	for(size_t i = 0; i < this->getHeight(); ++i)
	{
		for(size_t j = 0; j < this->getWidth(); ++j)
			std::cout << bmuHits[i*width+j] << "\t";
		
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
	//return std::sqrt( (this->getNeuron(pos) - v).dot(this->getNeuron(pos) - v) );
	
	int intPos = pos.getY()*width + pos.getX();
	
	return( euclidianWeightedDist(intPos, v, valid, weights) );
	
}

// Returns the Euclidian squared distance between neuron at position pos and vector v.
// Considers only dimensions where valid is true.
// All dimensions are weighed individually by weights vector
double Som::euclidianWeightedDist(
	const size_t &pos, const Eigen::VectorXf &v, 
	const std::vector<int> &valid, const std::vector<float> &weights) const
{
	// Don't look at dimensions with invalid values
	auto tmp = std::vector<float>(valid.size());
	
	Eigen::ArrayXf sM(sigmaMap[pos].size());
	Eigen::VectorXf validEigen;
	
	for(size_t i = 0; i<valid.size(); ++i)
	{
		tmp[i] = (float)valid[i]*weights[i];
	}
	
	for(size_t i = 0;i<static_cast<size_t>(sigmaMap[pos].size()); ++i)
	{
		sM(i) = std::max((double)sigmaMap[pos][i], 0.00001);
	}
	
	// validEigen is an Eigen vector with weight if valid and 0 if not
	validEigen = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(tmp.data(), valid.size());
	
	// (pos-v)*(pos-v)*weight*valid/sigma2
	return  (((this->getNeuron(pos) - v).array()/sM).matrix().dot( ((this->getNeuron(pos) - v).array()*validEigen.array()/sM).matrix() ) );
}

double Som::euclidianWeightedDist(
	const size_t &pos, const Eigen::VectorXf &v, 
	const Eigen::VectorXf &valid, const Eigen::VectorXf &weights) const
{
	Eigen::ArrayXf sM(sigmaMap[pos].size());
	Eigen::VectorXf validEigen;
	
	sM = (sigmaMap[pos].array() < 0.00001).select(0.00001, sigmaMap[pos]);

	// validEigen is an Eigen vector with weight if valid and 0 if not
	validEigen = valid.array()*weights.array();
	
	// (pos-v)*(pos-v)*weight*valid/sigma2
	return  (((this->getNeuron(pos) - v).array()/sM).matrix().dot( ((this->getNeuron(pos) - v).array()*validEigen.array()/sM).matrix() ) );
}

unsigned int Som::getHeight() const
{
	return height;
}

unsigned int Som::getWidth() const
{
	return width;
}

unsigned int Som::getIndex(SomIndex i) const
{
	return i.getY()*width+i.getX();
}

Eigen::VectorXf Som::getNeuron(SomIndex i) const
{
	return map[i.getY()*width+i.getX()];
}

Eigen::VectorXf Som::getNeuron(size_t i) const
{
	return map[i];
}

Eigen::VectorXf Som::getSigmaNeuron(SomIndex i) const
{
	return sigmaMap[i.getY()*width+i.getX()];
}

Eigen::VectorXf Som::getSigmaNeuron(size_t i) const
{
	return sigmaMap[i];
}

SomIndex Som::findBmu(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, const Eigen::VectorXf &weights) const
{
	auto minDist = this->euclidianWeightedDist(0, v, valid, weights);
	int minIndex = 0;
	
	for(size_t i = 0; i < height*width; ++i)
	{
		double currentDist;
		if( (currentDist = this->euclidianWeightedDist(i, v, valid, weights)) < minDist )
		{
			minDist = currentDist;
			minIndex = i;
		}
	}
	
	SomIndex returnIndex(minIndex % width, minIndex / width );
	
	return returnIndex;
}

// Search through whole map where bmuHits > minBmuHits and return index of neuron that is closest to v.
// Only considers valid dimensions with weights
SomIndex Som::findRestrictedBmu(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, 
	const int minBmuHits, const Eigen::VectorXf &weights) const
{
	auto minDist = this->euclidianWeightedDist(0, v, valid, weights);
	int minIndex = 0;
	
	for(size_t i = 0; i < height*width; ++i)
	{
		double currentDist;
		if( (currentDist = this->euclidianWeightedDist(i, v, valid, weights)) < minDist && bmuHits[i] >= minBmuHits )
		{
			minDist = currentDist;
			minIndex = i;
		}
	}
	
	
	SomIndex returnIndex(minIndex % width, minIndex / width );
	
	return returnIndex;
}

// Search through whole map where bmuHits > minBmuHits and return index of neuron that is closest to v.
// Only considers valid dimensions with weights
SomIndex Som::findRestrictedBmu(const Eigen::VectorXf &v, const std::vector<int> &valid, 
	const int minBmuHits, const std::vector<float> &weights) const
{
	auto minDist = this->euclidianWeightedDist(0, v, valid, weights);
	int minIndex = 0;
	
	for(size_t i = 0; i < height*width; ++i)
	{
		double currentDist;
		if( (currentDist = this->euclidianWeightedDist(i, v, valid, weights)) < minDist && bmuHits[i] >= minBmuHits )
		{
			minDist = currentDist;
			minIndex = i;
		}
	}
	
	
	SomIndex returnIndex(minIndex % width, minIndex / width );
	
	return returnIndex;
}

// Searches the neighbourhood of map starting from neuron lastBMU and returns index of neuron closest to v
SomIndex Som::findLocalBmu(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, 
	const size_t &lastBMUref, const Eigen::VectorXf &weights) const
{
	size_t lastBMU = lastBMUref;
	double minDist = this->euclidianWeightedDist(lastBMU, v, valid, weights);
	size_t minIndex = lastBMU;
	int firstSearchX[8] = {-1, 0, 1, 1, 1, 0, -1, -1};
	int firstSearchY[8] = {1, 1, 1, 0, -1, -1, -1, 0};
	
	size_t lastMeasured = lastBMU;
	
	int lastMeasuredX, lastMeasuredY, lastBMUX, lastBMUY;
	
	int currentX, currentY;
	int startX, endX;
	
	bool success = false;
	
	while(!success)
	{
		lastMeasuredX = lastMeasured % width;
		lastMeasuredY = lastMeasured / width;
		
		lastBMUX = lastBMU % width;
		lastBMUY = lastBMU / width;
		
		//std::cout << "lastMeasuredX:" << lastMeasuredX << "\tlastMeasuredY:" << lastMeasuredY << "\n";
		//std::cout << "lastBMUX:" << lastBMUX << "\tlastBMUY:" << lastBMUY << "\n";
		
		// If first try
		if( lastMeasured == lastBMU )
		{
			for(int i = 0; i < 8; ++i)
			{
				currentX = std::max((int)std::min(lastMeasuredX + firstSearchX[i], (int)width-1), 0);
				currentY = std::max((int)std::min(lastMeasuredY + firstSearchY[i], (int)height-1), 0);
				double currentDist;
				//std::cout << "X:" << currentX << "\tY:" << currentY << "\tindex:" << currentY*width + currentX << "\n";
				if( (currentDist = this->euclidianWeightedDist( currentY*width + currentX , v, valid, weights)) < minDist )
				{
					minDist = currentDist;
					minIndex = currentY*width + currentX;
				}
			}
			
			// If BMU has not changed. Do no more
			if( minIndex == lastBMU )
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
			if( lastMeasuredX - lastBMUX )
			{
				for(int i = -1; i < 2; ++i)
				{
					currentX = std::max((int)std::min(lastMeasuredX + lastMeasuredX - lastBMUX, (int)width-1), 0);
					currentY = std::max((int)std::min(( lastMeasuredY + i ), (int)height-1), 0);
					//std::cout << "X2:" << currentX << "\tY2:" << currentY << "\t2index:" << currentY*width + currentX << "\n";
					double currentDist;
					if( (currentDist = this->euclidianWeightedDist( currentY*width + currentX , v, valid, weights)) < minDist )
					{
						minDist = currentDist;
						minIndex = currentY*width + currentX;
					}
				}
			}
			
			// If we are moving in Y direction
			if( lastMeasuredY - lastBMUY )
			{
				// When moving in positive Y direction, we have already measured distance from X-point + 1
				if( lastMeasuredX - lastBMUX > 0 )
				{
					startX = -1;
					endX = 0;
				}
				// When moving in negative Y direction, we have already measured distance from X-point - 1
				else if( lastMeasuredX - lastBMUX < 0 )
				{
					startX = 0;
					endX = 1;
				}
				// When not moving at all in Y direction, we have to measure all three points in X range from -1 to 1
				else
				{
					startX = -1;
					endX = 1;
				}
				for(int i = startX; i < (endX + 1); ++i)
				{
					currentX = std::max((int)std::min(( lastMeasuredX + i ), (int)width-1), 0);
					currentY = std::max((int)std::min(lastMeasuredY + lastMeasuredY - lastBMUY, (int)height-1), 0);
					//std::cout << "X3:" << currentX << "\tY3:" << currentY << "\t3index:" << currentY*width + currentX << "\n";
					double currentDist;
					if( (currentDist = this->euclidianWeightedDist( currentY*width + currentX , v, valid, weights)) < minDist )
					{
						minDist = currentDist;
						minIndex = currentY*width + currentX;
					}
				}
			}
			
			//std::cout << "minIndex:" << minIndex << "\tlastMeasured:" << lastMeasured << "\n";
			
			// If BMU did not change since last try, do no more
			if( minIndex == lastMeasured )
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
	int minBmuHits, const Eigen::VectorXf &weights) const
{
	std::vector<double> dist(height*width, -1);
	
	// Normalization constant
	double C = 0;
	
	for(size_t i = 0; i < height*width; ++i)
	{
		// Restriction
		if(bmuHits[i] >= minBmuHits)
		{
			dist[i] = this->euclidianWeightedDist(i, v, valid, weights);
			
			// Transform euclidian distance into a probability measure
			dist[i] = (double)std::exp( -dist[i]*dist[i]/2 );
			
			// Integrate all probabilities
			C += dist[i];
		}
		else
			dist[i] = 0;
	}
	
	// Normalize with respect to C in order to get a distribution
	for(size_t i = 0; i < height*width; ++i)
		dist[i] /= C;
	
	return dist;
}

/*				find restricted Best Matching Distribution				*/
std::vector<double> Som::findRestrictedBmd(const Eigen::VectorXf &v, const std::vector<int> &valid, 
	int minBmuHits, const std::vector<float> &weights) const
{
	std::vector<double> dist(height*width, -1);
	
	// Normalization constant
	double C = 0;
	
	for(unsigned int i = 0; i < height*width; i++)
	{
		// Restriction
		if(bmuHits[i] >= minBmuHits)
		{
			dist[i] = this->euclidianWeightedDist(i, v, valid, weights);
			
			// Transform euclidian distance into a probability measure
			dist[i] = (double)std::exp( -dist[i]*dist[i]/2 );
			
			// Integrate all probabilities
			C += dist[i];
		}
		else
			dist[i] = 0;
	}
	
	// Normalize with respect to C in order to get a distribution
	for(unsigned int i = 0; i < height*width; i++)
		dist[i] /= C;
	
	return dist;
}

// Evaluates mean error of a map with dataset data
double Som::evaluate(const DataSet &data) const
{
	double error = 0;
	
	Eigen::VectorXf binaryError;
	Eigen::ArrayXf ones(map[0].size(), 1);

	auto continuous = data.getContinuousEigen();
	
	auto binaryEigen = data.getBinaryEigen().cast<float>();
	
	for(size_t i = 0; i<data.size(); i++)
	{
		auto val = (data.getValidityEigen(i).array()*continuous).cast<float>();
		
		SomIndex bmu = this->findBmu(data.getData(i), val.cast<float>(), data.getWeightsEigen());
		
		binaryError = (((this->getNeuron(bmu).array().log()*data.getData(i).array())+((ones-this->getNeuron(bmu).array()).log()*(ones-data.getData(i).array())))).array();
		
		// Replace NaNs and Infs with something big
		binaryError = (binaryError.array().isNaN() || binaryError.array().isInf()).select(-99999, binaryError);

		binaryError = (binaryError.array()*binaryEigen.array()*val.array()).matrix();
		/*	Ref.
		 * Incremental calculation of weighed mean and variance. Tony Finch. University of Cambrige Computing Service. February 2009.
		 * Equation 4
		 * Mean of binaryError + continious error */
		error += (double)1/(i+1)*( this->euclidianWeightedDist(bmu, data.getData(i), val.cast<float>(), data.getWeightsEigen()) + std::sqrt(binaryError.dot(binaryError)) - error );
		
		//std::cout << "X:" << bmu.getX() << "\tY:" << bmu.getY() << "\tbinError:" << binaryError[6] << "\tv:" << data->getData(i)[6] << "\tmap:" << this->getNeuron(bmu)[6] << "\tbin?:" << binaryEigen[6] << "\tval?:" << validEigen[6] << "\n";
		//std::cout << "\terror: " << error << "\ti+1:" << i+1 << "\tx_(i+1):" << this->euclidianWeightedDist(bmu, data->getData(i), val, data->getWeights()) << "\tbinary error:" << std::sqrt(binaryError.dot(binaryError)) << "\n";
		
		if( ( ( (data.size()) > 100 && (i % (int)((data.size())/100)) == 0 ) || (data.size()) < 100 ) && _verbose )
			std::cout << "\rEvaluating dataset:" << 100*i/(data.size()) << "%";
	}
	
	if(_verbose) std::cout << "\rEvaluating dataset:100%\n";
	
	return error;
}

size_t Som::variationalAutoEncoder(const DataSet *data, int minBmuHits) const
{
	std::random_device rd;
    std::mt19937 gen(rd());
	
	size_t modelVector{0};
	
	for( size_t i = 0; i < data->size(); ++i )
	{
		
		// Extract sample vector
		Eigen::VectorXf v = data->getData(i);
		
		auto val = data->getValidityEigen(i).cast<float>();
		
		auto w = data->getWeightsEigen();
		
		// Find distribution for sample vector
		auto probability = findRestrictedBmd(v, val, minBmuHits, w);
		
		std::discrete_distribution<size_t> d(probability.begin(), probability.end());
		
		// TODO: This only outputs a model vector for the last row of the dataset
		modelVector = d(gen);
		
		// std::cout << "Model vector: "<< modelVector << " is chosen with probability: " << probability[modelVector] << "\n";
		
		// Print in an Octave-friendly format
		//std::cout << "prob" << i << "= [";
		//for(unsigned int c = 0; c < width; c++)
		//{
		//	for(unsigned int d = 0; d < height; d++)
		//	{
		//		std::cout << probability[d*width + c] << " ";
		//	}
		//	if(c < width -1) std::cout << ";";
		//}
		//std::cout << "]\n";
	}
	
	return modelVector;
}

int Som::autoEncoder(const DataSet *data, int minBmuHits) const
{
	bool success = true;
	
	std::srand((unsigned)(time(NULL)+clock()));
	
	for( size_t i = 0; i < data->size(); i++ )
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
		
		int bmuInt = variationalAutoEncoder(data, minBmuHits);
		
		SomIndex bmu(bmuInt % width, bmuInt / width);
		
		// Find best matching unit (bmu)
		//SomIndex bmu = findRestrictedBmu(v, val, minBmuHits, w);
		
		
		if( _verbose )
		{
			std::cout << "X: " << bmu.getX() << "\tY: " << bmu.getY() << "\tID: " << i << "\n";
		}
		
		
		for( int n = 0; n < v.size(); n++ )
		{
			// These ifs have to be sorted out. Why did I condition it on _verbose this way?!?!
			if( !_verbose ) 
			{
				std::cout<< v(n)  << "\n";
			}
			
			if( _verbose )
			{	// This value should be sampled from a normal distribution with mean this->getNeuron(bmu)[n] and standard deviation this->getSigmaNeuron(bmu)[n] if continuous
				// and from a probability of (1-n)*this->getNeuron(bmu)[n] + n*(1-this->getNeuron(bmu)[n]) if binary
				//
				// Approximate normal distribution by sampling from L(0,1) and calculate logit(L)/1.6*sigma + mean = log(L/(1-L))/1.6*sigma + mean
				// 
				
				double L = (double)(std::rand() % (int)(1000))/1000;
				double N = std::log(L/(1-L))/1.6*this->getSigmaNeuron(bmu)[n] + this->getNeuron(bmu)[n];
				
				std::cout << data->getName(n) << "\t" << N << "\t" << "\n";
				//
				// As for now, we're just printing a deterministic value of this->getNeuron(bmu)[n]
				//std::cout << data->getName(n) << "\t" << this->getNeuron(bmu)[n] << "\n";
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
int Som::measureSimilarity(const DataSet *data, int numOfSigmas, int minBmuHits) const
{
	bool success = true;
	float maxValue = -99999999;
	int maxValueDataSetRow = 0;
	int last = 0;
	
	auto binaryEigen = data->getBinaryEigen();
	auto contEigen = data->getContinuousEigen();
	
	for( size_t i = 0; i < data->size() + 1; i++ )
	{
		if( i == data->size() )
		{
			i = maxValueDataSetRow;
			last = 1;
		}
		
		// Extract sample vector
		Eigen::VectorXf v = data->getData(i);
	
		// Find best matching unit (bmu)
		SomIndex bmu = findRestrictedBmu(v, data->getValidityEigen(i).cast<float>(), minBmuHits, data->getWeightsEigen());
		
		int pos = width*bmu.getY() + bmu.getX();

		// Modified sigma vector
		Eigen::VectorXf sM = (sigmaMap[pos].array() > 0.00001f).select(0.00001f, sigmaMap[pos]);
		
		// Make map binary in order to get a proper interval (min - max) when comparing binaries
		Eigen::ArrayXf binaryMap = (binaryEigen.array() > 0).select(0.0f, map[pos]);
		binaryMap = (binaryEigen.array() > 0 && this->getNeuron(bmu).array() >= 0.5)
					.select(Eigen::VectorXf::Constant(binaryMap.size(), 1.0f), 0.0f);
		
		
		auto validEigen = data->getValidityEigen(i);
		
		
		// Calculate delta vector
		//Eigen::VectorXf delta = (validEigen.array()*((contEigen.array()*(v - this->getNeuron(bmu)).array()/sM.array()/numOfSigmas) + 
		//(binaryEigen.array()*(v.array()*(this->getNeuron(bmu).array().log())+(ones-v.array())*((ones-this->getNeuron(bmu).array()).log())))) ).matrix();
		
		Eigen::ArrayXf delta = /*contEigen.array()*/(v - this->getNeuron(bmu)).array()/sM.array()/numOfSigmas;
		//Eigen::ArrayXf binOnLoss = binaryEigen.array()*( v.array()*(this->getNeuron(bmu).array().log()) );
		//Eigen::ArrayXf binOffLoss = binaryEigen.array()*( (v-this->getNeuron(bmu)).array()*((ones - this->getNeuron(bmu).array()).log()) );
		
		Eigen::ArrayXf min = ((/*contEigen.array()*/((this->getNeuron(bmu) - sM*numOfSigmas).array())));// + 
		//(binaryEigen.array()*(binaryMap - 0.0001*ones)));
		Eigen::ArrayXf max = ((/*contEigen.array()*/((this->getNeuron(bmu) + sM*numOfSigmas).array())));// + 
		//(binaryEigen.array()*(binaryMap + 0.0001*ones)));
		
		if( _verbose )
		{
			std::cout << "X: " << bmu.getX() << "\tY: " << bmu.getY() << "\tID: " << i << "\n";
		}
		
		//std::cout << "# BMU hits: " << bmuHits[pos] << "\tX: " << bmu.getX() << "\tY: " << bmu.getY() << "\n";
		//if( ( ( (data->size()) > 100 && (i % (int)((data->size())/100)) == 0 ) || (data->size()) < 100 ) && _verbose )
		//	std::cout << "\rMeasuring dataset:" << 100*i/(data->size()) << "%";
			
		if(last && _verbose)
		{
			//std::cout << "\rMeasuring dataset:100%\n";
			std::cout << "Measurement nr: " << i << " is chosen\n";
		}
		
		for( int n = 0; n < v.size(); n++ )
		{
			if( delta(n) > maxValue )
			{
				maxValue = fabs(delta(n));
				maxValueDataSetRow = i;
			}
			
			// Replace NaN and Inf with zero (0)
			delta(n) = std::isnan(delta(n)) || std::isinf(delta(n)) ? 0 : delta(n);
			//binOnLoss(n) = std::isnan(binOnLoss(n)) ? 0 : binOnLoss(n);
			//binOffLoss(n) = std::isnan(binOffLoss(n)) ? 0 : binOffLoss(n);
			
			if( ((v(n) < min(n) || v(n) > max(n)) && validEigen(n)) && _verbose ) 
			{
				std::cout << data->getName(n) << "\t( " << v(n) << " ) [" << min(n) << " - " << max(n) << "]\n";
				success = false;
			}
			
			if(!_verbose && validEigen(n) && last)
			{
				// TODO: Make this some kind of return instead
				// std::cout << data->getName(n) << "\t" << delta(n)/*+binOnLoss(n)+binOffLoss(n)*/ << "\t" << v(n) << "\t" << min(n) << "\t" << max(n) << "\n";
				
				if( ((v(n) < min(n) || v(n) > max(n)) && validEigen(n)) )
					success = false;
			}
		}
		
		
		if(last)
			break;
		
		// std::cout << "\n";
	}
	
	return success;
}

// Trains on a single data point, i.e. a single record vector v is used to update 
// the map weights with parameters valid, weights, eta, sigma and weightDecayFunction
// lastBMU is the index of the BMU that this specific record found last epoch. 
// It's used as a starting point if we're looking for a local BMU.
SomIndex Som::trainSingle(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, const Eigen::VectorXf &weights, 
	const double eta, const double sigma, size_t &lastBMU, const int weightDecayFunction)
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
	size_t startX = std::max( (int)(bmu.getX()-2.5*sigma), 0);
	size_t startY = std::max( (int)(bmu.getY()-2.5*sigma), 0);
	
	size_t endX = std::min( (int)(bmu.getX()+2.5*sigma), (int)width);
	size_t endY = std::min( (int)(bmu.getY()+2.5*sigma), (int)height);
	
	for(size_t j = startY; j < endY; j++)
	{
		for(size_t i = startX; i < endX; i++)
		{
			// Calculate delta vector
			auto delta = ((v - map[j*width + i]).array()*valid.array()).matrix();
			
			// Strength of neighbourhood. Calculated for each neuron
			auto neighbourhoodWeight = calculateNeighbourhoodWeight(i, j, bmu.getX(), bmu.getY(), sigma);
			
			/*	Ref.
			* Incremental calculation of weighed mean and variance. Tony Finch. University of Cambrige Computing Service. February 2009.
			* */
			
			// Exponential weight decay function
			if(weightDecayFunction == 0)
			{
				weightMap[j*width + i] += neighbourhoodWeight*eta;
				map[j*width + i] += neighbourhoodWeight*eta*delta;
			}
			// Inverse proportional weight decay function
			else
			{
				weightMap[j*width + i] += neighbourhoodWeight; // Equation 47
				
				// Protection against division by zero
				auto tempWeight = weightMap[j*width + i] == 0 ? 1.0 : neighbourhoodWeight/weightMap[j*width + i];
				
				map[j*width + i] += tempWeight*(v - map[j*width + i]); // Equation 53
			}
			
			// Protection against division by zero
			auto tempWeight = weightMap[j*width + i] == 0 ? 0.000001 : weightMap[j*width + i];
			
			SMap[j*width + i] += neighbourhoodWeight*(delta.array()*(v - map[j*width + i]).array()).matrix(); // Equation 68
			sigmaMap[j*width + i] = (SMap[j*width + i].array()/tempWeight).abs().sqrt().matrix(); // Equation 69
			
		}
	}
	// Return best matching unit
	return bmu;
}

double Som::calculateNeighbourhoodWeight(
	const size_t &currentX, const size_t &currentY, 
	const size_t &bmuX, const size_t &bmuY, 
	const double &currentSigma)
{
	// Calculate neighbourhood fundction
	if (currentSigma > 1)
	{
		return std::exp(-((double)(currentX - (double)bmuX) * (double)(currentX - (double)bmuX) / (double)2 / currentSigma / currentSigma +
						  (double)(currentY - (double)bmuY) * (double)(currentY - (double)bmuY) / (double)2 / currentSigma / currentSigma));
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
	
	for(size_t i = 0; i < width*height; i++)
	{
		// Initialize each element in each neuron to a value between -sigma and sigma
		for(int n = 0; n < map[i].size(); n++)
		{
			map[i](n) = (double)((std::rand() % (int)(2000*sigma)) - (1000*sigma))/1000;
		}
		//std::cout << "InitV[" << i << "]: " << map[i] << "\n";
	}
}

void Som::updateUMatrix(const DataSet *data)
{
	std::vector<double> U(width*height);
	auto val = Eigen::VectorXf::Constant(map[0].size(), 1);
	auto weights = data->getWeightsEigen();
	double diagonalFactor = 0.3;
	double min = 10000000, 
			max = 0;
			//span;
			
	
	// For each vector in the map, calculate mean distance to closest neighbouring vectors.
	// Put these mean distances in the U matrix
	for(size_t i = 0; i < height; ++i)
	{
		for(size_t j = 0; j < width; ++j)
		{
			if( j > 0 && i > 0 && j < (width - 1) && i < (height - 1) ) // Middle part of map
			{
				U[i*width+j] = (this->euclidianWeightedDist(i*width+j, map[(i+0)*width+j-1], val, weights) + 
								this->euclidianWeightedDist(i*width+j, map[(i+0)*width+j+1], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j-0], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j-0], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j-1], val, weights)*diagonalFactor +
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j-1], val, weights)*diagonalFactor +
								this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j+1], val, weights)*diagonalFactor +
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j+1], val, weights)*diagonalFactor )/8;
			}
			else if( i == 0 && j > 0 && j < (width - 1) ) // Upper edge
			{
				U[i*width+j] = (this->euclidianWeightedDist(i*width+j, map[(i+0)*width+j-1], val, weights) + 
								this->euclidianWeightedDist(i*width+j, map[(i+0)*width+j+1], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j-0], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j-1], val, weights)*diagonalFactor +
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j+1], val, weights)*diagonalFactor )/5;
			}
			else if( i == (height - 1) && j > 0 && j < (width - 1) ) // Lower edge
			{
				U[i*width+j] = (this->euclidianWeightedDist(i*width+j, map[(i+0)*width+j-1], val, weights) + 
									this->euclidianWeightedDist(i*width+j, map[(i+0)*width+j+1], val, weights) +
									this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j-0], val, weights) +
									this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j-1], val, weights)*diagonalFactor +
									this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j+1], val, weights)*diagonalFactor )/5;
			}
			else if( j == 0 && i > 0 && i < (height - 1) ) // Left edge
			{
				U[i*width+j] = (this->euclidianWeightedDist(i*width+j, map[(i+0)*width+j+1], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j-0], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j-0], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j+1], val, weights)*diagonalFactor +
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j+1], val, weights)*diagonalFactor )/5;
			}
			else if( j == (width - 1) && i > 0 && i < (height - 1) ) // Right edge
			{
				U[i*width+j] = (this->euclidianWeightedDist(i*width+j, map[(i+0)*width+j-1], val, weights) + 
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j-0], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j-0], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j-1], val, weights)*diagonalFactor +
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j-1], val, weights)*diagonalFactor )/5;
			}
			else if( j == 0 && i == 0 ) // Top left corner
			{
				U[i*width+j] = (this->euclidianWeightedDist(i*width+j, map[(i+0)*width+j+1], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j-0], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j+1], val, weights)*diagonalFactor )/3;
			}
			else if( j == (width - 1) && i == 0 ) // Top rigth corner
			{
				U[i*width+j] = (this->euclidianWeightedDist(i*width+j, map[(i+0)*width+j-1], val, weights) + 
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j-0], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i+1)*width+j-1], val, weights)*diagonalFactor )/3;
			}
			else if( j == 0 && i == (height - 1) ) // Bottom left corner
			{
				U[i*width+j] = (this->euclidianWeightedDist(i*width+j, map[(i+0)*width+j+1], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j-0], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j+1], val, weights)*diagonalFactor )/3;
			}
			else if( j == (width - 1) && i == (height - 1) ) // Bottom right corner
			{
				U[i*width+j] = (this->euclidianWeightedDist(i*width+j, map[(i+0)*width+j-1], val, weights) + 
								this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j-0], val, weights) +
								this->euclidianWeightedDist(i*width+j, map[(i-1)*width+j-1], val, weights)*diagonalFactor )/3;
			}
			else
			{
				U[i*width+j] = 0;
			}
			
			// Find upper and lower bounds in U-matrix
			if( U[i*width+j] > max )
				max = U[i*width+j];
			if( U[i*width+j] < min )
				min = U[i*width+j];
		}
	}
	
	// Normalize to 0-255 interval
	//span = max - min;
	//if(_verbose) std::cout << "Span: " << span << "\nMin: " << min << "\nMax: " << max << "\n";
	for( size_t i = 0; i < uMatrix.size(); ++i)
	{
		uMatrix[i] = U[i];
	}
}

// Loops through the data set for numberOfEpochs turns while changing eta and sigma accordingly and updates map weights on each turn by calling trainSingle
void Som::train(DataSet *data, size_t numberOfEpochs, double eta0, double etaDecay, double sigma0, double sigmaDecay, int weightDecayFunction)
{	
	for(size_t i = 0; i < numberOfEpochs; ++i)
	{
		auto eta = eta0*std::exp(-etaDecay*i);
		auto sigma = sigma0*std::exp(-sigmaDecay*i);
		
		if(sigma < 1.0) sigma = 1.0;
		
		std::cout << "Epoch: " << i+1 << "/" << numberOfEpochs << "\teta: " << eta << "\tsigma: " << sigma << "\n";

		auto weights = data->getWeightsEigen();
		
		for(size_t j = 0; j < data->size(); ++j)
		{
			auto pos = trainSingle(data->getData(j), data->getValidityEigen(j).cast<float>(), weights, eta, sigma, data->getLastBMU(j), weightDecayFunction);
			
			addBmu(pos);
			if( ( data->size() > 100 && (j % (int)(data->size()/100)) == 0 ) || data->size() < 100 )
				std::cout << "\rTraining SOM:" << 100*j/data->size() << "%";
		}
		std::cout << "\rTraining SOM:100%\n";
	}
}

void Som::addBmu(SomIndex pos)
{
	bmuHits[this->getIndex(pos)] += 1;
}

void Som::displayUMatrix() const
{
	std::cout << "\nU-matrix:\n[";
	for(unsigned int i = 0; i < height; i++)
	{
		for(unsigned int j = 0; j < width; j++)
		{
			std::cout << (double)uMatrix[i*width+j] << " ";
		}
		std::cout << "; ";
	}
	std::cout << "]\n";
}

// The SOM is saved in an octave/matlab friendly file that holds all variables needed to analyze what's been learned in octave/matlab
void Som::save(const char *fileName) const
{
	FILE *fp;
	
	if( (fp = fopen(fileName, "w")) == NULL )
	{
		std::cout << "Could not open file " << fileName << " for writing. Quitting...\n";
		exit(EXIT_FAILURE);
	}
	
	fprintf(fp, "# This file is generated by Som.exe\n# It contains all data needed to evaluate a sample according to the SOM trained by Som.exe\n");
	fprintf(fp, "# name: som\n");
	fprintf(fp, "# type: matrix\n");
	fprintf(fp, "# ndims: 3\n");
	fprintf(fp, " %u %u %lu\n", this->getWidth(), this->getHeight(), (long unsigned int)map[0].rows());
	//fprintf(fp, "Vector length:%lu\n", (long unsigned int)map[0].rows());
	//fprintf(fp, "Som width:%u\nSom height:%u\n", this->getWidth(), this->getHeight());
	//fprintf(fp, "Som data:\n");
	
	// Print map data
	for(int i = 0; i<map[0].size(); i++)
	{
		for(unsigned int j = 0; j<height*width; j++)
		{
			fprintf(fp, " %f\n", map[j](i));
		}
		//fprintf(fp, "\n");
	}
	
	fprintf(fp, "\n\n# name: sigmaSom\n");
	fprintf(fp, "# type: matrix\n");
	fprintf(fp, "# ndims: 3\n");
	fprintf(fp, " %u %u %lu\n", this->getWidth(), this->getHeight(), (long unsigned int)sigmaMap[0].size());
	//fprintf(fp, "Sigma SOM data:\n");
	// Print map data
	for(int i = 0; i<sigmaMap[0].size(); i++)
	{
		for(unsigned int j = 0; j<height*width; j++)
		{
			fprintf(fp, " %f\n", sigmaMap[j](i));
		}
		//fprintf(fp, "\n");
	}
	
	fprintf(fp, "\n\n# name: weightMap\n");
	fprintf(fp, "# type: matrix\n");
	fprintf(fp, "# rows: %u\n", this->getWidth());
	fprintf(fp, "# columns: %u\n", this->getHeight());
	//fprintf(fp, "Weight som data:\n");
	// Print bmuHits data
	for(unsigned int j = 0; j<width; j++)
	{
		for(unsigned int i = 0; i<height; i++)
			fprintf(fp, "%f\t", weightMap[j*height + i]);
		
		fprintf(fp, "\n");
	}
	
	fprintf(fp, "\n\n# name: bmuHits\n");
	fprintf(fp, "# type: matrix\n");
	fprintf(fp, "# rows: %u\n", this->getWidth());
	fprintf(fp, "# columns: %u\n", this->getHeight());
	//fprintf(fp, "\nBmuHits data:\n");
	// Print bmuHits data
	for(unsigned int j = 0; j<width; j++)
	{
		for(unsigned int i = 0; i<height; i++)
			fprintf(fp, "%d\t", bmuHits[j*height + i]);
		fprintf(fp, "\n");
	}
	
	fprintf(fp, "\n\n# name: U\n");
	fprintf(fp, "# type: matrix\n");
	fprintf(fp, "# rows: %u\n", this->getWidth());
	fprintf(fp, "# columns: %u\n", this->getHeight());
	//fprintf(fp, "\nU-matrix data:\n");
	// Print U-matrix data
	for(unsigned int j = 0; j<width; j++)
	{
		for(unsigned int i = 0; i<height; i++)
			fprintf(fp, "%f\t", uMatrix[j*height + i]);
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
	
	while ( std::getline (file,line) )
	{
		// Check vector length
		if( line.compare(0, 9, "# ndims: ") == 0 )
		{
			std::getline(file, line);
			found = line.find_last_of(" ");
			loadedVectorLength = std::stoul(line.substr(found+1), &ptr, 10);
			initV.resize(loadedVectorLength);
			//std::cout << "Found vector length:" << loadedVectorLength << "\n";
		}
		
		// Find width of som
		else if( line.compare(0, 8, "# rows: ") == 0 )
		{
			height = std::stoul(&(line[8]), &ptr, 10);
			continue;
		}
		
		// Find height of som
		else if( line.compare(0,11, "# columns: ") == 0 )
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
	
	while ( !file.eof() )
	{
		
		getline(file, line);
		//std::cout << "|" << line << "|\n";
		
		if( line[0] == '#' )
		{
		//std::cout << line << "\n";
			// Check vector name
			if( line.compare(0, 8, "# name: ") == 0 )
			{
				//std::cout << "Found name line: ";
				if( line.compare(8, 20, "som") == 0 )
				{
					readSomData = 1;
					readSigmaSomData = 0;
					readWeightMapData = 0;
					readBmuHitsData = 0;
					readUMatrixData = 0;
					//std::cout << "SOM\n";
					continue;
				}
				else if( line.compare(8, 20, "sigmaSom") == 0 )
				{
					readSomData = 0;
					readSigmaSomData = 1;
					readWeightMapData = 0;
					readBmuHitsData = 0;
					readUMatrixData = 0;
					//std::cout << "sigmaSom\n";
					continue;
				}
				else if( line.compare(8, 20, "weightMap") == 0 )
				{
					readSomData = 0;
					readSigmaSomData = 0;
					readWeightMapData = 1;
					readBmuHitsData = 0;
					readUMatrixData = 0;
					//std::cout << "weightMap\n";
					continue;
				}
				else if( line.compare(8, 20, "bmuHits") == 0 )
				{
					readSomData = 0;
					readSigmaSomData = 0;
					readWeightMapData = 0;
					readBmuHitsData = 1;
					readUMatrixData = 0;
					//std::cout << "bmuHits\n";
					continue;
				}
				else if( line.compare(8, 20, "U") == 0 )
				{
					readSomData = 0;
					readSigmaSomData = 0;
					readWeightMapData = 0;
					readBmuHitsData = 0;
					readUMatrixData = 1;
					//std::cout << "U\n";
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
			else if( line.compare(0, 8, "# type: ") == 0 )
			{
				//std::cout << "Found type line: ";
				if( line.compare(8, 20, "matrix" ) != 0 )
				{
					std::cout << "Incorrect type in SOM file!\n";
					exit(EXIT_FAILURE);
				}
				//std::cout << "matrix\n";
			}
		
			// Find vector dimension
			else if( line.compare(0, 9, "# ndims: ") == 0 )
			{
				// dimension = std::stoul(&(line[9]), &ptr, 10);
				//std::cout << "Found ndim line: " << dimension << "\n";
				
				getline(file, line);
				
				//std::cout << "|" << line << "\n";
				
				// To go: Do something with dimension length
				
				continue;
			}
			
			// Find 2D matrix height
			else if( line.compare(0, 8, "# rows: ") == 0 )
			{
				width = std::stoul(&(line[8]), &ptr, 10);
				//std::cout << "Found rows line: " << height << "\n";
				continue;
			}
			
			// Find 2D matrix width
			else if( line.compare(0, 11, "# columns: ") == 0 )
			{
				height = std::stoul(&(line[11]), &ptr, 10);
				//std::cout << "Found columns line: " << width << "\n";
				continue;
			}
		}
		
		
		// Save SOM data to Eigen vector
		else if( readSomData && (std::isdigit(line[1]) || line[1] == '-') )
		{
			map[d](i) = std::stod(line, &ptr);
			
			//if( i == 0 )
			//	std::cout << "Found som value:" << map[d](i) << " for element:" << i << "\tat pos:" << d << "\n";
			
			d++;
			
			if( d >= height*width )
			{
				d = 0;
				i++;
				if( i >= map[0].size() )
				{
					i = 0;
					readSomData = 0;
				}
			}
			
			continue;
		}
		
		// Save sigmaSOM data to Eigen vector
		else if( readSigmaSomData && (std::isdigit(line[1]) || line[1] == '-') )
		{
			sigmaMap[d](i) = std::stod(line, &ptr);
			//std::cout << "Found sigmaSom value:" << sigmaMap[i](d) << " for element:" << d << "\n";
			d++;
			
			if( d >= height*width )
			{
				d = 0;
				i++;
				if( i >= sigmaMap[0].size() )
				{
					i = 0;
					readSigmaSomData = 0;
				} 
			}
			
			continue;
		}
		
		// Save  weightMap data to Eigen vector
		else if( readWeightMapData && std::isdigit(line.c_str()[0]) )
		{
			for(unsigned int l = 0; l < width; l++)
			{
				columnMarker = (char*)line.c_str();
				
				for(unsigned int k = 0; k < height; k++)
				{
					//printf("|%s\n", columnMarker);
					weightMap(l*height + k) = std::strtod(columnMarker, &columnMarker);
					//columnMarker = (char*)ptr;
					//std::cout << weightMap(l*width + k) << "\t";
				}
				//std::cout << "\n";
				getline(file, line);
			}
			
			continue;
		}
		
		// Save bmu hits data to vector
		else if( readBmuHitsData && std::isdigit(line.c_str()[0]) )
		{
			for(unsigned int l = 0; l < width; l++)
			{
				columnMarker = (char*)line.c_str();
				
				for(unsigned int k = 0; k < height; k++)
				{
					//printf("|%s\n", columnMarker);
					bmuHits[l*height + k] = std::strtoul(columnMarker, &columnMarker, 10);
					//columnMarker = (char*)ptr;
					//std::cout << bmuHits[l*height + k] << "\t";
				}
				//std::cout << "\n";
				getline(file, line);
			}
			
			continue;
		}
		
		// Save U-matrix data to vector
		else if( readUMatrixData && std::isdigit(line.c_str()[0]) )
		{
			for(unsigned int l = 0; l < width; l++)
			{
				columnMarker = (char*)line.c_str();
				
				for(unsigned int k = 0; k < height; k++)
				{
					//printf("|%s\n", columnMarker);
					uMatrix[l*height + k] = std::strtod(columnMarker, &columnMarker);
					//columnMarker = (char*)ptr;
					//std::cout << uMatrix[l*width + k] << "\t";
				}
				//std::cout << "\n";
				getline(file, line);
			}
			
			continue;
		}
		else
		{
			;
		}
		
		
		//std::cout << line << '\n';
	}
	
	file.close();
	
}