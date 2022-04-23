#pragma once

#define _GLIBCXX_USE_C99 1

#include "SomIndex.hpp"
#include "DataSet.hpp"
#include "UMatrix.hpp"

#include <vector>
#include <atomic>
#include <mutex>
#include "Eigen/Dense"

#define VERSION 1.00

#define ARG_VERBOSE 1
#define ARG_SETTING 2

#define ARG_DB_FILE 3
#define ARG_SOM_FILE 4

// Training parameters
#define ARG_SOM_HEIGHT 5
#define ARG_SOM_WIDTH 6
#define ARG_SOM_ETA0 7
#define ARG_SOM_ETA_DEC 8
#define ARG_SOM_SIGMA0 9
#define ARG_SOM_SIGMA_DEC 10
#define ARG_SOM_EPOCHS 11
#define ARG_SOM_INIT_SIGMA 12
#define ARG_SOM_WEIGHT_DECAY_FUNCTION 13

// Measuring parameters
#define ARG_ALLOWED_STD_DEV 6
#define ARG_MIN_BMU_HITS 5

#define SIGMA_SWITCH_TO_LOCAL 1

class Som
{
protected:
	struct TrainingReturnValue
	{
		SomIndex bmu;
		Eigen::VectorXf residual;
	};
	struct Metrics
	{
		std::vector<float> MeanSquaredError;
		Metrics() : MeanSquaredError{} {}
		Metrics(size_t size) : MeanSquaredError(size) {}
	};
	std::vector<Eigen::VectorXf> map;
	std::vector<Eigen::VectorXf> sigmaMap;
	std::vector<Eigen::VectorXf> SMap;
	Metrics metrics;
	Eigen::VectorXf weightMap;
	std::vector<int> bmuHits;
	std::vector<double> uMatrix;
	std::atomic<bool> _isTraining;
	size_t height, width, depth;
	bool _verbose;

	// Lägg till att ta hänsyn till validitet i functionerna som nu tar del av den variabeln
public:
	enum class WeigthDecayFunction
	{
		Exponential,
		InverseProportional
	};
	std::mutex metricsMutex;

	Som(size_t width, size_t height, size_t depth, bool verbose = false);
	Som(const char *filename, bool verbose = false);
	Som(const Som &som) : map{som.map},
						  sigmaMap{som.sigmaMap},
						  SMap{som.SMap},
						  metrics{},
						  weightMap{som.weightMap},
						  bmuHits{som.bmuHits},
						  uMatrix{som.uMatrix},
						  _isTraining{},
						  height{som.height},
						  width{som.width},
						  depth{som.depth},
						  _verbose{som._verbose},
						  metricsMutex{}
	{
		_isTraining.store(som._isTraining.load(std::memory_order_seq_cst), std::memory_order_seq_cst);
	}
	~Som() = default;
	void train(DataSet &data, size_t numberOfEpochs, double eta0, double etaDecay, double sigma0, double sigmaDecay, WeigthDecayFunction weightDecayFunction, bool updateUMatrixAfterEpoch = false);
	double evaluate(const DataSet &dataset) const;
	TrainingReturnValue trainSingle(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, const Eigen::VectorXf &weights, const double eta, const double sigma, size_t &lastBMU, const WeigthDecayFunction weightDecayFunction);
	int measureSimilarity(const DataSet *dataset, int numberOfSigmas, int minBmuHits) const;
	int autoEncoder(const DataSet *dataset, int minBmuHits) const;
	size_t variationalAutoEncoder(const DataSet *dataset, int minBmuHits) const;
	SomIndex findBmu(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, const Eigen::VectorXf &weights) const;
	SomIndex findLocalBmu(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, const size_t &lastBMUref, const Eigen::VectorXf &weights) const;
	SomIndex findRestrictedBmu(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, const int minBmuHits, const Eigen::VectorXf &weights) const;
	std::vector<double> findRestrictedBmd(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, int minBmuHits, const Eigen::VectorXf &weights) const;
	double euclidianWeightedDist(
		const SomIndex &pos,
		const Eigen::VectorXf &v,
		const Eigen::VectorXf &valid,
		const Eigen::VectorXf &weights) const;
	double euclidianWeightedDist(
		const size_t &pos,
		const Eigen::VectorXf &v,
		const Eigen::VectorXf &valid,
		const Eigen::VectorXf &weights) const;
	void display() const;
	void displayUMatrix() const;
	UMatrix getUMatrix() const noexcept;
	Eigen::VectorXf getWeigthMap() const noexcept;
	std::vector<int> getBmuHits() const noexcept;
	size_t getHeight() const noexcept;
	size_t getWidth() const noexcept;
	size_t getIndex(SomIndex index) const noexcept;
	Eigen::VectorXf getNeuron(SomIndex index) const noexcept;
	Eigen::VectorXf getNeuron(size_t index) const noexcept;
	Eigen::VectorXf getSigmaNeuron(SomIndex index) const noexcept;
	Eigen::VectorXf getSigmaNeuron(size_t index) const noexcept;
	float getMaxValueOfFeature(size_t modelVectorIndex) const;
	float getMinValueOfFeature(size_t modelVectorIndex) const;
	float getMaxSigmaOfFeature(size_t modelVectorIndex) const;
	float getMinSigmaOfFeature(size_t modelVectorIndex) const;
	Metrics getMetrics() const noexcept;
	bool isTraining() const noexcept;
	void randomInitialize(int seed, float sigma);
	void addBmu(SomIndex position);
	void updateUMatrix(const Eigen::VectorXf &weights);
	void save(const char *filename) const;
	void load(const char *filename);
	Eigen::VectorXf getSizeFromFile(const char *filename);

	double static calculateNeighbourhoodWeight(
		const size_t &currentX, const size_t &currentY,
		const size_t &bmuX, const size_t &bmuY,
		const double &currentSigma);

	Som &operator=(const Som &other)
	{
		map = other.map;
		sigmaMap = other.sigmaMap;
		SMap = other.SMap;
		weightMap = other.weightMap;
		bmuHits = other.bmuHits;
		uMatrix = other.uMatrix;
		_isTraining.store(other._isTraining.load(std::memory_order_seq_cst), std::memory_order_seq_cst);
		height = other.height;
		width = other.width;
		depth = other.depth;
		_verbose = other._verbose;

		return *this;
	}
};
