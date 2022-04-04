#pragma once

#define _GLIBCXX_USE_C99 1

#include "SomIndex.hpp"
#include "DataSet.hpp"

#include <vector>
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
	std::vector<Eigen::VectorXf> map;
	std::vector<Eigen::VectorXf> sigmaMap;
	std::vector<Eigen::VectorXf> SMap;
	Eigen::VectorXf weightMap;
	std::vector<int> bmuHits;
	std::vector<double> uMatrix;
	size_t height, width;
	bool _verbose;

	// Lägg till att ta hänsyn till validitet i functionerna som nu tar del av den variabeln
public:
	Som(size_t, size_t, size_t, bool verbose = false);
	Som(const char *, bool verbose = false);
	~Som() = default;
	void train(DataSet *, size_t, double, double, double, double, int);
	double evaluate(const DataSet&) const;
	SomIndex trainSingle(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, const Eigen::VectorXf &weights, const double eta, const double sigma, size_t &lastBMU, const int weightDecayFunction);
	int measureSimilarity(const DataSet *, int, int) const;
	int autoEncoder(const DataSet *, int) const;
	size_t variationalAutoEncoder(const DataSet *, int) const;
	SomIndex findBmu(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, const Eigen::VectorXf &weights) const;
	SomIndex findLocalBmu(const Eigen::VectorXf &v, const Eigen::VectorXf &valid, const size_t &lastBMUref, const Eigen::VectorXf &weights) const;
	// TODO:
	SomIndex findRestrictedBmu(const Eigen::VectorXf&, const std::vector<int>&, const int, const std::vector<float>&) const;
	// TODO:
	std::vector<double> findRestrictedBmd(const Eigen::VectorXf&, const std::vector<int>&, int, const std::vector<float>&) const;
	double euclidianWeightedDist(
		const SomIndex &pos, 
		const Eigen::VectorXf &v, 
		const Eigen::VectorXf &valid, 
		const Eigen::VectorXf &weights) const;
	// TO be removed
	double euclidianWeightedDist(
		const size_t &pos, 
		const Eigen::VectorXf&, 
		const std::vector<int>&, 
		const std::vector<float>&) const;
	double euclidianWeightedDist(
		const size_t &pos, 
		const Eigen::VectorXf &v, 
		const Eigen::VectorXf &valid, 
		const Eigen::VectorXf &weights) const;
	void display() const;
	void displayUMatrix() const;
	unsigned int getHeight() const;
	unsigned int getWidth() const;
	unsigned int getIndex(SomIndex) const;
	Eigen::VectorXf getNeuron(SomIndex) const;
	Eigen::VectorXf getNeuron(size_t) const;
	Eigen::VectorXf getSigmaNeuron(SomIndex) const;
	Eigen::VectorXf getSigmaNeuron(size_t) const;
	void randomInitialize(int, float);
	void addBmu(SomIndex);
	void updateUMatrix(const DataSet *data);
	void save(const char *) const;
	void load(const char *);
	Eigen::VectorXf getSizeFromFile(const char *);

	double static calculateNeighbourhoodWeight(
		const size_t &currentX, const size_t &currentY,
		const size_t &bmuX, const size_t &bmuY,
		const double &currentSigma);
};
