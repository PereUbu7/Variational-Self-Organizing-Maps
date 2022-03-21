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
		unsigned int height, width;
		bool _verbose;
	
	// Lägg till att ta hänsyn till validitet i functionerna som nu tar del av den variabeln
	public:
		Som(size_t, size_t, size_t, bool verbose = false);
		Som(const char*, bool verbose = false);
		~Som() = default;
		void train(DataSet*, int, double, double, double, double, int);
		double evaluate(const DataSet*);
		SomIndex trainSingle(Eigen::VectorXf, std::vector<int>, std::vector<double>, double, double, int*, int);
		int measureSimilarity(const DataSet*, int, int);
		int autoEncoder(const DataSet*, int);
		size_t variationalAutoEncoder(const DataSet*, int);
		SomIndex findBmu(Eigen::VectorXf, std::vector<int>, std::vector<double>) const;
		SomIndex findLocalBmu(Eigen::VectorXf, std::vector<int>, int, std::vector<double>) const;
		SomIndex findRestrictedBmu(Eigen::VectorXf, std::vector<int>, int, std::vector<double>) const;
		std::vector<double> findRestrictedBmd(Eigen::VectorXf, std::vector<int>, int, std::vector<double>) const;
		double euclidianWeightedDist(SomIndex, Eigen::VectorXf, std::vector<int>, std::vector<double>) const;
		double euclidianWeightedDist(int pos, Eigen::VectorXf, std::vector<int>, std::vector<double>) const;
		void display() const;
		void displayUMatrix();
		unsigned int getHeight() const;
		unsigned int getWidth() const;
		unsigned int getIndex(SomIndex) const;
		Eigen::VectorXf getNeuron(SomIndex) const;
		Eigen::VectorXf getNeuron(int) const;
		Eigen::VectorXf getSigmaNeuron(SomIndex) const;
		Eigen::VectorXf getSigmaNeuron(int) const;
		void randomInitialize(int, float);
		void addBmu(SomIndex);
		void updateUMatrix(const DataSet *data);
		void save(const char*);
		void load(const char*);
		Eigen::VectorXf getSizeFromFile(const char*);
};
