#ifndef SOM_HPP_INCLUDED
#define SOM_HPP_INCLUDED

#define _GLIBCXX_USE_C99 1

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <cctype>
#include <vector>
#include <algorithm>
#include <random>
#include "Eigen/Dense"
#include "sqlite/sqlite3.h"

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

#define TRUE 1
#define FALSE 0

extern int verbose;

int saveElement(void*, int, char**, char**);

class DataBase
{
	protected:
		sqlite3 *db;
	public:
		DataBase();
		~DataBase();
		int open(const char*);
		void close();
		int getElement(char*, int, std::string);
		int getMax(char*, std::string);
		int getRecord(char*[], int, std::string);
		int rows();
		int minId();
		int maxId();
		int doesExist(int);
		int startTransaction();
		int endTransaction();
};

class SomIndex
{
	protected:
		int x,y;
	
	public:
		SomIndex(int, int);
		~SomIndex();
		//int getSomIndex(Som);
		unsigned int getX() const;
		unsigned int getY() const;
		void setX(int);
		void setY(int);
};

class DataSet
{
	protected:
		std::vector<Eigen::VectorXf> data;
		std::vector<std::vector<int>> valid;
		std::vector<unsigned int> index;
		std::vector<std::string> variableNames;
		std::vector<int> binary;
		std::vector<int> continuous;
		std::vector<double> weight;
		std::vector<int> lastBMU;
		unsigned int depth;
		int n;
	
public:
		DataSet(const char*);
		~DataSet();
		Eigen::VectorXf getData(int) const;
		std::vector<int> getValidity(int) const;
		std::vector<int> getBinary() const;
		std::vector<int> getContinuous() const;
		std::vector<double> getWeights() const;
		int *getLastBMU(int);
		unsigned int size() const;
		void addVector(Eigen::VectorXf);
		void addSet(DataSet d);
		void loadTextFile(const char*);
		void display() const;
		void loadDataBase(DataBase *db);
		void loadIcanFilter(const char *fileName);
		unsigned int vectorLength() const;
		
		std::string getName(int) const;
		
		void shuffle();
};

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
	
	// Lägg till att ta hänsyn till validitet i functionerna som nu tar del av den variabeln
	public:
		Som(unsigned int, unsigned int, unsigned int);
		Som(const char*);
		~Som();
		void train(DataSet*, int, double, double, double, double, int);
		double evaluate(const DataSet*);
		SomIndex trainSingle(Eigen::VectorXf, std::vector<int>, std::vector<double>, double, double, int*, int);
		int measureSimilarity(const DataSet*, int, int);
		int autoEncoder(const DataSet*, int);
		int variationalAutoEncoder(const DataSet*, int);
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


#endif
