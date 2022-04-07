#pragma once

#include "Database.hpp"

#include <vector>
#include <string>
#include "Eigen/Dense"

class DataSet
{
	protected:
		std::vector<Eigen::VectorXf> data;
		std::vector<std::vector<int>> valid;
		std::vector<size_t> index;
		std::vector<std::string> variableNames;
		std::vector<int> binary;
		std::vector<int> continuous;
		std::vector<float> weight;
		std::vector<size_t> lastBMU;
		size_t depth, n;
		bool _verbose;
	
public:
		DataSet(const char*, bool verbose = false);
		~DataSet() = default;
		Eigen::VectorXf getData(size_t) const;
		const Eigen::VectorXi getValidity(size_t index) const;
		const Eigen::ArrayXi getBinary() const;
		const Eigen::ArrayXi getContinuous() const;
		const Eigen::VectorXf getWeights() const;
		std::string getName(size_t) const;
		size_t &getLastBMU(size_t);
		size_t size() const;
		void addVector(Eigen::VectorXf);
		void addSet(const DataSet &d);
		void loadTextFile(const char*);
		void display() const;
		void loadDataBase(DataBase *db);
		void loadIcanFilter(const char *fileName);
		size_t vectorLength() const;
		
		
		void shuffle();
};