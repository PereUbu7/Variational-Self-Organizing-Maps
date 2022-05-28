#pragma once

#include <string>
#include <optional>
#include "Eigen/Dense"


struct RowData
{
    Eigen::VectorXf values;
	std::vector<int> valid;
};


class IDataLoader
{
	public:
		virtual ~IDataLoader() = default;
		IDataLoader() = default;
		IDataLoader(const IDataLoader&) = default;
		IDataLoader& operator=(const IDataLoader&) = default;
		IDataLoader(IDataLoader&&) = default;
		IDataLoader& operator=(IDataLoader&&) = default;
		virtual size_t load(std::optional<size_t> maxCount = std::nullopt) = 0;
		virtual bool open(const char *path) = 0;
		std::vector<RowData> data;

		virtual float &getWeight(size_t index) = 0;
		virtual const std::vector<float> &getWeights() const noexcept = 0;
		virtual const std::vector<int> &getBinary() const noexcept = 0;
		virtual const std::vector<int> &getContinuous() const noexcept = 0;
		virtual std::string getName(size_t index) const noexcept = 0;
		virtual const std::vector<std::string> &getNames() const noexcept = 0;
		virtual size_t getDepth() const noexcept = 0;
};