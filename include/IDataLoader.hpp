#pragma once

#include "ColumnSpec.hpp"

#include <string>
#include <vector>
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
		IDataLoader(std::optional<size_t> maxLoadCount = std::nullopt) 
			: m_maxLoadCount{maxLoadCount}, m_currentIndex{0} {};
		IDataLoader(const IDataLoader&) = default;
		IDataLoader& operator=(const IDataLoader&) = default;
		IDataLoader(IDataLoader&&) = default;
		IDataLoader& operator=(IDataLoader&&) = default;
		virtual size_t load() = 0;
		virtual std::vector<RowData> getPreview(size_t count) = 0;
		virtual bool open(const char *path) = 0;
		virtual std::vector<std::string> findAllColumns() = 0;
		std::vector<RowData> data;

		virtual void setColumnSpec(const std::vector<ColumnSpec> columnSpec) noexcept = 0;
		virtual const std::vector<ColumnSpec> getColumnSpec() noexcept = 0;
		virtual float getWeight(size_t index) = 0;
		virtual const std::vector<float> getWeights() const noexcept = 0;
		virtual const std::vector<int> &getBinary() const noexcept = 0;
		virtual const std::vector<int> &getContinuous() const noexcept = 0;
		virtual std::string getName(size_t index) const noexcept = 0;
		virtual const std::vector<std::string> getNames() const noexcept = 0;
		virtual size_t getDepth() const noexcept = 0;
		virtual bool isAtStartOfDataStream() const noexcept = 0;
	
	protected:
		std::optional<size_t> m_maxLoadCount;
		size_t m_currentIndex;
};