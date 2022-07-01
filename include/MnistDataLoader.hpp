#pragma once

#include "IDataLoader.hpp"

#include "mnistReader/mnist_reader.hpp"

#include <string>

class MnistDataLoader : public IDataLoader
{
protected:
    std::vector<float> _weights;
    std::vector<int> _isBinary;
    std::vector<int> _isContinuous;
    std::vector<std::string> _names;
    std::string _filePath;
    bool _verbose;

    void close();
    int getMax(char *, std::string);
    int getRecord(char *[], int, std::string) const;
    void generateNames();

public:
    MnistDataLoader(std::optional<size_t> maxLoadCount = std::nullopt, bool verbose = false) : 
        IDataLoader{maxLoadCount},
        _weights(28*28 + 10, 1.0f),
        _isBinary(28*28 + 10, 0),
        _isContinuous(28*28 + 10, 1),
        _names{},
        _filePath{}, 
        _verbose{verbose} 
    {
        generateNames();
    };
    // ~MnistDataLoader() = default override;

    size_t load() override;
    std::vector<RowData> getPreview(size_t count) override;
    bool open(const char *path) override;
    float &getWeight(size_t index) override;
    const std::vector<float> &getWeights() const noexcept override;
    const std::vector<int> &getBinary() const noexcept override;
    const std::vector<int> &getContinuous() const noexcept override;
    std::string getName(size_t index) const noexcept override;
    const std::vector<std::string> &getNames() const noexcept override;
    size_t getDepth() const noexcept override;
    bool isAtStartOfDataStream() const noexcept override 
	{ 
		return m_currentIndex == 0;
	}
};