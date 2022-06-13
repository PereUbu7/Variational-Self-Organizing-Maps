#include "MnistDataLoader.hpp"

#include "Eigen/Dense"

#include <vector>
#include <cstdint>

bool MnistDataLoader::open(const char *path)
{
    _filePath = path;

    return true;
}

size_t MnistDataLoader::load()
{
    auto dataset = mnist::read_dataset<std::vector, std::vector, uint8_t, uint8_t>(_filePath, m_currentIndex, m_maxLoadCount.value_or(0), m_maxLoadCount.value_or(0));

    auto numberOfSamples = dataset.training_labels.size();

    if(numberOfSamples == 0) m_currentIndex = 0;
    else if(numberOfSamples > 0) m_currentIndex += numberOfSamples;

    data = std::vector<RowData>();
    data.reserve(numberOfSamples);

    for(size_t index{0}; index < numberOfSamples; ++index)
    {
        Eigen::VectorXf imageValues = Eigen::Map<Eigen::Matrix<uint8_t, Eigen::Dynamic, 1>>(
            dataset.training_images.at(index).data(), 
            dataset.training_images.at(index).size()).cast<float>();
        data.emplace_back(
            imageValues,
            std::vector<int>(28*28, 1) 
        );
    }

    return numberOfSamples;
}

float &MnistDataLoader::getWeight(size_t index)
{
    return _weights.at(index);
}

const std::vector<float> &MnistDataLoader::getWeights() const noexcept
{
    return _weights;
}
const std::vector<int> &MnistDataLoader::getBinary() const noexcept
{
    return _isBinary;
}
const std::vector<int> &MnistDataLoader::getContinuous() const noexcept
{
    return _isContinuous;
}
std::string MnistDataLoader::getName(size_t index) const noexcept
{
    return _names.at(index);
}
const std::vector<std::string> &MnistDataLoader::getNames() const noexcept
{
    return _names;
}
size_t MnistDataLoader::getDepth() const noexcept
{
    return _names.size();
}

void MnistDataLoader::generateNames()
{
    _names.reserve(28 * 28);

    for (size_t index{0}; index < 28 * 28; ++index)
    {
        auto x_value = index % 28;
        auto y_value = index / 28;

        std::ostringstream name;
        name << std::to_string(x_value) << "x" << std::to_string(y_value);
        _names.emplace_back(name.str());
    }
}