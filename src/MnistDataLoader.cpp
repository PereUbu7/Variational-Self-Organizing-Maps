#include "MnistDataLoader.hpp"

#include "Eigen/Dense"

#include <vector>
#include <cstdint>
#include <algorithm>

bool MnistDataLoader::open(const char *path)
{
    _filePath = path;

    return true;
}

std::vector<RowData> MnistDataLoader::getPreview(size_t count)
{
    auto dataset = mnist::read_dataset<std::vector, std::vector, uint8_t, uint8_t>(_filePath, 0, count, count);
    auto previewData = std::vector<RowData>();
    previewData.reserve(dataset.training_images.size());

    std::transform(
        dataset.training_images.begin(),
        dataset.training_images.end(), 
        dataset.training_labels.begin(),
        std::back_inserter(previewData),
        [](auto& image, const auto& label)
        {
            auto oneHotEncodedLabel = std::vector<uint8_t>(10, 0);
            if(label >= 0 && label < 10) oneHotEncodedLabel[label] = 1;

            /* Concatenate image with one-hot-encoded label */
            image.reserve(image.size() + oneHotEncodedLabel.size());
            image.insert(std::end(image), std::begin(oneHotEncodedLabel), std::end(oneHotEncodedLabel));
            Eigen::VectorXf values = Eigen::Map<Eigen::Matrix<uint8_t, Eigen::Dynamic, 1>>(
                image.data(), 
                image.size()).cast<float>();
            
            return RowData {
                values,
                std::vector<int>(values.size(), 1)
            };
        });
    return previewData;
}

size_t MnistDataLoader::load()
{
    auto dataset = mnist::read_dataset<std::vector, std::vector, uint8_t, uint8_t>(_filePath, m_currentIndex, m_maxLoadCount.value_or(0), m_maxLoadCount.value_or(0));

    auto numberOfSamples = dataset.training_labels.size();

    if(numberOfSamples == 0) m_currentIndex = 0;
    else if(numberOfSamples >= 60000) m_currentIndex = 0;
    else if(numberOfSamples > 0) m_currentIndex += numberOfSamples;

    data = std::vector<RowData>();
    data.reserve(numberOfSamples);

    std::transform(
        dataset.training_images.begin(),
        dataset.training_images.end(), 
        dataset.training_labels.begin(),
        std::back_inserter(data),
        [](auto& image, const auto& label)
        {
            auto oneHotEncodedLabel = std::vector<uint8_t>(10, 0);
            if(label >= 0 && label < 10) oneHotEncodedLabel[label] = 1;

            /* Concatenate image with one-hot-encoded label */
            image.reserve(image.size() + oneHotEncodedLabel.size());
            image.insert(std::end(image), std::begin(oneHotEncodedLabel), std::end(oneHotEncodedLabel));
            Eigen::VectorXf values = Eigen::Map<Eigen::Matrix<uint8_t, Eigen::Dynamic, 1>>(
                image.data(), 
                image.size()).cast<float>();
            
            return RowData {
                values,
                std::vector<int>(values.size(), 1)
            };
        });

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
    _names.reserve(28 * 28 + 10);

    for (size_t index{0}; index < 28 * 28; ++index)
    {
        auto x_value = index % 28;
        auto y_value = index / 28;

        std::ostringstream name;
        name << std::to_string(x_value) << "x" << std::to_string(y_value);
        _names.emplace_back(name.str());
    }
    for (size_t index{0}; index < 10; ++index)
    {
        std::ostringstream name;
        name << "label:" << std::to_string(index);
        _names.emplace_back(name.str());
    }
}