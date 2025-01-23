#pragma once

#include "Eigen/Dense"
#include <functional>

#include <string>
#include <sstream>


struct Transformation
{
    using T = Eigen::VectorXf;
    /* Standard transform by default TODO: Use dispersion and weight */
    std::function<T(const T &value, const T &model, const T &dispersion, const T &valueWeight)> Comparer {
        [](const T &value, const T &model, const T &dispersion, const T &valueWeight)
    { return model - value; }
    };

    std::function<T(const T &value, const T &model, const T &valueWeight)> Stepper {
        [](const T &value, const T &model, const T &valueWeight)
    { return value - model; }
    };

    std::vector<std::string> names{};
    std::function<std::vector<std::string>(const T &model)> Displayer{
        [&names = names](const T &model)
        {
            return names;
        }};

    std::string Name{"Standard transformation"};

    static Transformation Standard(const std::vector<std::string> &columnNames);
    static Transformation CombinatorialLinearRegression(const std::vector<std::string> &columnNames);
};