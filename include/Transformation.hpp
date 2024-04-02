#pragma once

#include <functional>

#include <string>

template<typename T>
struct Transformation
{
    /* Standard transform by default TODO: Use dispersion and weight */
    std::function<T(const T& value, const T& model, const T& dispersion, const T& valueWeight)> Comparer = 
        [](const T &value, const T &model, const T &dispersion, const T &valueWeight) { return model - value; };

    std::function<T(const T& value, const T& model, const T& valueWeight)> Stepper = 
        [](const T &value, const T &model, const T& valueWeight) { return value - model;  };

    std::function<T(const T& model)> Displayer = 
        [](const T &model) { return model; };

    std::string Name = std::string{"Standard transformation"};
};