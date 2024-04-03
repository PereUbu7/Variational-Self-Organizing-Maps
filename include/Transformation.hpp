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

    static Transformation<const Eigen::VectorXf> CombinatorialLinearRegression()
    {
        return Transformation<const Eigen::VectorXf>{
                .Comparer = [](const Eigen::VectorXf &value, const Eigen::VectorXf &model, const Eigen::VectorXf &dispersion, const Eigen::VectorXf &valueWeight)
                { 
                    /* The model vector is split into two parts, A and B, where A*x' + B = y'
                     * x' and y' together constitutes all combinations of the elements of 'value',
                     * not considering direction of the combination with itself (so only the upper half matrix) */
                Eigen::VectorXf A = model.head(model.size()/2);
                Eigen::VectorXf B = model.tail(model.size()/2);

                /* Create all combinations of value's element (n*(n-1)/2 combinations)*/
                Eigen::VectorXf xPrime(A.size());
                Eigen::VectorXf yPrime(A.size());

                int index = 0;
                for (int i = 0; i < value.size(); ++i) {
                    for (int j = i + 1; j < value.size(); ++j) {
                        xPrime[index] = value[i];
                        yPrime[index] = value[j];
                        index++;
                    }
                }

                /* Calculate residual */
                Eigen::VectorXf res = A.array() * xPrime.array() + B.array() - yPrime.array();

                return res; },
            .Stepper = [](const Eigen::VectorXf &value, const Eigen::VectorXf &model, const Eigen::VectorXf &valueWeight)
            { 
                    /* The model vector is split into two parts, A and B, where A*x' + B = y'
                     * x' and y' together constitutes all combinations of the elements of 'value',
                     * not considering direction of the combination with itself (so only the upper half matrix) */
                Eigen::VectorXf A = model.head(model.size()/2);
                Eigen::VectorXf B = model.tail(model.size()/2);

                /* Create all combinations of value's element (n*(n-1)/2 combinations)*/
                Eigen::VectorXf xPrime(A.size());
                Eigen::VectorXf yPrime(A.size());

                int index = 0;
                for (int i = 0; i < value.size(); ++i) {
                    for (int j = i + 1; j < value.size(); ++j) {
                        xPrime[index] = value[i];
                        yPrime[index] = value[j];
                        index++;
                    }
                }

                /* Calculate residual (or the inner gradient when deriving A*x' + B - y')*/
                Eigen::VectorXf inner = A.array() * xPrime.array() + B.array() - yPrime.array();

                /* Using squared error loss and doing partial derivates w.r.t. A and B 
                 * (A*x' + B - Y')² => 
                 * dE/dA = 2*(A*x' + B - y')*x' 
                 * dE/dB = 2*(A*x' + B - y') */
                Eigen::VectorXf aDelta = 2*inner.array()*xPrime.array();
                Eigen::VectorXf bDelta = 2*inner;

                /* Concatenate dE/dA and dE/dB back into original form of the model vector */
                Eigen::VectorXf delta(model.size());
                delta << aDelta , bDelta;

                return delta; },
            .Displayer = [](const Eigen::VectorXf &model) 
            { 
                return model; 
            },
            .Name = "Linear regression"};
    }
};