#include "Transformation.hpp"

Transformation Transformation::Standard(const std::vector<std::string> &columnNames)
{
    return Transformation {
        .Comparer =
            [](const T &value, const T &model, const T &dispersion, const T &valueWeight)
                { return model - value; },

        .Stepper =
            [](const T &value, const T &model, const T &valueWeight)
                { return value - model; },

        .names{columnNames},
        .Displayer =
            [columnNames](const T &model)
                {
                auto disp = std::vector<std::string>{};
                disp.reserve(columnNames.size());

                for (size_t index{0}; index < columnNames.size(); ++index)
                {
                    std::stringstream ss;
                    ss << columnNames[index] << " = " << model[index];

                    disp.emplace_back(ss.str());
                }
                return disp;
            },

        .Length =
        [](size_t vectorLength)
        {
            return vectorLength;
        },

        .Name = "Standard transformation"
    };
}

Transformation Transformation::StandardMedianEstimator(const std::vector<std::string> &columnNames)
{
    return Transformation {
        .Comparer =
            [](const T &value, const T &model, const T &dispersion, const T &valueWeight)
                { return model - value; },

        .Stepper =
            [](const T &value, const T &model, const T &valueWeight)
                { return (value - model).array().sign(); },

        .names{columnNames},
        .Displayer =
            [columnNames](const T &model)
                {
                auto disp = std::vector<std::string>{};
                disp.reserve(columnNames.size());

                for (size_t index{0}; index < columnNames.size(); ++index)
                {
                    std::stringstream ss;
                    ss << columnNames[index] << " = " << model[index];

                    disp.emplace_back(ss.str());
                }
                return disp;
            },

        .Length =
        [](size_t vectorLength)
        {
            return vectorLength;
        },

        .Name = "Standard median estimator transformation"
    };
}

Transformation Transformation::CombinatorialLinearRegression(const std::vector<std::string> &columnNames)
    {
        return Transformation{
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
                Eigen::VectorXf aDelta = -2*inner.array()*xPrime.array();
                Eigen::VectorXf bDelta = -2*inner;

                /* Concatenate dE/dA and dE/dB back into original form of the model vector */
                Eigen::VectorXf delta(model.size());
                delta << aDelta , bDelta;

                return delta; },
            .names{columnNames},
            .Displayer = [columnNames](const Eigen::VectorXf &model)
            { 
                auto disp = std::vector<std::string>{};
                disp.reserve(columnNames.size());

                const long tailIndex = static_cast<long>(model.size())/2;

                size_t index = 0;
                for (size_t i{0}; i < columnNames.size() && static_cast<long>(index) < tailIndex; ++i) {
                    for (size_t j{i + 1}; j < columnNames.size(); ++j) {
                        std::stringstream ss;
                        ss << columnNames[i] << " = " << model[index] << '*' << columnNames[j] << " + " << model[tailIndex + index];

                        disp.emplace_back(ss.str());
                        index++;
                    }
                }
                return disp; },
            .Length = [](size_t vectorLength)
            {
                return vectorLength*(vectorLength - 1u);
            },
            .Name = "Linear regression"};
    }