#include "Transformation.hpp"
#include "SqliteDataLoader.hpp"

#include <iostream>
#include <string>

#include "Eigen/Dense"

#include "doctest/doctest.h"

TEST_CASE("dummies")
{
    /* Call method in SqliteDataLoad.cpp for correct linkage */
    auto a = SqliteDataLoader("path");
    a.open();
}

TEST_CASE("chatGTP")
{
     Eigen::Vector4f x;
    x << 1, 2, 3, 4;

    // Calculate the size of the resulting vectors
    int size_x_prime = static_cast<int>(x.size()) * (static_cast<int>(x.size()) - 1) / 2;

    // Initialize vectors x_prime and y_prime with the correct sizes
    Eigen::VectorXf x_prime(size_x_prime);
    Eigen::VectorXf y_prime(size_x_prime);

    // Construct x_prime and y_prime
    int index = 0;
    for (int i = 0; i < x.size(); ++i) {
        for (int j = i + 1; j < x.size(); ++j) {
            x_prime[index] = x[i];
            y_prime[index] = x[j];
            index++;
        }
    }

    // Output x_prime and y_prime
    std::cout << "x_prime = " << x_prime.transpose() << std::endl;
    std::cout << "y_prime = " << y_prime.transpose() << std::endl;
}


TEST_CASE("Fakes")
{
    Eigen::VectorXf P(8);
    P << 1, 2, 3, 4, 5, 6, 7, 8;

    Eigen::VectorXf X(4);
    X << 1, 2, 3, 4;

    auto names = std::vector<std::string>{};

    auto transformation = Transformation{
        .Comparer =
        [](const Eigen::VectorXf& X, const Eigen::VectorXf& P, const Eigen::VectorXf &dispersion, const Eigen::VectorXf &valid)
        { 
            auto A = P.head(P.size()/2);
            auto B = P.tail(P.size()/2);

            Eigen::VectorXf Y = (A.array() * X.array() + B.array()).matrix();

            return Y;
        },
        .Stepper =
        [](const Eigen::VectorXf& X, const Eigen::VectorXf& P, const Eigen::VectorXf& valid)
        { 
            auto A = P.head(P.size()/2);
            auto B = P.tail(P.size()/2);

            Eigen::VectorXf Y = (A.array() * X.array() + B.array()).matrix();

            return Y;
        },
        .Displayer =
        [&names](const Eigen::VectorXf& v)
        { 
            // auto A = v.head(v.size()/2);

            return names;
        }
        };

    std::cout << P.transpose() << "\n" << X.transpose() << "\n";
    std::cout << transformation.Comparer(X, P, P, X).transpose() << '\n';
    std::cout << transformation.Stepper(X, P, X).transpose() << '\n';
    // std::cout << transformation.Displayer(names) << '\n';
}