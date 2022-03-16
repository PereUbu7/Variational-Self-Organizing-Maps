#include "SOM.hpp"

#include "doctest/doctest.h"
#include <chrono>
#include <iostream>

int verbose = 0;
static constexpr size_t numberOfRuns = 100;

TEST_CASE("Complicated integration tests go here")
{
    CHECK(1 == 1);
}

TEST_CASE("Performance tests")
{
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    auto sut = Som{100, 100, 100};

    auto startTime = high_resolution_clock::now();

    SUBCASE("SOM::randomInitialize()")
    {
        for (int i{0}; i < numberOfRuns; ++i)
            sut.randomInitialize(1, 1);

        auto endTime = high_resolution_clock::now();

        auto ms_int = duration_cast<milliseconds>(endTime - startTime);

        std::cout << "Duration: " << ms_int.count() << "ms\n";

        CHECK(ms_int.count() < 1000);
    }
    SUBCASE("SOM:findBmu()")
    {
        sut.randomInitialize(1, 1);
        auto modelVector = Eigen::VectorXf{100};
        auto validityVector = std::vector<int>{100, 1};
        auto weightVector = std::vector<double>{100, 1.0};
        
        for (int i{0}; i < numberOfRuns; ++i)
            sut.findBmu(modelVector, validityVector, weightVector);

        auto endTime = high_resolution_clock::now();

        auto ms_int = duration_cast<milliseconds>(endTime - startTime);

        std::cout << "Duration: " << ms_int.count() << "ms\n";

        CHECK(ms_int.count() < 1000);
    }
}