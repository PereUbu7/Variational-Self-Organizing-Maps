#include "SOM.hpp"

#include "doctest/doctest.h"
#include <chrono>
#include <iostream>

int verbose = 0;
static constexpr size_t numberOfRuns = 1000;


TEST_CASE( "Complicated integration tests go here" )
{
    CHECK(1 == 1);
}

TEST_CASE("Performance test - SOM::randomInitialize")
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    auto sut = Som{100, 100, 100};

    auto startTime = high_resolution_clock::now();

    for(int i{0}; i<numberOfRuns; ++i)
        sut.randomInitialize(1, 1);

    auto endTime = high_resolution_clock::now();

    auto ms_int = duration_cast<milliseconds>(endTime - startTime);

    std::cout << "Duration: " << ms_int.count() << "ms\n";

    CHECK(ms_int.count() < 1000);
}