#include "SOM.hpp"
#include "dateTimeUtils.hpp"

#include <chrono>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>

int verbose = 0;

namespace Perftests
{
    using namespace std::chrono_literals;
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    void printResultToFile(auto duration, std::string functionName)
    {
        std::ofstream myfile;
        myfile.open("example.txt", std::ios::out | std::ios::app);
        myfile << DateTime::Utils::CurrentTimeStr() << '\t' << functionName << '\t' << duration_cast<std::chrono::milliseconds>(duration).count() << "ms\n";
        myfile.close();
    }

    class SomTests
    {
    private:
        Som sut;
        size_t numberOfRuns;
        void setup() { sut = Som{100, 100, 100}; };
        void teardown();

    public:
        SomTests() :
            sut{1, 1, 1},
            numberOfRuns{100} 
        {}

        auto test_randomInitialize()
        {
            setup();

            auto startTime = high_resolution_clock::now();

            for (int i{0}; i < numberOfRuns; ++i)
                sut.randomInitialize(1, 1);

            auto endTime = high_resolution_clock::now();

            return (endTime - startTime) / numberOfRuns;
        };
        auto test_findBmu()
        {
            setup();

            auto startTime = high_resolution_clock::now();

            sut.randomInitialize(1, 1);
            auto modelVector = Eigen::VectorXf{100};
            auto validityVector = std::vector<int>{100, 1};
            auto weightVector = std::vector<double>{100, 1.0};

            for (int i{0}; i < numberOfRuns; ++i)
                sut.findBmu(modelVector, validityVector, weightVector);

            auto endTime = high_resolution_clock::now();

            return (endTime - startTime) / numberOfRuns;
        };
    };

}
int main()
{
    using namespace Perftests;
    auto tester = SomTests{};

    printResultToFile(tester.test_randomInitialize(), "Som::randomInitialize()");
    printResultToFile(tester.test_findBmu(), "Som::findBmu()");
}
