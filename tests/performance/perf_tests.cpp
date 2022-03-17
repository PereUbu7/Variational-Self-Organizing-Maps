#include "SOM.hpp"

#include <chrono>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>

int verbose = 0;
static constexpr size_t numberOfRuns = 100;

using namespace std::chrono_literals;

std::string FormatTime(std::chrono::system_clock::time_point tp)
{
    std::stringstream ss;
    auto t = std::chrono::system_clock::to_time_t(tp);
    auto tp2 = std::chrono::system_clock::from_time_t(t);
    if (tp2 > tp)
        t = std::chrono::system_clock::to_time_t(tp - std::chrono::seconds(1));
    ss << std::put_time(std::localtime(&t), "%Y-%m-%d %T")
       << "." << std::setfill('0') << std::setw(3)
       << (std::chrono::duration_cast<std::chrono::milliseconds>(
               tp.time_since_epoch())
               .count() %
           1000);
    return ss.str();
}

std::string CurrentTimeStr()
{
    return FormatTime(std::chrono::system_clock::now());
}

void printResultToFile(auto duration, std::string functionName)
{
    std::ofstream myfile;
    myfile.open("example.txt", std::ios::out | std::ios::app);
    myfile << CurrentTimeStr() << '\t' << functionName << '\t' << duration_cast<std::chrono::milliseconds>(duration).count() << "ms\n";
    myfile.close();
}


auto test_randomInitialize()
{
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    auto sut = Som{100, 100, 100};

    auto startTime = high_resolution_clock::now();

    for (int i{0}; i < numberOfRuns; ++i)
        sut.randomInitialize(1, 1);

    auto endTime = high_resolution_clock::now();

    return (endTime - startTime)/numberOfRuns;
}
auto test_findBmu()
{
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;

    auto sut = Som{100, 100, 100};

    auto startTime = high_resolution_clock::now();

    sut.randomInitialize(1, 1);
    auto modelVector = Eigen::VectorXf{100};
    auto validityVector = std::vector<int>{100, 1};
    auto weightVector = std::vector<double>{100, 1.0};

    for (int i{0}; i < numberOfRuns; ++i)
        sut.findBmu(modelVector, validityVector, weightVector);

    auto endTime = high_resolution_clock::now();

    return (endTime - startTime)/numberOfRuns;
}

int main()
{
    std::cout << CurrentTimeStr() << std::endl;

    printResultToFile(test_randomInitialize(), "Som::randomInitialize()");
    printResultToFile(test_findBmu(), "Som::findBmu()");
}
