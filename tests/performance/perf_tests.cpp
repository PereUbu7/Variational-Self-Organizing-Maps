#include "SOM.hpp"
#include "Database.hpp"
#include "DataSet.hpp"

#include "dateTimeUtils.hpp"

#include <chrono>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>

int verbose = 1;

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
        SomTests() : sut{1, 1, 1},
                     numberOfRuns{100}
        {
        }

        auto test_train(const char *dbPath, const char *columnSpecPath, int weightDecayFunction)
        {
            auto db = DataBase{};
            auto dbOpenResult = db.open(dbPath);
            assert(dbOpenResult != 0);

            auto trainingSet = DataSet(columnSpecPath);

            trainingSet.loadDataBase(&db);

            sut = Som{100, 100, trainingSet.vectorLength()};

            sut.randomInitialize(std::time(NULL), 1);

            auto numOfEpochs = int{100};
            auto eta0 = double{0.01};
            auto etaDec = double{0.01};
            auto sigma0 = double{20.0};
            auto sigmaDec = double{0.01};

            auto startTime = high_resolution_clock::now();
            sut.train(&trainingSet, numOfEpochs, eta0, etaDec, sigma0, sigmaDec, weightDecayFunction);
            auto endTime = high_resolution_clock::now();
            return (endTime - startTime) / numOfEpochs;
        }

        template<size_t NumOfRuns>
        auto test_trainSingle(const int weightDecayFunction)
        {
            setup();

            sut.randomInitialize(std::time(NULL), 1);

            auto positions = std::array<size_t, NumOfRuns>();
            for(auto& element : positions)
                element = std::rand() % 100*100;
            
            auto samples = std::array<Eigen::VectorXf, NumOfRuns>();
            for(auto& sample : samples)
                sample = Eigen::VectorXf::Random(100);

            auto validityVector = Eigen::VectorXf::Ones(100);
            auto weightVector = Eigen::VectorXf::Ones(100);

            auto startTime = high_resolution_clock::now();
            
            for (size_t i{0}; i < NumOfRuns; ++i)
                sut.trainSingle(samples[i], validityVector, weightVector, 0.1, 50, positions[i], weightDecayFunction);

            auto endTime = high_resolution_clock::now();

            return endTime - startTime;
        }

        auto test_evaluate(const char *dbPath, const char *columnSpecPath)
        {
            auto db = DataBase{};
            auto dbOpenResult = db.open(dbPath);
            assert(dbOpenResult != 0);

            auto trainingSet = DataSet(columnSpecPath);

            trainingSet.loadDataBase(&db);

            sut = Som{100, 100, trainingSet.vectorLength()};

            sut.randomInitialize(std::time(NULL), 1);

            auto startTime = high_resolution_clock::now();
            for(size_t run = int{}; run < numberOfRuns; ++run)
                sut.evaluate(trainingSet);
            auto endTime = high_resolution_clock::now();
            return (endTime - startTime) / numberOfRuns;
        }
        auto test_randomInitialize()
        {
            setup();

            auto startTime = high_resolution_clock::now();

            for (size_t i{0}; i < numberOfRuns; ++i)
                sut.randomInitialize(1, 1);

            auto endTime = high_resolution_clock::now();

            return (endTime - startTime) / numberOfRuns;
        };
        template<size_t NumOfRuns>
        auto test_findBmu()
        {
            setup();

            sut.randomInitialize(std::time(NULL), 1);
            auto modelVector = Eigen::VectorXf::Random(100);
            auto validityVector = Eigen::VectorXf::Ones(100);
            auto weightVector = Eigen::VectorXf::Ones(100);

            auto startTime = high_resolution_clock::now();
            
            for (size_t i{0}; i < NumOfRuns; ++i)
                sut.findBmu(modelVector, validityVector, weightVector);

            auto endTime = high_resolution_clock::now();

            return endTime - startTime;
        }
        template<size_t NumOfRuns>
        auto test_findLocalBmu()
        {
            setup();

            sut.randomInitialize(std::time(NULL), 1);

            auto positions = std::array<size_t, NumOfRuns>();
            for(auto& element : positions)
                element = std::rand() % 100*100;

            auto modelVector = Eigen::VectorXf::Random(100);
            auto validityVector = Eigen::VectorXf::Ones(100);
            auto weightVector = Eigen::VectorXf::Ones(100);

            auto startTime = high_resolution_clock::now();
            
            for (size_t i{0}; i < NumOfRuns; ++i)
                sut.findLocalBmu(modelVector, validityVector, positions[i], weightVector);

            auto endTime = high_resolution_clock::now();

            return endTime - startTime;
        }
        template<size_t NumOfRuns>
        auto test_findRestrictedBmu()
        {
            setup();

            sut.randomInitialize(std::time(NULL), 1);

            auto positions = std::array<size_t, NumOfRuns>();
            for(auto& element : positions)
                element = std::rand() % 100*100;

            auto modelVector = Eigen::VectorXf::Random(100);
            auto validityVector = Eigen::VectorXf::Ones(100);
            auto weightVector = Eigen::VectorXf::Ones(100);

            auto startTime = high_resolution_clock::now();
            
            for (size_t i{0}; i < NumOfRuns; ++i)
                sut.findRestrictedBmu(modelVector, validityVector, 1, weightVector);

            auto endTime = high_resolution_clock::now();

            return endTime - startTime;
        }
        template<size_t NumOfRuns>
        auto test_findRestrictedBmd()
        {
            setup();

            sut.randomInitialize(std::time(NULL), 1);

            auto positions = std::array<size_t, NumOfRuns>();
            for(auto& element : positions)
                element = std::rand() % 100*100;

            auto modelVector = Eigen::VectorXf::Random(100);
            auto validityVector = Eigen::VectorXf::Ones(100);
            // auto validityVector = std::vector<int>(100);
            auto weightVector = Eigen::VectorXf::Ones(100);
            // auto weightVector = std::vector<float>(100);

            auto startTime = high_resolution_clock::now();
            
            for (size_t i{0}; i < NumOfRuns; ++i)
                sut.findRestrictedBmd(modelVector, validityVector, 0, weightVector);

            auto endTime = high_resolution_clock::now();

            return endTime - startTime;
        }
        template<size_t NumOfRuns>
        auto test_euclidianWeightedDist()
        {
            setup();

            sut.randomInitialize(std::time(NULL), 1);

            auto positions = std::array<size_t, NumOfRuns>();
            for(auto& element : positions)
                element = std::rand() % 100*100;

            auto modelVector = Eigen::VectorXf(100);
            auto validityVector = Eigen::VectorXf::Ones(100);
            auto weightVector = Eigen::VectorXf::Ones(100);

            double res = 0;

            auto startTime = high_resolution_clock::now();
            for (size_t i{0}; i < NumOfRuns; ++i)
                res += sut.euclidianWeightedDist(positions[i], modelVector, validityVector, weightVector)/1000000;

            auto endTime = high_resolution_clock::now();

            // Just to ensure compiler does not optimize 'res' away
            std::cout << res << '\n';

            return endTime - startTime;
        }

        auto test_measureSimilarity(const char *dbPath, const char *columnSpecPath)
        {
            auto db = DataBase{};
            auto dbOpenResult = db.open(dbPath);
            assert(dbOpenResult != 0);

            auto trainingSet = DataSet(columnSpecPath);

            trainingSet.loadDataBase(&db);

            sut = Som{100, 100, trainingSet.vectorLength()};

            sut.randomInitialize(std::time(NULL), 1);

            auto startTime = high_resolution_clock::now();
            for(size_t run = int{}; run < numberOfRuns; ++run)
                sut.measureSimilarity(&trainingSet, 1, 1);
            auto endTime = high_resolution_clock::now();
            return (endTime - startTime) / numberOfRuns;
        }

        auto test_updateUMatrix(const char *dbPath, const char *columnSpecPath)
        {
            auto db = DataBase{};
            auto dbOpenResult = db.open(dbPath);
            assert(dbOpenResult != 0);

            auto trainingSet = DataSet(columnSpecPath);

            trainingSet.loadDataBase(&db);

            sut = Som{100, 100, trainingSet.vectorLength()};

            sut.randomInitialize(std::time(NULL), 1);

            auto startTime = high_resolution_clock::now();
            for(size_t run = int{}; run < numberOfRuns; ++run)
                sut.updateUMatrix(trainingSet.getWeights());
            auto endTime = high_resolution_clock::now();
            
            return (endTime - startTime) / numberOfRuns;
        }

        auto test_variationalAutoEncoder(const char *dbPath, const char *columnSpecPath)
        {
            auto db = DataBase{};
            auto dbOpenResult = db.open(dbPath);
            assert(dbOpenResult != 0);

            auto trainingSet = DataSet(columnSpecPath);

            trainingSet.loadDataBase(&db);

            sut = Som{100, 100, trainingSet.vectorLength()};

            sut.randomInitialize(std::time(NULL), 1);

            auto startTime = high_resolution_clock::now();
            for(size_t run = int{}; run < numberOfRuns; ++run)
                sut.variationalAutoEncoder(&trainingSet, 0);
            auto endTime = high_resolution_clock::now();
            
            return (endTime - startTime) / numberOfRuns;
        }
    };

}
int main()
{
    using namespace Perftests;
    auto tester = SomTests{};

    printResultToFile(tester.test_randomInitialize(), "Som::randomInitialize()");
    printResultToFile(tester.test_train("./data/testDb.sq3", "./data/columnSpec.txt", 0), "Som::train(exponentialWeightDecay)");
    printResultToFile(tester.test_train("./data/testDb.sq3", "./data/columnSpec.txt", 1), "Som::train(inverseProportionalWeightDecay)");
    printResultToFile(tester.test_evaluate("./data/testDb.sq3", "./data/columnSpec.txt"), "Som::evaluate()");
    printResultToFile(tester.test_measureSimilarity("./data/testDb.sq3", "./data/columnSpec.txt"), "Som::measureSimilarity()");
    printResultToFile(tester.test_updateUMatrix("./data/testDb.sq3", "./data/columnSpec.txt"), "Som::updateUMatrix()");
    printResultToFile(tester.test_variationalAutoEncoder("./data/testDb.sq3", "./data/columnSpec.txt"), "Som::variationalAutoEncoder()");

    printResultToFile(tester.test_trainSingle<1000>(0), "Som:trainSingle(exponentialWeightDecay) per thousand");
    printResultToFile(tester.test_trainSingle<1000>(1), "Som:trainSingle(inverseProportionalWeightDecay) per thousand");
    printResultToFile(tester.test_euclidianWeightedDist<1000000>(), "Som:euclidianWeightedDist() per milion");
    printResultToFile(tester.test_findBmu<1000>(), "Som::findBmu() per thousand");
    printResultToFile(tester.test_findLocalBmu<1000000>(), "Som::findLocalBmu() per milion");
    printResultToFile(tester.test_findRestrictedBmu<1000>(), "Som::findRestrictedBmu() per thousand");
    printResultToFile(tester.test_findRestrictedBmd<1000>(), "Som::findRestrictedBmd() per thousand");
}
