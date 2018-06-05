#include "SOM.hpp"

void displayHelp();

using namespace std;

int verbose = 0;

int main(int argc, char *argv[])
{
	char *ptr;
	
	unsigned int numOfEpochs = 0;
	unsigned int somWidth = 0;
	unsigned int somHeight = 0;
	
	unsigned int allowedNumberOfStdDev = 0;
	unsigned int minBmuHits = 0;
	
	int success = FALSE;
	
	double eta0 = 0;
	double etaDec = 0;
	double sigma0 = 0;
	double sigmaDec = 0;
	double initSigma = 0;
	
	int weightDecayFuncion = 0;
	
	int trainingSOM = 0;
	int evaluatingSOM = 0;
	int measuringSOM = 0;
	int autoencoderSOM = 0;
	int variationalAutoEncoderSOM = 0;
	DataBase db;
	
	if( argc < (ARG_SETTING + 1) )
	{
		displayHelp();
		exit(-1);
	}
	else if( std::strcmp(argv[ARG_VERBOSE], "-v") == 0 )
	{
		std::cout << "Setting verbose flag\n";
		verbose = 1;
	}
	else if( std::strcmp(argv[ARG_VERBOSE], "-q") == 0 )
		verbose = 0;
	else
	{
		displayHelp();
		exit(-1);
	}
	
	// Set task to "train"
	if( std::strcmp(argv[ARG_SETTING], "-t") == 0 )
	{
		if(verbose) std::cout << "Training SOM\n";
		trainingSOM = 1;
	}
	// Set task to "evaluate"
	else if( std::strcmp(argv[ARG_SETTING], "-e") == 0 )
	{
		if(verbose) std::cout << "Evaluating SOM\n";
		evaluatingSOM = 1;
	}
	// Set task to "measure"
	else if( std::strcmp(argv[ARG_SETTING], "-m") == 0 )
	{
		if(verbose) std::cout << "Measuring SOM\n";
		measuringSOM = 1;
	}
	else if( std::strcmp(argv[ARG_SETTING], "-ae") == 0 )
	{
		if(verbose) std::cout << "AutoEncoder SOM\n";
		autoencoderSOM = 1;
	}
	else if( std::strcmp(argv[ARG_SETTING], "-vae") == 0 )
	{
		if(verbose) std::cout << "Variational AutoEncoder SOM\n";
		variationalAutoEncoderSOM = 1;
	}
	else
	{
		displayHelp();
		exit(-1);
	}
	
	// Check that we have correct number of arguments for corresponding task
	if( (trainingSOM && argc != (ARG_SOM_WEIGHT_DECAY_FUNCTION + 1 )) || (evaluatingSOM  && argc != (ARG_SOM_FILE + 1)) || (measuringSOM && argc != (ARG_ALLOWED_STD_DEV + 1)) )
	{
		displayHelp();
		exit(-1);
	}
	
	if( trainingSOM )
	{
		// Set values
		somHeight = std::strtoul(argv[ARG_SOM_HEIGHT], &ptr, 10);
		somWidth = std::strtoul(argv[ARG_SOM_WIDTH], &ptr, 10);
		
		eta0 = std::strtod(argv[ARG_SOM_ETA0], &ptr);
		etaDec = std::strtod(argv[ARG_SOM_ETA_DEC], &ptr);
		sigma0 = std::strtod(argv[ARG_SOM_SIGMA0], &ptr);
		sigmaDec = std::strtod(argv[ARG_SOM_SIGMA_DEC], &ptr);
		
		numOfEpochs = std::strtoul(argv[ARG_SOM_EPOCHS], &ptr, 10);
		initSigma = std::strtod(argv[ARG_SOM_INIT_SIGMA], &ptr);
		
		weightDecayFuncion = std::strtoul(argv[ARG_SOM_WEIGHT_DECAY_FUNCTION], &ptr, 10);
		
		// Check of any value is zero or not set
		if( !somHeight || !somWidth || !eta0 || !etaDec || !sigma0 || !sigmaDec || !numOfEpochs || !initSigma )
		{
			displayHelp();
			exit(-1);
		}
		
		if( db.open(argv[ARG_DB_FILE]) == 0 )
			return 0;
			
		DataSet trainingSet("icanNames.txt");
		//trainingSet.loadTextFile("7246051 efter rep2.txt");
		//trainingSet.loadTextFile("7246037 20s.txt");
		
		trainingSet.loadDataBase(&db);
		
		trainingSet.shuffle();
		
		Som test(somWidth, somHeight, trainingSet.vectorLength());
		
		test.randomInitialize(std::time(NULL), initSigma);
		
		if(verbose) trainingSet.display();
		
		test.train(&trainingSet, numOfEpochs, eta0, etaDec, sigma0, sigmaDec, weightDecayFuncion);
		
		if(verbose) test.display();
		
		// Fixa så att denna klarar 1xX map eller Xx1. 
		// Lägg till en inputparameter för huruvida man vill skapa denna
		if( somWidth > 2 && somHeight > 2 )
			test.updateUMatrix(&trainingSet);
		
		if(verbose) test.displayUMatrix();
		
		test.save(argv[ARG_SOM_FILE]);
	}
	else if( evaluatingSOM )
	{
		if( db.open(argv[ARG_DB_FILE]) == 0 )
			return 0;
			
		DataSet devSet("icanNames.txt");
		
		//devSet.loadTextFile("7246051 efter rep.txt");
		devSet.loadDataBase(&db);
		
		Som test(argv[ARG_SOM_FILE]);
		
		if(verbose) test.display();
		if(verbose) test.displayUMatrix();
		
		printf("Error on dataset: %f\n", test.evaluate(&devSet));
	}
	else if( measuringSOM )
	{
		allowedNumberOfStdDev = std::strtoul(argv[ARG_ALLOWED_STD_DEV], &ptr, 10);
		
		minBmuHits = std::strtoul(argv[ARG_MIN_BMU_HITS], &ptr, 10);
		
		if( db.open(argv[ARG_DB_FILE]) == 0 )
			return 0;
			
		DataSet testSet("icanNames.txt");
		
		testSet.loadDataBase(&db);
		
		Som test(argv[ARG_SOM_FILE]);
		
		success = test.measureSimilarity(&testSet, allowedNumberOfStdDev, minBmuHits);
	}
	else if( autoencoderSOM )
	{
		if( db.open(argv[ARG_DB_FILE]) == 0 )
			return 0;
		
		minBmuHits = std::strtoul(argv[ARG_MIN_BMU_HITS], &ptr, 10);
		
		DataSet testSet("icanNames.txt");
		
		testSet.loadDataBase(&db);
		
		Som test(argv[ARG_SOM_FILE]);
		success = test.autoEncoder(&testSet, minBmuHits);
	}
	else if( variationalAutoEncoderSOM )
	{
		if( db.open(argv[ARG_DB_FILE]) == 0 )
			return 0;
		
		minBmuHits = std::strtoul(argv[ARG_MIN_BMU_HITS], &ptr, 10);
		
		DataSet testSet("icanNames.txt");
		
		testSet.loadDataBase(&db);
		
		Som test(argv[ARG_SOM_FILE]);
		success = test.variationalAutoEncoder(&testSet, minBmuHits);
	}
	
	
	
	/*data.loadTextFile("7246037 20s.txt");
	data.loadTextFile("7246055 20s.txt");
	data.loadTextFile("7246040 20s.txt");
	data.loadTextFile("7246052 20s.txt");
	data.loadTextFile("7246036 20s.txt");
	*/
	
	/*std::cout << "X: " << test.findBmu(trainingSet.getData(0)).getX() << "\n" <<
				 "Y: " << test.findBmu(trainingSet.getData(0)).getY() << "\n" <<
				 "D: " << test.euclidianDist(test.findBmu(trainingSet.getData(0)), trainingSet.getData(0)) << "\n";
	
	printf("Error on training set: %f\n", test.evaluate(devSet));
	
	*/
	
	//Som bla("test.txt");
	//bla.display();
	//bla.displayUMatrix();
	 
    //system("PAUSE");
    return success;
}

void displayHelp()
{
	std::cout << "Variational Self Organizing Maps implemented by Martin Noring 2017 version " << VERSION << "\n\n";
	
	std::cout << "Program call:\n";
	std::cout << "som.exe -(verbose flag) -(task flag) [file path to database] [file path to SOM] [eventual parameters][][]...\n\n";
	
	std::cout << "Verbose flag:\n";
	std::cout << "\t-v\tverbose mode\n";
	std::cout << "\t-q\tquiet mode\n\n";
	
	std::cout << "Task Flag:\n";
	std::cout << "\t-t\ttrain SOM\n";
	std::cout << "\t-e\tevaluate SOM\n";
	std::cout << "\t-m\tmeasure SOM\n\n";
	
	std::cout << "Training:\n";
	std::cout << "In order to train a SOM, nine(9) different parameters need to be specified\n";
	std::cout << "\tSom height <unsigned integer>\n";
	std::cout << "\tSom width <unsigned integer>\n";
	std::cout << "\tInitial learning rate <float>\n";
	std::cout << "\tLearning rate exponential decay <float>\n";
	std::cout << "\tInitial neighbourhood width, sigma <float>\n";
	std::cout << "\tNeighbourhood width exponential decay <float>\n";
	std::cout << "\tNumber of training epochs <unsigned integer>\n";
	std::cout << "\tSOM random initialization sigma <float>\n";
	std::cout << "\tWeight decay function to be used <int> (0 = exponential, 1 = inverse proportional (non-parametric)\n";
	std::cout << "Those parameters must be specified in the above order.\n\n";
	
	std::cout << "Evaluating:\n";
	std::cout << "\tThis command does not need any extra parameters\n";
	std::cout << "\tThe evaluation procedure results in a sum over all residual Euclidian distances of the entire evaluation set\n\n";
	
	std::cout << "Measuring:\n";
	std::cout << "In order to measure on a SOM, three(3) additional parameters need to be specified\n";
	std::cout << "\tMinimum allowed BMU hits in map in order to consider a point. Set to zero(0) to consider all points <unsigned integer>\n";
	std::cout << "\tNumber of allowed standard deviations <unsigned integer>\n";
	std::cout << "\tThe measuring procedure prints the residuals of all elements over the entire measuring set\n\n";
	
	std::cout << "Example calls:\n";
	std::cout << "\tTrain SOM: -v -t trainingSet.db learnedSOM.txt 12 12 0.9 0.5 6 0.3 20 1 1\n";
	std::cout << "\tEvaluate SOM: -v -e devSet.db learnedSOM.txt\n";
	std::cout << "\tMeasure SOM: -q -m cloudPass.db learnedSOM.txt 10 3\n";
	std::cout << "\tAutoEncoder SOM: -q -ae cloudPass.db learnedSOM.txt 10\n";
	std::cout << "\tVariationalAutoEncoder SOM: -v -vae cloudPass.db learnedSOM.txt 10\n\n";
}