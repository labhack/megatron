#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
//#include "fann.h"
//#include "floatfann.h"

int datapoints;
#define GOLDENRATIO 1.61803398875
#define DATARATE 480.0

void smooth(double sigma, double *data, int n)
{
	const double coefficient = 1.0 / (sqrt(2.0*M_PI)* sigma);
	const double denominator = 2.0*sigma*sigma;

	int half = (1.0/coefficient)+0.5;
	if(half<0)
	{
		half=0;
	}
	else if(half>(n-1)/2)
	{
		half=(n-1)/2;
	}
	int filterSize = 2*half+1;

	double *filter = new double[filterSize];
	double *tempFilteringArray = new double[n];

	double sum = 0.0;
	// Left half
	for(int i=0;i<half;i++)
	{
		filter[i] = exp(-((i-half)*(i-half))/denominator);
		sum += filter[i];
	}

	// Middle value
	filter[half]=1.0;

	// Right half
	for(int i=half+1; i<filterSize; i++)
	{
		filter[i]=filter[filterSize-i-1];
	}
	sum=2.0*sum+1.0;

	// Normalize filter
	for(int i=0; i<filterSize; i++)
	{
		filter[i]=filter[i]/sum;
	}

	// Smoothing
	for(int j=0; j<n; j++)
	{
		double value = 0.0;
		for(int k=0; k<filterSize; k++)
		{
			int nCol = j-half+k;

			if(nCol < 0)
			{
				nCol = 0;
			}
			else if(nCol > n-1)
			{
				nCol = n-1;
			}

			value += data[nCol]*filter[k];
		}

		tempFilteringArray[j] = value;
	}

	for(int i=0; i<n; i++)
	{
		data[i]=tempFilteringArray[i];
	}

	delete [] filter;
	delete [] tempFilteringArray;
}

int sign(double a)
{
	if(a < 0)
	{
		return -1;
	}
	else
	{
		return 1;
	}
}

void detectSaccades(double* arr, double* second, double** output, int mode)
{
	const double dt = 1.0/DATARATE;
	for(int i = 1; i < datapoints; i++)
	{
		double firstderivative = (arr[i] - arr[i-1])/dt;
		if(fabs(firstderivative) > 0.8)
		{
			// Save value of i
			int j = i;

			// Since derivative is high, move back to when it was lower to record amplitude
			double tempDeriv = (arr[i] - arr[i-1])/dt;
			while(fabs(tempDeriv) > 0.1)
			{
				i--;
				tempDeriv = (arr[i] - arr[i-1])/dt;
			}
			double startAmplitude = arr[i-1];

			// Move back to original trigger point
			i = j;

			// Analyzing data to find difference between common amplitudes
			double max = arr[i];
			double secondMax = second[i];
			double scanningFirstDerivative = (arr[i] - arr[i-1])/dt;
			bool run = true;
			int counter = 0;
			while(run && i < datapoints -1)
			{
				counter++;
				i++;
				scanningFirstDerivative = (arr[i] - arr[i-1])/dt;
				if(sign(firstderivative) != sign(scanningFirstDerivative) && fabs(startAmplitude) > fabs(arr[i]))
				{
					run = false;
				}

				if(max < arr[i])
				{
					max = arr[i];
				}

				if(secondMax < second[i])
				{
					secondMax = second[i];
				}
			}

			// Look at max amplitude between the delta time
			if(counter <= 96 )
			{
				//std::cerr << "Blink at " << i << " time=" << counter*dt << std::endl;
			}
			else
			{
				double angle = 0.0;
				if(mode == 0)
				{
					angle = atan2(max, secondMax);
				}
				else
				{
					angle = atan2(secondMax, max);
				}

				double magnitude = sqrt(max*max + secondMax*secondMax);

				std::cerr << "Index," << i << ",Angle," << angle << ",Magnitude," << magnitude << std::endl;
				output[i-1][0] = angle;
				output[i-1][1] = magnitude;
			}
		}
	}
}

int main(int argc, char** argv)
{
	datapoints = 0;

	if(argc < 2)
	{
		std::cerr << "No file name argument! Error!" << std::endl;
	}

	double sigma = 5.0;
	if(argc == 3)
	{
		sigma = atof(argv[2]);
	}
	std::cerr << "Using sigma of " << sigma << std::endl;

	// Read in lines from file
	double *veog;
	double *heog;

	char line[500];
	// Read Image
	std::fstream file(argv[1]);
	if(file.is_open())
	{
		std::cerr << "File open successful!" << std::endl;

		// Count lines in file
		file.getline(line, 500);
		std::string l;
		while(!file.eof())
		{
			file.getline(line, 500);
			datapoints++;
		}
		file.close();

		std::fstream file2(argv[1]);
		veog = new double[datapoints];
		heog = new double[datapoints];

		// Read header
		file2.getline(line, 137);
		for(int i = 0; i < datapoints; i++)
		{
			for(int j = 0; j < 9; j++)
			{
				std::string next;
				file2 >> next;
				if(j == 7)
				{
					veog[i] = atof(next.c_str());
				}
				else if (j == 8)
				{
					heog[i] = atof(next.c_str());
				}
			}
		}

		file2.close();
	}
	else
	{
		std::cerr << "Could not open file " << argv[1] << std::endl;
	}

	// Filter
	std::cerr << "Filtering.." << std::endl;
	smooth(sigma, veog, datapoints);
	smooth(sigma, heog, datapoints);

	// Detect Saccades
	std::cerr << "Detecting Saccades" << std::endl;
	double **veogOut = new double*[datapoints];
	double **heogOut = new double*[datapoints];
	for(int i = 0; i < datapoints; i++)
	{
		veogOut[i] = new double[2];
		heogOut[i] = new double[2];
		veogOut[i][0] = 0;
		veogOut[i][1] = 0;
		heogOut[i][0] = 0;
		heogOut[i][1] = 0;
	}
	detectSaccades(veog, heog, veogOut, 0);
	detectSaccades(heog, veog, heogOut, 1);

	// Write out Filtered Data
	std::cerr << "Writing filtered data!" << std::endl;
	std::ofstream outfile("filtered_output.csv");
	if(outfile.is_open())
	{
		outfile << "filtered_veog,filtered_heog,veog_angle,veog_mag,heog_angle,heog_mag" << std::endl;
		for(int i = 0; i < 10000; i++)
		{
			outfile << veog[i] << "," << heog[i] << "," << veogOut[i][0] << "," << veogOut[i][1] << "," << heogOut[i][0] << "," << heogOut[i][1] << std::endl;
		}
		outfile.close();
	}

	/*
	// Write out Filtered Data for neural nets
	std::cerr << "Writing filtered data!" << std::endl;
	std::ofstream outfile("filtered_output.csv");
	if(outfile.is_open())
	{
		outfile << "21 480 1" << std::endl;
		for(int i = 0; i < 21; i++)
		{
			for(int j = 0; j < DATARATE; j++)
			{
				outfile << veog[(int)(i*DATARATE) + j] << " ";
			}
			outfile << std::endl << "" << std::endl;
		}

		outfile.close();
	}
*/

	/*
	// Testing out ANNs
	const unsigned int num_input = 480;
	const unsigned int num_output = 1;
	const unsigned int num_layers = 3;
	const unsigned int num_neurons_hidden = 240;
	const float desired_error = (const float) 0.1;
	const unsigned int max_epochs = 25000;
	const unsigned int epochs_between_reports = 100;

	struct fann *ann = fann_create_standard(num_layers, num_input, num_neurons_hidden, num_output);

	fann_set_training_algorithm(ann, FANN_TRAIN_BATCH);

	fann_set_activation_function_hidden(ann, FANN_SIGMOID_SYMMETRIC);
	fann_set_activation_function_output(ann, FANN_SIGMOID_SYMMETRIC);

	fann_set_activation_steepness_hidden(ann, 0.1);

	fann_train_on_file(ann, "nn21.data", max_epochs, epochs_between_reports, desired_error);

	fann_save(ann, "nn_float.net");

	fann_destroy(ann);

	struct fann *annRun = fann_create_from_file("nn_float.net");

	float small[480];
	float* output;
	for(int i = 2; i < 200; i++)
	{
		for(int j = 0; j < DATARATE; j++)
		{
			small[i*480 + j] = (float)(veog[i*480 + j]);
		}

		output = fann_run(annRun, small);
        printf("NN test -> %f\n", output[0]);
	}

	fann_destroy(annRun);
*/
	for(int i = 0; i < datapoints; i++)
	{
		delete [] veogOut[i];
		delete [] heogOut[i];
	}
	delete [] veog;
	delete [] heog;
	delete [] veogOut;
	delete [] heogOut;

	std::cerr << "Exiting!" << std::endl;
	return 0;
}
