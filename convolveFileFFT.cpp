
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <malloc.h>
#include <cstdint>
#include <iostream>
#include <cstdio>
#include <ctime>
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr


using namespace std;


// char *format;
// int chunk_size;
// char *chunk_id;
// char *subchunk1_id;
int size1;
// short block_align;
// short bits_per_sample;
// short audio_format;
// short num_channels;
// int sample_rate;
char *subchunk2_id;
int data_size;
short* file_data;
// int subchunk1_size;
// int byte_rate;

typedef struct {
	char *chunk_id;
	int chunk_size;
	char *format;
	char *subchunk1_id;
	int subchunk1_size;
	short audio_format;
	short num_channels;
	int sample_rate;
	int byte_rate;
	short block_align;
	short bits_per_sample;
} Wavheader;

Wavheader header;

//function definitions
void wavWrite(char *fileName, int numSamples, float *signal);
float* wavRead(char *fileName, float *signal, int *Thesize);
void convolve(float x[], int N, float h[], int M, float y[], int P);
void four1(double data[], int nn, int isign);
void FFTScale (double signal[], int N);



int main(int argc, char* args[])
{
	std::clock_t start;
    double duration;
    start = std::clock();
	if (argc!= 4){ //check if we hve 3 command line arguments
		
		printf("Please enter an input file, and IR file and an output file name." );
		return 0;
	}

	char *inputFileName = args[1];
	char *IRFileName = args[2];
	char *outputFileName = args[3];
	float *inFileSignal;
	int inFileSignalSize;
	float *IRFileSignal;
	int IRFileSignalSize;
	inFileSignal = wavRead(inputFileName, inFileSignal, &size1 );
	inFileSignalSize = size1;
	IRFileSignal = wavRead(IRFileName, IRFileSignal, &size1);
	IRFileSignalSize = size1;
	int outFileSignalSize = inFileSignalSize + IRFileSignalSize - 1;
	float *outFileSignal = new float[outFileSignalSize];
	printf("Convolving...");
	convolve(inFileSignal, inFileSignalSize, IRFileSignal, IRFileSignalSize, outFileSignal, outFileSignalSize);
	//scale output below
	float min = 0, max = 0;
	int i = 0;

	for(i = 0; i < outFileSignalSize; i++)
	{
		if(outFileSignal[i] > max)
			max = outFileSignal[i];
		if(outFileSignal[i] < min)
			min = outFileSignal[i];
	}

	min = min * -1;
	if(min > max)
		max = min;
	for(i = 0; i < outFileSignalSize; i++)
	{
		outFileSignal[i] = outFileSignal[i] / max;
	}
	wavWrite(outputFileName, outFileSignalSize, outFileSignal);
	printf("Written to file: %s\n", outputFileName);
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	printf("The duration of the program is: ");
	std::cout<< duration <<" seconds";
	return 0;
}

void wavWrite(char *fileName, int numSamples, float *signal)
{
	ofstream outFile( fileName, ios::out | ios::binary);
	//  Calculate the total number of bytes for the data chunk  
	int chunk_size = header.num_channels * numSamples * (header.bits_per_sample / 8);
	header.chunk_id = "RIFF";
	outFile.write( header.chunk_id, 4);
	outFile.write( (char*) &chunk_size, 4);
	header.format = "WAVE";
	outFile.write( header.format, 4);
	outFile.write( header.subchunk1_id, 4);
	header.subchunk1_size = 16;
	outFile.write( (char*) &header.subchunk1_size, 4);
	outFile.write( (char*) &header.audio_format, 2);
	outFile.write( (char*) &header.num_channels, 2);
	outFile.write( (char*) &header.sample_rate, 4);
	outFile.write( (char*) &header.byte_rate, 4);
	outFile.write( (char*) &header.block_align, 2);
	outFile.write( (char*) &header.bits_per_sample, 2);
	outFile.write( subchunk2_id, 4);
	data_size = numSamples * 2;
	outFile.write( (char*)&data_size, 4);
	short val;
	for(int i = 0; i < numSamples; i++)
	{
		val = (short)(signal[i] * (pow(2,15) - 1));
		outFile.write((char*)&val, 2);
	}
	outFile.close();
}


float* wavRead(char *fileName, float *signal, int *Thesize)
{
	ifstream inputFile( fileName, ios::in | ios::binary);
	inputFile.seekg(ios::beg);
	header.chunk_id = new char[4];
	inputFile.read( header.chunk_id, 4);
	inputFile.read( (char*) &header.chunk_size, 4);
	header.format = new char[4];
	inputFile.read( header.format, 4);
	header.subchunk1_id = new char[4];
	inputFile.read( header.subchunk1_id, 4);
	inputFile.read( (char*) &header.subchunk1_size, 4);
	inputFile.read( (char*) &header.audio_format, 2);
	inputFile.read( (char*) &header.num_channels, 2);
	inputFile.read( (char*) &header.sample_rate, 4);
	inputFile.read( (char*) &header.byte_rate, 4);
	inputFile.read( (char*) &header.block_align, 2);
	inputFile.read( (char*) &header.bits_per_sample, 2);

	if(header.subchunk1_size == 18)
	{
		char *garbage;
		garbage = new char[2];
		inputFile.read( garbage, 2);
	}

	subchunk2_id = new char[4];
	inputFile.read( subchunk2_id, 4);
	//DataSize
	inputFile.read( (char*)&data_size, 4);
	//GetData
	*Thesize = data_size / 2;
	int size = data_size / 2;
	file_data = new short[size];
	for(int j = 0 ; j < size; j++)
	{
		inputFile.read((char*) &file_data[j], 2);
	}

	short val;
	signal = new float[size];
	for(int i = 0; i < size; i++)
	{
		val = file_data[i];
		signal[i] = (val * 1.0) / (pow(2,15) - 1);
		if(signal[i] < -1.0)
			signal[i] = -1.0;

	}
	inputFile.close();
	return signal;
}


//uses FFT algorithm to convolve
void convolve(float x[], int N, float h[], int M, float y[], int P)
{
	int myArraySize = 1;
	int i = 0;
	// For FFT we need array size of a power of 2
	while (myArraySize < P) {
		myArraySize *= 2;
	}

	double *paddedInput = new double[2 * myArraySize];
	for (i = 0; i < (N * 2); i+=2) {
		paddedInput[i] = x[i/2];
		paddedInput[i+1] = 0;
	}
	for (; i <myArraySize; i++) {
		paddedInput[i] = 0;
	}
	
	double *paddedImpulseResponse = new double[2 * myArraySize];
	for (i = 0; i < (M * 2); i+=2) {
		paddedImpulseResponse[i] = h[i/2];
		paddedImpulseResponse[i+1] = 0;
	}
	for (; i < myArraySize; i++) {
		paddedImpulseResponse[i] = 0;
	}
	
	double *paddedOutput = new double[2 * myArraySize];
	for (i = 0; i < myArraySize; i++) {
		paddedOutput[i] = 0;
	}
	four1((paddedInput - 1), myArraySize, 1);
	four1((paddedImpulseResponse - 1), myArraySize, 1);
	for (i = 0; i < (myArraySize * 2); i+=2) {
		paddedOutput[i] = (paddedInput[i] * paddedImpulseResponse[i]) - (paddedInput[i+1] * paddedImpulseResponse[i+1]);
		paddedOutput[i+1] = (paddedInput[i+1] * paddedImpulseResponse[i]) + (paddedInput[i] * paddedImpulseResponse[i+1]);
	}
	four1((paddedOutput - 1), myArraySize, -1);
	
	// FFT scaling.. we need to scale as per class notes
	FFTScale(paddedOutput, myArraySize);
	
	// removing padding
	for (i = 0; i < P; i++) {
		y[i] = paddedOutput[i*2];
	}
}

//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below)

void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    float wtemp, wr, wpr, wpi, wi, theta;
    float tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
		if (j > i) {
			SWAP(data[j], data[i]);
			SWAP(data[j+1], data[i+1]);
		}
		m = nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
    }

    mmax = 2;
    while (n > mmax) {
		istep = mmax << 1;
		theta = isign * (6.28318530717959 / mmax); //changed just this
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr * data[j] - wi * data[j+1];
				tempi = wr * data[j+1] + wi * data[j];
				data[j] = data[i] - tempr;
				data[j+1] = data[i+1] - tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr = (wtemp = wr) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
    }
}

//the output from either the FFTs or the IFFT (but not both) will have to be scaled by dividing by each data point by N
//this algorithm is taken from class Notes
void FFTScale (double x[], int N)
{
	int k;
	int i;
	for (k = 0, i = 0; k < N; k++, i+=2) {
		x[i] /= (float)N;
		x[i+1] /= (float)N;
	}
}



