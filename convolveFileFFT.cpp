
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


// Name: Anish Pokhrel
// UCID: 30115576
// CPSC 501 (Fall 2023) Assignment 4


int size1;
char subchunk2_id[4];
int data_size;
short* file_data;

const int subchunk1_size_max = 18;

typedef struct {
	char *chunk_id;
	int chunk_size;
	char *format;
	char subchunk1_id[4];
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
void convolve(float x[], int N, float h[], int M, float y[], int P);
void four1(double data[], int nn, int isign);
void writeWavFile(char *file_name, int num_samples, float *signal);
float* readWavFile(char *file_name, float *signal, int *combined_size);




int main(int argc, char* argv[]){

	// begin timer
	clock_t startingTime = clock();

	// Process the command line arguments
	if (argc <= 3){
		printf("Error: must enter an input file, an, IR file, and an output file name.\n");
		fprintf(stderr, "Usage:  %s input_file ouput_file\n", argv[0]);
		exit(-1);
	}

	// read input and IR files
	char *input_file = argv[1];
	char *IR_file = argv[2];
	char *ouput_file = argv[3];


	// Initializing file signals
	float *file_signal, *file_signal_IR, *outfile_signal;
	int file_signal_size, file_signal_IR_size, outfile_signal_size;

	// reding raw file
	file_signal = readWavFile(input_file, file_signal, &size1 );
	file_signal_size = size1;

	// reading IR file 
	file_signal_IR = readWavFile(IR_file, file_signal_IR, &size1);
	file_signal_IR_size = size1;

	outfile_signal_size = file_signal_size + file_signal_IR_size - 1;
	outfile_signal = new float[outfile_signal_size];

	// begin convolving
	convolve(file_signal, file_signal_size, file_signal_IR, file_signal_IR_size, outfile_signal, outfile_signal_size);
	
	float smallest_size = 0;
	float largest_size = 0;

	int i = 0;
	while (i < outfile_signal_size) {
		if(outfile_signal[i] > largest_size) {
			largest_size = outfile_signal[i];
		}
			
		if(outfile_signal[i] < smallest_size) {
			smallest_size = outfile_signal[i];
		}
		i++;
	}	

	smallest_size = smallest_size * -1;
	if(smallest_size > largest_size)
		largest_size = smallest_size;
	for(int i = 0; i < outfile_signal_size; i++){
		outfile_signal[i] = outfile_signal[i] / largest_size;
	}

	// write to file 
	writeWavFile(ouput_file, outfile_signal_size, outfile_signal);
	printf("Written to file: %s\n", ouput_file);

	// Elapsed time of the program
	double time_elapsed = clock() - startingTime;
	printf("FFT convolution program finished in %f seconds\n\n", time_elapsed/CLOCKS_PER_SEC);
	return 0;
}

void writeWavFile(char *file_name, int num_samples, float *signal) {
	ofstream outputFileStream( file_name, ios::out | ios::binary);

	//  Calculate the total number of bytes for the data chunk  
	int chunk_size = header.num_channels * num_samples * (header.bits_per_sample / 8);
	header.chunk_id = "RIFF";
	outputFileStream.write( header.chunk_id, 4);
	outputFileStream.write( (char*) &chunk_size, 4);
	header.format = "WAVE";
	outputFileStream.write( header.format, 4);
	outputFileStream.write( header.subchunk1_id, 4);
	header.subchunk1_size = 16;
	outputFileStream.write( (char*) &header.subchunk1_size, 4);
	outputFileStream.write( (char*) &header.audio_format, 2);
	outputFileStream.write( (char*) &header.num_channels, 2);
	outputFileStream.write( (char*) &header.sample_rate, 4);
	outputFileStream.write( (char*) &header.byte_rate, 4);
	outputFileStream.write( (char*) &header.block_align, 2);
	outputFileStream.write( (char*) &header.bits_per_sample, 2);
	outputFileStream.write( subchunk2_id, 4);
	data_size = num_samples * 2;
	outputFileStream.write( (char*)&data_size, 4);
	short temp;

	for(int i = 0; i < num_samples; i++) {
		temp = (short)(signal[i] * (pow(2,15) - 1));
		outputFileStream.write((char*)&temp, 2);
	}
	outputFileStream.close();
}


float* readWavFile(char *file_name, float *signal, int *combined_size) {
	ifstream inputFileStream( file_name, ios::in | ios::binary);
	inputFileStream.seekg(ios::beg);
	header.chunk_id = new char[4];
	inputFileStream.read( header.chunk_id, 4);
	inputFileStream.read( (char*) &header.chunk_size, 4);
	header.format = new char[4];
	inputFileStream.read( header.format, 4);

	inputFileStream.read( header.subchunk1_id, 4);
	inputFileStream.read( (char*) &header.subchunk1_size, 4);
	inputFileStream.read( (char*) &header.audio_format, 2);
	inputFileStream.read( (char*) &header.num_channels, 2);
	inputFileStream.read( (char*) &header.sample_rate, 4);
	inputFileStream.read( (char*) &header.byte_rate, 4);
	inputFileStream.read( (char*) &header.block_align, 2);
	inputFileStream.read( (char*) &header.bits_per_sample, 2);

	// remove bytes if subchunk1 size is 18
	char extra_bytes[2];
	if(header.subchunk1_size == subchunk1_size_max) {
		inputFileStream.read( extra_bytes, 2);
	}


	// ****** The following code is from TA
	inputFileStream.read( subchunk2_id, 4);
	inputFileStream.read( (char*)&data_size, 4);
	
	*combined_size = data_size / 2;
	int size = data_size / 2;
	file_data = new short[size];

	for(int j = 0 ; j < size; j++) {
		inputFileStream.read((char*) &file_data[j], 2);
	}

	signal = new float[size];
	short temp;
	for(int i = 0; i < size; i++) {
		temp = file_data[i];
		signal[i] = (temp * 1.0) / (pow(2,15) - 1);
		if(signal[i] < -1.0)
			signal[i] = -1.0;

	}
	inputFileStream.close();
	return signal;

	// *********** TA code ends
}

void padding(double paddedArray[], int some_M, float some_x[]) {
	int i;
	for (i = 0; i < (some_M * 2); i+=2) {
		paddedArray[i] = some_x[i/2];
		paddedArray[i+1] = 0;
	}
}

// Function to pad an input array with zeros
double* padArray(float input[], int inputSize, double* paddedArray, int paddedSize) {
    for (int i = 0; i < paddedSize; i += 2) {
        if (i / 2 < inputSize) {
            paddedArray[i] = input[i / 2];
            paddedArray[i + 1] = 0;
        } else {
            paddedArray[i] = 0;
            paddedArray[i + 1] = 0;
        }
    }
	return paddedArray;
}

void pad_zeros_to(double *arr, int current_array_size) {
		for (int i = 0; i < current_array_size; ++i) {
        	arr[current_array_size + i] = 0.0;
    	}
}

void increase_padding(double *arr_output, double *arr_inputfile, double *arr_paddedIR, int arr_size){
	for (int i = 0; i < (arr_size * 2); i+=2) {
		arr_output[i] = (arr_inputfile[i] * arr_paddedIR[i]) - (arr_inputfile[i+1] * arr_paddedIR[i+1]);
		arr_output[i+1] = (arr_inputfile[i+1] * arr_paddedIR[i]) + (arr_inputfile[i] * arr_paddedIR[i+1]);
	}
}

//uses FFT algorithm to convolve
void convolve(float x[], int N, float h[], int M, float y[], int P) {
	int arr_size;

	arr_size = 1;
	
	// Build array of size 2^n
	while (arr_size < P) {
		arr_size *= 2;
	}


	double *inputfile_padded = new double[2 * arr_size];
	double *padded_IR = new double[2 * arr_size];

	for (int i = 0; i < (N * 2); i+=2) {
		inputfile_padded[i] = x[i/2];
		inputfile_padded[i+1] = 0;
	}

	pad_zeros_to(inputfile_padded, arr_size);
	

	for (int i = 0; i < (M * 2); i+=2) {
		padded_IR[i] = h[i/2];
		padded_IR[i+1] = 0;
	}

	pad_zeros_to(padded_IR, arr_size);
	

	double *padded_output = new double[2 * arr_size];
	pad_zeros_to(padded_output, arr_size);


	// call four1 in padded input
	four1((inputfile_padded - 1), arr_size, 1);

	// call four1 on padded impulse
	four1((padded_IR - 1), arr_size, 1);

	// adding paddings
	increase_padding(padded_output, inputfile_padded, padded_IR, arr_size);

	
	four1((padded_output - 1), arr_size, -1);
	

	// Scaling FFT by padding
	for (int k = 0, i = 0; k < arr_size; k++, i+=2) {
		padded_output[i] /= (float)arr_size;
		padded_output[i+1] /= (float)arr_size;
	}
	
	// Removing the padding
	for (int i = 0; i < P; i++) {
		y[i] = padded_output[i*2];
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

// Code is from TA Kimiya Saadat

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




