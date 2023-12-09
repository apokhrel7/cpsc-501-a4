
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
// Baseline file in C++


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
	printf("Baseline convolution program finished in %f seconds\n\n", time_elapsed/CLOCKS_PER_SEC);
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


/*
	The convolve() function is from TA Kimiya Saadat
	
    The function convolve takes six arguments: 
        Two input arrays x[] and h[], their respective sizes N and M, and an output array y[] with size P.

    The first loop initializes the output array y[] to zero. 
        This is necessary because the convolution operation involves accumulating values in y[].

    The second loop (outer loop) iterates over each element of the input array x[].

    The third loop (inner loop) iterates over each element of the array h[]. 
        For each pair of elements x[n] and h[m], it adds their sum to the corresponding element in y[].
*/
void convolve(float x[], int N, float h[], int M, float y[], int P){
    int n,m;

	//  Ensure output buffer is the right size: P = N + M - 1  //
    if (P != (N + M - 1)) {
        printf("Output signal vector is the wrong size\n");
        printf("It is %-d, but should be %-d\n", P, (N + M - 1));
        printf("Aborting convolution\n");
        return;
    }

    /* Clear Output Buffer y[] */
    for (n=0; n < P; n++) {
        y[n] = 0.0;
    }

    /* Outer Loop: process each input value x[n] in turn */
    for (n=0; n<N; n++){
        /* Inner loop: process x[n] with each sample of h[n] */
        for (m=0; m<M; m++){
            y[n+m] += x[n] * h[m];
        }
    }
}

