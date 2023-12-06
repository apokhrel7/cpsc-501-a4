/**

Unfinished code for a convolution of two input .wav files. The first file being a sample tone and the second being an impulse tone.
A few functions for reading data and writing data have been left out for students to try

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include <stdint.h>

// Function definitions (:
float bytesToFloat(char firstByte, char secondByte) ;
void convolve(float x[], int N, float h[], int M, float y[], int P);
void readTone(char *sampleTone, char *impulseTone);

// struct to hold all data up until the end of subchunk1
// you still have 8 bytes to read before actual data, namely being subchunk2ID & subChunk2Size
// note it may actually be more than 8 bytes as subchunk1size may not be 16!
typedef struct {
    char chunk_id[4];
    int chunk_size;
    char format[4];
    char subchunk1_id[4];
    int subchunk1_size;
    short audio_format;
    short num_channels;
    int sample_rate;
    int byte_rate;
    short block_align;
    short bits_per_sample;
} WavHeader;

// main line of execution
int main (int argc, char *argv[]){
    char *sampleTone = NULL;
    char *impulseTone = NULL;

    /*  Process the command line arguments  */
    if (argc == 3) {
        /*  Set a pointer to the output filename  */
        sampleTone = argv[1]; impulseTone = argv[2];
    }
    else {
        /*  The user did not supply the correct number of command-line
            arguments.  Print out a usage message and abort the program.  */
        fprintf(stderr, "Usage:  %s sampleTone impulseTone\n", argv[0]); exit(-1);
    }
    readTone(sampleTone, impulseTone);
}

/**
Read the tones, and call convolve on them
*/
void readTone(char *sampleTone, char *impulseTone){
    FILE *sampleFileStream = fopen(sampleTone, "rb");
    FILE *impulseFileStream = fopen(impulseTone, "rb");
    FILE *outputFileStream = fopen("output.wav", "wb");

    WavHeader header_sample;
    WavHeader header_impulse;
    // read the header subchunk 1, write the header into a new file
    fread(&header_sample, sizeof(header_sample), 1, sampleFileStream);
    fread(&header_impulse, sizeof(header_impulse), 1, impulseFileStream);

    if (header_sample.subchunk1_size != 16){
        // eliminate Null Bytes
        int remainder = header_sample.subchunk1_size -16;
        char randomVar[remainder];
        fread(randomVar, remainder, 1, sampleFileStream);
    }
    
    if (header_sample.subchunk1_size != 16){
        // eliminate Null Bytes
        int remainder = header_impulse.subchunk1_size -16;
        char randomVar[remainder];
        fread(randomVar, remainder, 1, impulseFileStream);
    }
    char subchunk2_id_sample[4];
    char subchunk2_id_impulse[4];
    int subchunk2_size_sample; // an integer is 4 bytes
    int subchunk2_size_impulse; // an integer is 4 bytes
    fread(&subchunk2_id_sample, sizeof(subchunk2_id_sample), 1, sampleFileStream);
    fread(&subchunk2_size_sample, sizeof(subchunk2_size_sample), 1, sampleFileStream);
    fread(&subchunk2_id_impulse, sizeof(subchunk2_id_impulse), 1, impulseFileStream);
    fread(&subchunk2_size_impulse, sizeof(subchunk2_size_impulse), 1, impulseFileStream);

    int num_samples = subchunk2_size_sample / (header_sample.bits_per_sample / 8); // number of data points in the sample
    int num_impulse = subchunk2_size_impulse / (header_impulse.bits_per_sample / 8); // number of data points in the impulse
    
    /*
            Now Please Try the following:
            1. Read the two files, extract from each of them a float array, and a size. 
                You may find it helpful to use the other code provided in this tutorial for file reading.
            2. Call Convolve on the outputs
            3. Write the output.wav, making sure to write the header information as well.
    */
}


// Function to convert two bytes to one float in the range -1 to 1
// This is used as .wav files store data in short format (typically 16 bits, can also be extracted from the bits_per_sample header)
// assumes the data is read in bytes
float bytesToFloat(char firstByte, char secondByte) {
    // Convert two bytes to one short (little endian)
    short s = (secondByte << 8) | firstByte;
    // Convert to range from -1 to (just below) 1
    return s / 32768.0;
}

/*
    The function convolve takes six arguments: 
        Two input arrays x[] and h[], their respective sizes N and M, and an output array y[] with size P.

    The first loop initializes the output array y[] to zero. 
        This is necessary because the convolution operation involves accumulating values in y[].

    The second loop (outer loop) iterates over each element of the input array x[].

    The third loop (inner loop) iterates over each element of the array h[]. 
        For each pair of elements x[n] and h[m], it adds their sum to the corresponding element in y[].
*/
void convolve(float x[], int N, float h[], int M, float y[], int P)
{
    int n,m;

    /* Clear Output Buffer y[] */
    for (n=0; n < P; n++)
    {
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