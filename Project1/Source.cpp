#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
#define _USE_MATH_DEFINES
#include <iostream>
#include "MFCC.h"
#include "FFT.h"
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <vector>
#include "Audio.h"
#include "hmm.h"
#include "FinalReEstimate.h"
#include "Recognizer.h"
//#define _CRT_SECURE_NO_DEPRECATE
#pragma warning (disable : 4996)
using namespace std;

void printfhmm(HMM *hmm)
{
	for (int i = 0; i < hmm->states-2; i++)
	{
		printf("\nState %d:\n",hmm->state[i].state_nr);
		printf("Mean\n");
		for (int j = 0; j < hmm->vector_size; j++)
		{
			printf("%f\t",hmm->state[i].mean[j]);

		}
		printf("\nVariance\n");
		for (int j = 0; j < hmm->vector_size; j++)
		{
			printf("%f\t",hmm->state[i].var[j]);

		}
		printf("\nGConst\t%f",hmm->state[i].g_const);

	}
	printf("\n");
	for (int i = 0; i <hmm->states; i++)
	{
		for (int j = 0; j < hmm->states; j++)
		{
			printf("%f ",exp(hmm->transition[i][j]));
		}
		printf("\n");
	}
}



int main()
{
	
	
	

	Audio audio;
	ParamAudio *pm=audio.ParamAudioFile("train1.wav","train1.lab");
	ParamAudio *pm1=audio.ParamAudioFile("train1.wav",NULL);
	ParamAudio *pm2=audio.ParamAudioFile("test1.wav",NULL);

	HMM *hmm_s = new HMM(36,5,"sil");
	HMM *hmm_a = new HMM(36,5,"A");
	HMM *hmm_b = new HMM(36,5,"B");
	FinalReEstimate *fre = new FinalReEstimate();
	Recognizer *rec = new Recognizer();


	hmm_s->Initialise(pm,10);
	hmm_s->minVar = 0.05;
	hmm_s->ReEstimate(pm,10);
	
	hmm_a->Initialise(pm,10);
	hmm_a->minVar = 0.05;
	hmm_a->ReEstimate(pm,10);

	hmm_b->Initialise(pm,10);
	hmm_b->minVar = 0.05;
	hmm_b->ReEstimate(pm,10);

	fre->AddHmm(hmm_s);
	fre->AddHmm(hmm_a);
	fre->AddHmm(hmm_b);
	fre->ForwardBackward(pm);
	fre->UpdateModels();
	fre->ForwardBackward(pm);
	fre->UpdateModels();
	fre->ForwardBackward(pm);
	fre->UpdateModels();

	printfhmm(hmm_s);
	printfhmm(hmm_a);
	printfhmm(hmm_b);


	rec->AddHmm(hmm_s);
	rec->AddHmm(hmm_a);
	rec->AddHmm(hmm_b);
	rec->DoRecognition(pm1);

	hmm_s->SaveHmm("S");
	hmm_a->SaveHmm("A");
	hmm_b->SaveHmm("B");
	delete pm;
	delete pm1;
	delete hmm_s;
	delete hmm_a;
	delete hmm_b;
	delete fre;
	
	
	/*
	for(int i=0; i<300; i++)
			{
				printf("buffer_array[%d] = buffer_array[%d];\n",i,400 - 300 + i);
			}
	
	ofstream myFile;
    myFile.open("mizi.txt");
	MFCC * test = new MFCC();
	FFT * fft = new FFT();
	float data_in[512],data_in1[513],real_out[512],img_out[512], data_out[13], processed[256]; 
	srand(time(NULL));

	FILE * infile = fopen("a1.wav","rb");		
	ofstream output;	
	output.open("invi1_spectrum11.txt", std::ofstream::out | std::ofstream::app);
	int count = 0, windows=0;
	short buffer;
	short int tmp;
	float tmp1;
	unsigned short int tmp2;
	float last_sample = 1;
	float preemphasis_param = 0.97;
	int nb;			
	const float ONEOVERSHORTMAX = 3.0517578125e-5f; 
	for(int i=0;i<256;i++)
     processed[i]=0;
	data_in1[0]=0.0;
	if (infile)
	{

		int i=0;
		nb = fread(&HeaderInfo,1,sizeof(HeaderInfo),infile);
		
		printf("%d\n",(int)floor(HeaderInfo.Subchunk2Size/240));
		printf("%d",(int)floor(HeaderInfo.Subchunk2Size/240)*240);
		float** mfcc_data = new float*[HeaderInfo.Subchunk2Size/160];
		for(int i = 0; i < HeaderInfo.Subchunk2Size/160; ++i)
			mfcc_data[i] = new float[12];

		float** deltas = new float*[HeaderInfo.Subchunk2Size/160];
		for(int i = 0; i < HeaderInfo.Subchunk2Size/160; ++i)
			deltas[i] = new float[12];

		float** acceleration = new float*[HeaderInfo.Subchunk2Size/160];
		for(int i = 0; i < HeaderInfo.Subchunk2Size/160; ++i)
			acceleration[i] = new float[12];

		while ((nb = fread(&buffer,1,sizeof(buffer),infile))>0)
		{
			  
				data_in[i++] = (float)buffer;
			    data_in1[i]   = (float)buffer;

			   if (i>=400)
			   {
				   
				   
				  
				   
				   test->Compute(i,12,data_in,mfcc_data[windows]);
				   
				   
				  
				
				   windows++;
				   for(int i=0;i<240;i++)
							data_in[i]=data_in[i+160];
					   i=240;
					
			   }

			
		}
		test -> AddRegression(mfcc_data,deltas,windows,12,2);
		test -> AddRegression(deltas,acceleration,windows,12,2);
		for (int j = 0; j < windows; j++)
		{
			for (int i = 0; i < 12; i++)
					{
						printf("MFCC[%d][%d] = %f\n",j,i,mfcc_data[j][i]);	
					}
			for (int i = 0; i < 12; i++)
					{
									printf("Deltas[%d][%d] = %f\n",j,i,deltas[j][i]);	
									
					}
			for (int i = 0; i < 12; i++)
					{
									printf("Acceleration[%d][%d] = %f\n",j,i,acceleration[j][i]);	
									
					}
			getchar();
		}
		
		

	}
	
	

	output.close();
	
	*/
	
	//getchar();
	_CrtDumpMemoryLeaks();
	return 0;
}


