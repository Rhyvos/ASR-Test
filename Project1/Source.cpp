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
#include "Config.h"
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

	Config *cf = new Config("conf.ini");
	cf->Exist("STATES");
	int states, iterations;
	Audio audio(cf);
	ParamAudio *pm=audio.ParamAudioFile("train1.wav","train1.lab");
	ParamAudio *pm1=audio.ParamAudioFile("train1.wav",NULL);
	ParamAudio *pm2=audio.ParamAudioFile("test1.wav",NULL);

	if(cf->Exist("STATES"))
		states = cf->GetConfig("STATES");
	else
		states=5;

	if(cf->Exist("ITERATIONS"))
		iterations = cf->GetConfig("ITERATIONS");
	else
		iterations=10;


	HMM *hmm_s = new HMM(36,cf,"sil");
	HMM *hmm_a = new HMM(36,cf,"A");
	HMM *hmm_b = new HMM(36,cf,"B"); 
	FinalReEstimate *fre = new FinalReEstimate(cf);
	Recognizer *rec = new Recognizer(cf);



	hmm_s->Initialise(pm,iterations);
	hmm_s->ReEstimate(pm,iterations);
	
	hmm_a->Initialise(pm,iterations);
	hmm_a->ReEstimate(pm,iterations);

	hmm_b->Initialise(pm,iterations);
	hmm_b->ReEstimate(pm,iterations);

	fre->AddHmm(hmm_s);
	fre->AddHmm(hmm_a);
	fre->AddHmm(hmm_b);
	fre->ForwardBackward(pm);
	fre->UpdateModels();
	fre->ForwardBackward(pm);
	fre->UpdateModels();
	fre->ForwardBackward(pm);
	fre->UpdateModels();

	//printfhmm(hmm_s);
	//printfhmm(hmm_a);
	//printfhmm(hmm_b);


	rec->AddHmm(hmm_s);
	rec->AddHmm(hmm_a);
	rec->AddHmm(hmm_b);
	rec->DoRecognition(pm1);
	//rec->DoRecognition(pm1);

	hmm_s->SaveHmm("S");
	hmm_a->SaveHmm("A");
	hmm_b->SaveHmm("B");
	delete pm;
	delete pm1;
	delete pm2;

	delete fre;
	delete rec;
	delete hmm_s;
	delete hmm_a;
	delete hmm_b;
	delete cf;
	
	_CrtDumpMemoryLeaks();
	return 0;
}


