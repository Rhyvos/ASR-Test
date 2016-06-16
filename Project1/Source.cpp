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
#include "Audio.h"
#include "HMM.h"
#include "FinalReEstimate.h"
#include "Recognizer.h"
#include "Config.h"
#include <fstream>
//#define _CRT_SECURE_NO_DEPRECATE
#pragma warning (disable : 4996)
using namespace std;


std::vector<std::string> splitstring(char * str,const char * delimiter) {
  std::vector<std::string> internal;
  char * buff;
  buff = strtok(str,delimiter);
  while(buff != NULL) {
    internal.push_back(buff);
	buff = strtok(NULL,delimiter);
  }
  
  return internal;
}


int main(int argc, char * argv[])
{
	Config *cf = NULL;
	
	vector<ParamAudio *> pa;
	vector<HMM *> hmms;
	vector<string> tmp;
	vector<string> tmp1;
	int iterations = 10;
	for (int i = 1; i < argc; i++)
	{
		if(string(argv[i]).compare("-c") == 0)
		{
			i++;
			cf = new Config(argv[i]);
			if(cf->Exist("ITERATIONS"))
				iterations = cf->GetConfig("ITERATIONS");
		}


		if(string(argv[i]).compare("-s") == 0)
		{
			i++;
			Audio *audio = new Audio(cf);
			ifstream input(argv[i],ifstream::in);
			char buffer[256];
			while (input.good()) 
			{
				input.getline(buffer,256);
				tmp = splitstring(buffer," \t");
				if(tmp.size() == 2)
					pa.push_back(audio->ParamAudioFile(tmp[0].c_str(),tmp[1].c_str()));
				else if(tmp.size() == 1)
					pa.push_back(audio->ParamAudioFile(tmp[0].c_str(),""));
				else
					fprintf(stderr,"wrong script line: %s\n",string(buffer).c_str());
			}
			input.close();
			delete audio;
		}

		if(string(argv[i]).compare("-h") == 0)
		{
			i++;
			ifstream input(argv[i],ifstream::in);
			char buffer[256];
			while (input.good()) 
			{
				input.getline(buffer,256);
				tmp = splitstring(buffer," \t");
				if(tmp.size() == 1)
					hmms.push_back(new HMM(tmp[0].c_str(),cf));
				else
					fprintf(stderr,"wrong hmm list line: %s\n",string(buffer).c_str());
			}
			input.close();
		}


		if(string(argv[i]).compare("-g") == 0)
		{
			i++;
			ifstream input(argv[i],ifstream::in);
			char buffer[256];
			while (input.good()) 
			{
				input.getline(buffer,256);
				tmp = splitstring(buffer," \t");
				

				if(tmp.size() == 1)
				{
					string s;
					size_t st,en;
					s = tmp[0];
					st = s.find_last_of("/");
					en = s.find_last_of(".");
					if(st == string::npos)
						st = 0;
					if( en == string::npos)
						en = s.size();
					s = s.substr(st+1,en-st-1);
					hmms.push_back(new HMM(cf,s.c_str())); 
				}
				else
					fprintf(stderr,"wrong hmm list line: %s\n",string(buffer).c_str());
			}
			input.close();
		}

		if(string(argv[i]).compare("-i") == 0)
		{
			for (int j = 0; j < pa.size(); j++)
			{
				for (int k = 0; k < hmms.size(); k++)
				{
					hmms[k]->Initialise(pa[j],iterations);
				}
			}
		}
		else if(string(argv[i]).compare("-t") == 0)
		{
			for (int j = 0; j < pa.size(); j++)
			{
				for (int k = 0; k < hmms.size(); k++)
				{
					hmms[k]->ReEstimate(pa[j],iterations);
				}
			}

		}
		else if(string(argv[i]).compare("-f") == 0)
		{
			FinalReEstimate *fr = new FinalReEstimate(cf);
			for (int k = 0; k < hmms.size(); k++)
				fr->AddHmm(hmms[k]);

			for (int j = 0; j < pa.size(); j++)
				fr->ForwardBackward(pa[j]);

			delete fr;
		}
		else if(string(argv[i]).compare("-r") == 0)
		{
			Recognizer *rec = new Recognizer(cf);
			for (int k = 0; k < hmms.size(); k++)
				rec->AddHmm(hmms[k]);

			for (int j = 0; j < pa.size(); j++)
				rec->DoRecognition(pa[j]);
			delete rec;
		}

	}
	

	for (int i = 0; i < pa.size(); i++)
	{
		delete pa[i];
	}

	for (int i = 0; i < hmms.size(); i++)
	{
		hmms[i]->SaveHmm();
		delete hmms[i];
	}

	if(cf != NULL)
	{
		delete cf;
	}
	tmp.clear();
	pa.clear();
	hmms.clear();

	
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();  
#endif // _DEBUG

	return 0;
	
}


