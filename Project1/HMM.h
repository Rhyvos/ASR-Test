#pragma once
#include "Structs.h"
class HMM
{
public:
	~HMM(void);
//private:
	HMM(void);
	State * state;
	float ** transition;
	int states;
	int vector_size;
public:
	HMM(char * hmm_src);
	HMM(int vector_size, int states);
	void Initialise(ParamAudio * pa, int iteration, std::string label_name);
	void GetMean(ParamAudio * pa, std::string label_name);
	void GetVariance(ParamAudio * pa, std::string label_name);
};

