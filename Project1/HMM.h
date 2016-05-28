#pragma once
#include "Structs.h"

#define LZERO  (-1.0E10)   /* ~log(0) */
#define LSMALL (-0.5E10)   /* log values < LSMALL are set to LZERO */
#define MINEARG (-708.3)   /* lowest exp() arg  = log(MINLARG) */
#define MINLARG 2.45E-308  /* lowest log() arg  = exp(MINEARG) */

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
	float minVar;
public:
	HMM(char * hmm_src);
	HMM(int vector_size, int states);
	void Initialise(ParamAudio * pa, int iteration, std::string label_name);
	void GetMean(ParamAudio * pa, std::string label_name);
	void GetVariance(ParamAudio * pa, std::string label_name);
	void FindGConst(void);
	float OutP(ObservationSegment * os, int fr_number, int state_nr);
	float ViterbiAlg(ObservationSegment * os, int * state_vec, int * mixes);

	void ComputeTraceBack(int frames_in_seg, int state, int * states_vec, short ** trace_back);
};

