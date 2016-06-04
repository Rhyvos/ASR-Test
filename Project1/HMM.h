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
	std::string name;
	State * state;
	float ** transition;
	float * obs_count;				//counts occupation of the state int trans
	float ** trans_count;
	int states;
	int vector_size;
	float minVar;
	float epsilon;
	double ** alpha;
	double ** beta;
	double minLogExp;
	int minimum_duration;
public:
	HMM(std::string hmm_src);
	HMM(int vector_size, int states, std::string l_name);
	void Initialise(ParamAudio * pa, int iteration);
	void GetMean(ParamAudio * pa);
	void GetVariance(ParamAudio * pa);
	void FindGConst(void);
	float OutP(ObservationSegment * os, int fr_number, int state_nr);
	float ViterbiAlg(ObservationSegment * os, int * state_vec, int * mixes);

	void ComputeTraceBack(int frames_in_seg, int state, int * states_vec, short ** trace_back);
	void UpdateCounts(ObservationSegment * os, int * state_vec);
	void ResetOldParams(void);
	void UpdateMean(void);
	void UpdateVar(void);
	void UpdateTransition(void);
	void ReEstimate(ParamAudio * pa, int iterations);
	float ** GetProbability(ObservationSegment * os);
	double GetAlpha(float ** prob, int frames);
	double LAdd(double x, double y);
	double GetBeta(float ** prob, int frames);
	void UpdateRestCount(float ** prob,ObservationSegment * os,double pr);
	void SetMinDuration(void);
	void FindSO(int * so, int * d, int s);
};

