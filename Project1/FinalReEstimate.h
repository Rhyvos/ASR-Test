#pragma once
#include <vector>
#include "HMM.h"
class FinalReEstimate
{
public:
	FinalReEstimate(void);
	~FinalReEstimate(void);
	std::vector<HMM *> hmms;
	HMM ** hmm_seq;							//list of hmms ordered by apperence at label file
	
	float *** outprob;
	double *** alpha;
	double *** beta;
	int frames;
	int seq_num;
	double pruneThresh;
	void AddHmm(HMM *  hmm);
	void LoadHmm(std::string hmm_src);
	void ForwardBackward(ParamAudio * pa);
	int *qHi;
	int *qLo;
	void ListHmms(ParamAudio * pa);
	void SetBeamTaper(int q, int t);
	float *** GetProbability(ParamAudio * pa);
	float OutP(ParamAudio * pa, int fr_number, int segment, int state);
	double Beta(void);
	double Alpha(void);
	double LAdd(double x, double y);
};

