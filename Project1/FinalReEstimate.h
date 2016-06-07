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
	double ** alphat;
	double ** alphat1;
	double *** beta;
	int frames;
	int seq_num;
	double pruneThresh;
	double pr;
	float minFrwdP; 
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
	double MaxModelProb(int q, int t, int minq);
	void StepAlpha(int t, int * start, int * end, int Q, int T, double pr);
	void ZeroAlpha(int qlo, int qhi);
};

